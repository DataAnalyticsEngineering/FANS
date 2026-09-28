// libfans_neml2.so: the NEML2 material plugin (fans_plugin.h).
//
// Mandel order: FANS is [11,22,33,12,13,23], NEML2 [11,22,33,23,13,12].
// History: every input "X~1" with a matching output "X"; SR2 ones are stored
// in FANS order too.

#include "fans_plugin.h"

#include "neml2/csrc/aoti/Model.h"
#include <ATen/ATen.h>
#include <ATen/Parallel.h>
#include <c10/util/accumulate.h>

#include <algorithm>
#include <cstdio>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct HistoryVar {
    std::string          name;  // output name; the input is name + "~1"
    std::vector<int64_t> shape; // {-1, per-point shape...}
    int64_t              size, offset;
    bool                 sr2;
};

const at::TensorOptions    kF64 = at::TensorOptions().dtype(at::kDouble);
const std::vector<int64_t> sr2{6}; // per-point shape of a symmetric tensor

} // namespace

struct FANSPluginModel {
    std::unique_ptr<neml2::aoti::Model> model;
    std::vector<HistoryVar>             history;
    std::string                         strain_name, stress_name;
    int64_t                             n_state{0};
    bool                                orientation{false};
    at::Device                          device{at::kCPU};
    at::Tensor                          perm = at::tensor({0, 1, 2, 5, 4, 3}, at::kLong);
};

extern "C" {

FANSPluginModel *fans_plugin_load(const char *spec, const char *device, int *n_state,
                                  int *wants_orientation, char *msg, size_t msglen)
{
    try {
        // FANS parallelises with MPI; one torch thread per rank
        static const bool threads_set = [] {
            at::set_num_threads(1);
            at::set_num_interop_threads(1);
            return true;
        }();
        (void) threads_set;

        auto h    = std::make_unique<FANSPluginModel>();
        h->device = at::Device(device);
        h->model  = std::make_unique<neml2::aoti::Model>(spec, h->device, at::kDouble);

        const auto &in = h->model->input_names(), &out = h->model->output_names();
        const auto &in_shape = h->model->input_base_shapes(), &out_shape = h->model->output_base_shapes();
        auto        find = [](const auto &v, const std::string &n) {
            const auto it = std::find(v.begin(), v.end(), n);
            return it == v.end() ? -1 : int(it - v.begin());
        };

        // 'strain'/'stress' for a single model, 'forces/E'/'state/S' for a composed one
        const int si = std::max(find(in, "strain"), find(in, "forces/E"));
        const int oi = std::max(find(out, "stress"), find(out, "state/S"));
        if (si < 0 || in_shape[si] != sr2)
            throw std::runtime_error("artifact needs an SR2 strain input named 'strain' or 'forces/E'");
        if (oi < 0 || out_shape[oi] != sr2)
            throw std::runtime_error("artifact needs an SR2 stress output named 'stress' or 'state/S'");
        h->strain_name = in[si];
        h->stress_name = out[oi];

        for (size_t i = 0; i < in.size(); ++i) {
            if (int(i) == si)
                continue;
            const std::string &name = in[i];
            if (name == "orientation" && in_shape[i] == std::vector<int64_t>{3, 3}) {
                h->orientation = true;
                continue;
            }
            const int bi = name.ends_with("~1") ? find(out, name.substr(0, name.size() - 2)) : -1;
            if (bi < 0 || out_shape[bi] != in_shape[i])
                throw std::runtime_error("artifact input '" + name + "' is neither the strain, an R2 'orientation', nor " +
                                         "a history variable 'X~1' with a matching output 'X'; FANS has no source for it");
            std::vector<int64_t> shape{-1};
            shape.insert(shape.end(), in_shape[i].begin(), in_shape[i].end());
            const int64_t size = c10::multiply_integers(in_shape[i]);
            h->history.push_back({out[bi], shape, size, h->n_state, in_shape[i] == sr2});
            h->n_state += size;
        }
        *n_state           = int(h->n_state);
        *wants_orientation = h->orientation;

        std::string vars;
        for (const HistoryVar &v : h->history)
            vars += (vars.empty() ? "[\"" : ", [\"") + v.name + "\", " + std::to_string(v.size) + "]";
        vars = "[" + vars + "]";
        if (vars.size() >= msglen)
            throw std::runtime_error("too many history variables to describe");
        std::snprintf(msg, msglen, "%s", vars.c_str());
        return h.release();
    } catch (const std::exception &e) {
        std::snprintf(msg, msglen, "%s", e.what());
        return nullptr;
    }
}

int fans_plugin_evaluate(FANSPluginModel *m, size_t n_points, const double *strain,
                         const double *orientation, const double *state_old, double *stress,
                         double *state_new, char *err, size_t errlen)
{
    try {
        const int64_t                     n = int64_t(n_points);
        std::map<std::string, at::Tensor> inputs;
        inputs[m->strain_name] = at::from_blob(const_cast<double *>(strain), {n, 6}, kF64).index_select(1, m->perm).to(m->device);
        if (m->orientation)
            inputs["orientation"] = at::from_blob(const_cast<double *>(orientation), {n, 3, 3}, kF64).to(m->device);

        const at::Tensor old_state = at::from_blob(const_cast<double *>(state_old), {n, m->n_state}, kF64);
        const at::Tensor new_state = at::from_blob(state_new, {n, m->n_state}, kF64);
        for (const HistoryVar &h : m->history) {
            at::Tensor v          = old_state.narrow(1, h.offset, h.size);
            inputs[h.name + "~1"] = (h.sr2 ? v.index_select(1, m->perm) : v).reshape(h.shape).to(m->device);
        }

        const auto outputs = m->model->forward(inputs);

        at::Tensor stress_out = at::from_blob(stress, {n, 6}, kF64);
        at::index_select_out(stress_out, outputs.at(m->stress_name).cpu(), 1, m->perm);
        for (const HistoryVar &h : m->history) {
            at::Tensor v = outputs.at(h.name).reshape({n, h.size}).cpu();
            new_state.narrow(1, h.offset, h.size).copy_(h.sr2 ? v.index_select(1, m->perm) : v);
        }
        return 0;
    } catch (const std::exception &e) {
        std::snprintf(err, errlen, "%s", e.what());
        return 1;
    }
}

void fans_plugin_free(FANSPluginModel *m)
{
    delete m;
}

} // extern "C"
