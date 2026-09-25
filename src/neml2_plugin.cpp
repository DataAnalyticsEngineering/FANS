// libfans_neml2.so: the NEML2 material plugin for FANS, implementing
// fans_plugin.h. Built only with -DFANS_NEML2=ON, as its own library,
// so FANS never links libtorch.
//
// Two translations happen here:
//
//   Component order. Both sides use sqrt(2)-scaled Mandel storage but
//   disagree on which slot holds the 12 and 23 shears -- FANS is
//   [11,22,33,12,13,23], NEML2 [11,22,33,23,13,12] -- so the mapping swaps
//   slots 3 and 5, and is its own inverse.
//
//   History. NEML2 names the previous-step value of a state variable by
//   appending "~1", so every input "X~1" with a matching output "X" is one
//   history variable, packed into the opaque per-point block FANS carries.
//
// An R2 input named "orientation" receives FANS's per-grain crystal-to-sample
// rotation matrices, e.g. for the orientation_matrix of crystal plasticity.

#include "fans_plugin.h"

#include "neml2/csrc/aoti/Model.h"
#include <ATen/ATen.h>
#include <ATen/Parallel.h>

#include <algorithm>
#include <cstdio>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct HistoryVar {
    std::string          old_name, new_name;
    std::vector<int64_t> shape; //!< per-point shape
    int64_t              size, offset;
};

const at::TensorOptions kF64 = at::TensorOptions().dtype(at::kDouble);

} // namespace

struct FANSPluginModel {
    std::unique_ptr<neml2::aoti::Model> model;
    std::vector<HistoryVar>             history;
    std::string                         strain_name, stress_name;
    int64_t                             n_state{0};
    bool                                orientation{false};
    at::Tensor                          perm = at::tensor({0, 1, 2, 5, 4, 3}, at::kLong);
};

extern "C" {

FANSPluginModel *fans_plugin_load(const char *spec, int *n_state, int *wants_orientation,
                                  char *err, size_t errlen)
{
    try {
        // libtorch would start a thread per core in every MPI rank; FANS
        // already parallelises by domain decomposition.
        static const bool threads_set = [] {
            at::set_num_threads(1);
            at::set_num_interop_threads(1);
            return true;
        }();
        (void) threads_set;

        auto h   = std::make_unique<FANSPluginModel>();
        h->model = std::make_unique<neml2::aoti::Model>(spec, at::kCPU, at::kDouble);

        const auto &in = h->model->input_names(), &out = h->model->output_names();
        const auto &in_shape = h->model->input_base_shapes(), &out_shape = h->model->output_base_shapes();
        auto        find = [](const auto &v, const std::string &n) {
            const auto it = std::find(v.begin(), v.end(), n);
            return it == v.end() ? -1 : int(it - v.begin());
        };

        // A single-model artifact names these 'strain'/'stress', a composed
        // one 'forces/E'/'state/S'. (neml2-compile's --rename-input/-output
        // are a no-op in 3.1.0, so they cannot be normalised there.)
        const std::vector<int64_t> sr2{6};
        const int                  si = std::max(find(in, "strain"), find(in, "forces/E"));
        const int                  oi = std::max(find(out, "stress"), find(out, "state/S"));
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
            const bool        old  = name.size() > 2 && name.compare(name.size() - 2, 2, "~1") == 0;
            const std::string base = old ? name.substr(0, name.size() - 2) : "";
            const int         bi   = old ? find(out, base) : -1;
            if (bi < 0 || out_shape[bi] != in_shape[i])
                throw std::runtime_error("artifact input '" + name + "' is neither the strain, an R2 'orientation', nor " +
                                         "a history variable 'X~1' with a matching output 'X'; FANS has no source for it");
            int64_t size = 1;
            for (int64_t d : in_shape[i])
                size *= d;
            h->history.push_back({name, base, in_shape[i], size, h->n_state});
            h->n_state += size;
        }
        *n_state           = int(h->n_state);
        *wants_orientation = h->orientation;
        return h.release();
    } catch (const std::exception &e) {
        std::snprintf(err, errlen, "%s", e.what());
        return nullptr;
    }
}

int fans_plugin_evaluate(FANSPluginModel *m, size_t n_points, const double *strain,
                         const double *orientation, const double *state_old, double *stress,
                         double *state_new, char *err, size_t errlen)
{
    try {
        const int64_t n = int64_t(n_points);
        // Views on FANS's buffers; index_select both reorders and copies.
        std::map<std::string, at::Tensor> inputs;
        inputs[m->strain_name] = at::from_blob(const_cast<double *>(strain), {n, 6}, kF64).index_select(1, m->perm);
        if (m->orientation)
            inputs["orientation"] = at::from_blob(const_cast<double *>(orientation), {n, 3, 3}, kF64);

        at::Tensor old_state, new_state;
        if (m->n_state) {
            old_state = at::from_blob(const_cast<double *>(state_old), {n, m->n_state}, kF64);
            new_state = at::from_blob(state_new, {n, m->n_state}, kF64);
        }
        for (const HistoryVar &h : m->history) {
            std::vector<int64_t> shape{n};
            shape.insert(shape.end(), h.shape.begin(), h.shape.end());
            inputs[h.old_name] = old_state.narrow(1, h.offset, h.size).reshape(shape);
        }

        // A local Newton solve that fails throws here, losing the whole
        // evaluation: NEML2 does not report which points failed.
        const auto outputs = m->model->forward(inputs);

        at::Tensor stress_out = at::from_blob(stress, {n, 6}, kF64);
        at::index_select_out(stress_out, outputs.at(m->stress_name), 1, m->perm);
        for (const HistoryVar &h : m->history)
            new_state.narrow(1, h.offset, h.size).copy_(outputs.at(h.new_name).reshape({n, h.size}));
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
