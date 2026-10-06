// libfans_neml2.so: the NEML2 material plugin (fans_plugin.h).
//
// Material properties: "artifact" (from neml2-compile), "gradient" and "flux",
// the names of the model's gradient input and flux output (e.g. 'forces/E' and
// 'state/S' in small strain, 'forces/F' and 'state/P' in large strain),
// "device" ("cpu" by default, "cuda", ...) and "batch_size" (points per
// evaluation, by default 1024 on the CPU and 65536 on a GPU).
// The model's other inputs are
//   time     't' and 't~1' of time-integrated models
//   history  every 'X~1' with a matching output 'X', plus its Newton initial
//            guess 'X' if the model has no predictor
//   fields   every other input, e.g. an R2 'orientation'
// Every axis of length 6 is a Mandel index (an SR2's, both of an SSR4 such as a
// stiffness), ordered [11,22,33,23,13,12] in NEML2 and [11,22,33,12,13,23] in
// FANS; everything else has the same layout in both.

#include "fans_plugin.h"

#include "neml2/csrc/aoti/Model.h"
#include <ATen/ATen.h>
#include <ATen/Parallel.h>
#include <c10/util/accumulate.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstdio>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// A model variable FANS passes as `size` doubles per point
struct Var {
    std::string          name;             // history: the output, the input being name + "~1"
    std::vector<int64_t> shape;            // {-1, per-point shape...}
    int64_t              size;             // doubles per point
    int64_t              offset{0};        // history: position among the history doubles
    at::Tensor           table;            // field: its rows, on the device
    bool                 per_phase{false}; // field: a row per phase, else per voxel
};

Var make_var(const std::string &name, const std::vector<int64_t> &base_shape)
{
    std::vector<int64_t> shape{-1};
    shape.insert(shape.end(), base_shape.begin(), base_shape.end());
    return {name, shape, c10::multiply_integers(base_shape)};
}

const at::TensorOptions kF64 = at::TensorOptions().dtype(at::kDouble);

} // namespace

struct FANSPluginModel {
    std::unique_ptr<neml2::aoti::Model> model;
    at::Device                          device{at::kCPU};
    Var                                 gradient, flux;
    Var                                 tangent; // d flux / d gradient: the flux's shape, then the gradient's
    std::vector<Var>                    history, fields;
    int64_t                             n_history{0};                                     // history doubles per point
    at::Tensor                          perm = at::tensor({0, 1, 2, 5, 4, 3}, at::kLong); // FANS <-> NEML2 Mandel order

    // v in the variable's shape, every Mandel axis swapped between FANS and NEML2
    // order (both ways, the swap undoes itself)
    at::Tensor reorder(const at::Tensor &v, const Var &var) const
    {
        at::Tensor t = v.reshape(var.shape);
        for (size_t d = 1; d < var.shape.size(); ++d)
            if (var.shape[d] == 6)
                t = t.index_select(int64_t(d), perm);
        return t;
    }
    // FANS -> NEML2: [n][size] doubles to the variable, on the device
    at::Tensor to_neml2(const at::Tensor &v, const Var &var) const
    {
        return reorder(v, var).to(device);
    }
    // NEML2 -> FANS: the variable into `out`, [n][size] doubles on the host
    void to_fans(const at::Tensor &v, const Var &var, at::Tensor out) const
    {
        out.copy_(reorder(v.cpu(), var).reshape(out.sizes()));
    }
};

extern "C" {

FANSPluginModel *fans_plugin_load(const char *config, char *msg, size_t msglen)
{
    try {
        // FANS parallelises with MPI; one torch thread per rank, set once per process
        static std::once_flag threads_set;
        std::call_once(threads_set, [] {
            at::set_num_threads(1);
            at::set_num_interop_threads(1);
        });

        const auto props = nlohmann::json::parse(config);
        auto       h     = std::make_unique<FANSPluginModel>();
        h->device        = at::Device(props.value("device", "cpu"));
        h->model         = std::make_unique<neml2::aoti::Model>(props.at("artifact").get<std::string>(), h->device, at::kDouble);

        const auto &in = h->model->input_names(), &out = h->model->output_names();
        const auto &in_shape = h->model->input_base_shapes(), &out_shape = h->model->output_base_shapes();
        auto        find = [](const auto &v, const std::string &n) {
            const auto it = std::find(v.begin(), v.end(), n);
            return it == v.end() ? -1 : int(it - v.begin());
        };

        if (!props.contains("gradient") || !props.contains("flux"))
            throw std::runtime_error("the material properties \"gradient\" and \"flux\" must name the model's gradient input and "
                                     "flux output, e.g. 'forces/E' and 'state/S' in small strain");
        const std::string gradient = props["gradient"], flux = props["flux"];
        const int         gi = find(in, gradient), fi = find(out, flux);
        if (gi < 0)
            throw std::runtime_error("the artifact has no input '" + gradient + "' (the \"gradient\")");
        if (fi < 0)
            throw std::runtime_error("the artifact has no output '" + flux + "' (the \"flux\")");
        h->gradient                        = make_var(gradient, in_shape[gi]);
        h->flux                            = make_var(flux, out_shape[fi]);
        std::vector<int64_t> tangent_shape = out_shape[fi];
        tangent_shape.insert(tangent_shape.end(), in_shape[gi].begin(), in_shape[gi].end());
        h->tangent = make_var("tangent", tangent_shape);
        if (h->flux.size != h->gradient.size)
            throw std::runtime_error("the gradient '" + gradient + "' and the flux '" + flux + "' differ in size");

        for (size_t i = 0; i < in.size(); ++i) {
            const std::string &name = in[i];
            if (name == gradient || name == "t" || name == "t~1")
                continue;
            if (find(in, name + "~1") >= 0)
                continue; // initial guess, fed the latest trial history
            Var       var = make_var(name, in_shape[i]);
            const int bi  = name.ends_with("~1") ? find(out, name.substr(0, name.size() - 2)) : -1;
            if (bi >= 0 && out_shape[bi] == in_shape[i]) { // history
                var.name   = out[bi];
                var.offset = h->n_history;
                h->n_history += var.size;
                h->history.push_back(var);
            } else {
                h->fields.push_back(var);
            }
        }

        nlohmann::json description = {{"gradient", h->gradient.size},
                                      {"batch_size", props.value("batch_size", h->device.is_cpu() ? 1024 : 65536)},
                                      {"history", nlohmann::json::array()},
                                      {"fields", nlohmann::json::array()}};
        for (const Var &v : h->history)
            description["history"].push_back(nlohmann::json::array({v.name, v.size}));
        for (const Var &v : h->fields)
            description["fields"].push_back(nlohmann::json::array({v.name, v.size}));
        const std::string text = description.dump();
        if (text.size() >= msglen)
            throw std::runtime_error("too many history variables and fields to describe");
        std::snprintf(msg, msglen, "%s", text.c_str());
        return h.release();
    } catch (const std::exception &e) {
        std::snprintf(msg, msglen, "%s", e.what());
        return nullptr;
    }
}

int fans_plugin_set_table(FANSPluginModel *m, size_t field, int per_phase, size_t n_rows, const double *table, char *err,
                          size_t errlen)
{
    try {
        Var             &var  = m->fields.at(field);
        const at::Tensor rows = at::from_blob(const_cast<double *>(table), {int64_t(n_rows), var.size}, kF64);
        var.table             = m->to_neml2(rows, var).clone(); // our own copy: FANS may free `table`
        var.per_phase         = per_phase;
        return 0;
    } catch (const std::exception &e) {
        std::snprintf(err, errlen, "%s", e.what());
        return 1;
    }
}

int fans_plugin_evaluate(FANSPluginModel *m, size_t n_points, double t_old, double t, const double *gradient,
                         const int *voxel, const int *phase, const double *history_old, double *flux,
                         double *history_new, double *tangent, char *err, size_t errlen)
{
    try {
        // A FANS array as an [n][size] tensor, without copying
        const int64_t n    = int64_t(n_points);
        auto          view = [n](const double *p, int64_t size) { return at::from_blob(const_cast<double *>(p), {n, size}, kF64); };

        std::map<std::string, at::Tensor> inputs;
        inputs[m->gradient.name] = m->to_neml2(view(gradient, m->gradient.size), m->gradient);
        inputs["t"]              = at::full({n}, t, kF64.device(m->device));
        inputs["t~1"]            = at::full({n}, t_old, kF64.device(m->device));
        for (const Var &f : m->fields) { // each point's row of the table
            const auto rows = at::from_blob(const_cast<int *>(f.per_phase ? phase : voxel), {n}, at::kInt).to(m->device, at::kLong);
            inputs[f.name]  = f.table.index_select(0, rows);
        }

        const at::Tensor old_history = view(history_old, m->n_history), new_history = view(history_new, m->n_history);
        for (const Var &h : m->history) {
            inputs[h.name + "~1"] = m->to_neml2(old_history.narrow(1, h.offset, h.size), h);
            inputs[h.name]        = m->to_neml2(new_history.narrow(1, h.offset, h.size), h); // initial guess
        }

        std::map<std::string, at::Tensor> outputs;
        if (tangent) { // the same evaluation, with its Jacobian
            neml2::aoti::VariablePairJacobian jacobian;
            std::tie(outputs, jacobian) = m->model->jacobian(inputs);
            m->to_fans(jacobian.at(m->flux.name).at(m->gradient.name), m->tangent, view(tangent, m->tangent.size));
        } else {
            outputs = m->model->forward(inputs);
        }

        m->to_fans(outputs.at(m->flux.name), m->flux, view(flux, m->flux.size));
        for (const Var &h : m->history)
            m->to_fans(outputs.at(h.name), h, new_history.narrow(1, h.offset, h.size));
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
