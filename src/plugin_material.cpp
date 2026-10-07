// A material evaluated by a plugin loaded at run time (fans_plugin.h), e.g.
// NEML2 in libfans_neml2.so, in thermal, small- or large-strain problems, and
// the solver's batched evaluation for it. Kept out of the headers so the
// native models' machine code is unaffected.
//
// Batching: before each sweep over the elements, evaluate_batch gets the flux
//   (the stress) at all the model's Gauss points from the plugin in large
//   batches (fast on a GPU); the sweep then takes these cached values.
// History: the plugin's internal variables at every Gauss point: `converged`,
//   the state a step starts from, and `trial`, the state of the step itself.
//   Once a step has converged, the solver evaluates it a last time (for the
//   output) and `trial` becomes `converged`; `previous` keeps the state that
//   step started from, for its tangent.
// Fields: further model inputs, e.g. an orientation, read from the
//   microstructure file per voxel or per phase (read_field).
// Linear: "linear": true in the material properties says that the flux is
//   linear in the gradient; if all materials are linear, FANS solves with its
//   linear CG, without a line search. FANS does not check it.

#include "general.h"
#include "matmodel.h"
#include "LargeStrainMechModel.h"
#include "solver.h"
#include "fans_plugin.h"

#include <algorithm>
#include <cctype>
#include <dlfcn.h>

namespace {

// The functions of fans_plugin.h, found in the plugin's library at run time
struct Plugin {
    decltype(&fans_plugin_load)      load;
    decltype(&fans_plugin_set_table) set_table;
    decltype(&fans_plugin_evaluate)  evaluate;
    decltype(&fans_plugin_free)      free_model;
};

// "NEML2" -> libfans_neml2.so, found next to libFANS.so through its RPATH
Plugin load_plugin(const string &name)
{
    string lib = "libfans_";
    for (char c : name)
        lib += char(std::tolower(static_cast<unsigned char>(c)));
    lib += ".so";

    void *h = dlopen(lib.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (!h)
        throw std::runtime_error("Could not load the '" + name + "' material plugin: " + dlerror() +
                                 "\nBuild FANS with -DFANS_" + name + "=ON.");

    Plugin p;
    p.load       = reinterpret_cast<decltype(p.load)>(dlsym(h, "fans_plugin_load"));
    p.set_table  = reinterpret_cast<decltype(p.set_table)>(dlsym(h, "fans_plugin_set_table"));
    p.evaluate   = reinterpret_cast<decltype(p.evaluate)>(dlsym(h, "fans_plugin_evaluate"));
    p.free_model = reinterpret_cast<decltype(p.free_model)>(dlsym(h, "fans_plugin_free"));
    if (!p.load || !p.set_table || !p.evaluate || !p.free_model)
        throw std::runtime_error(lib + " lacks a function of fans_plugin.h");
    return p;
}

// A model input FANS supplies from a dataset of the microstructure file
struct Field {
    bool           per_voxel; // data is [local voxel][size], else [phase][size]
    vector<double> data;
};

// Dataset `dset` of the microstructure file, next to the microstructure unless an
// absolute path; per voxel if it is [Z][Y][X][...] on the grid, else per phase [n_phase][...]
Field read_field(const Reader &reader, const string &name, const string &dset, int size)
{
    const string          path      = dset[0] == '/' ? dset : reader.MSGroup() + dset;
    const vector<hsize_t> shape     = Reader::DataShape(reader.ms_filename, path);
    const bool            per_voxel = shape.size() >= 3 && shape[0] == hsize_t(reader.dims[2]) && shape[1] == hsize_t(reader.dims[1]) &&
                                      shape[2] == hsize_t(reader.dims[0]);
    hsize_t               entry     = 1; // values per voxel or phase
    for (size_t i = per_voxel ? 3 : 1; i < shape.size(); ++i)
        entry *= shape[i];
    if (shape.empty() || entry != hsize_t(size))
        throw std::runtime_error("Field '" + name + "' needs " + std::to_string(size) + " values per voxel [Z][Y][X][...] or " +
                                 "per phase [n_phase][...], but dataset " + path + " has another shape.");

    Field f{per_voxel, {}};
    if (per_voxel) {
        f.data.resize(size_t(reader.local_n0) * reader.dims[1] * reader.dims[2] * size);
        reader.ReadSlab(f.data.data(), vector<int>(shape.begin() + 3, shape.end()), reader.ms_filename, path);
    } else {
        f.data = Reader::ReadData<double>(reader.ms_filename, path);
    }
    return f;
}

// Base is ThermalModel, SmallStrainMechModel or LargeStrainMechModel; the
// gradient and flux have n_str values per Gauss point
template <class Base, int howmany, int n_str>
class PluginModel : public Base {
  public:
    using Base::B;
    using Base::g0;
    using Base::n_gp;
    using Base::sigma;

    PluginModel(const Reader &reader)
        : Base(reader), plugin(load_plugin(reader.matmodel)), linear(reader.materialProperties.value("linear", false))
    {
        char msg[FANS_PLUGIN_MSGLEN] = {0};
        model                        = plugin.load(reader.materialProperties.dump().c_str(), msg, sizeof(msg));
        if (!model)
            throw std::runtime_error("Could not load the " + reader.matmodel + " material: " + msg);

        const json description = json::parse(msg);
        if (description.at("gradient") != n_str)
            throw std::runtime_error("The " + reader.matmodel + " material's gradient has " + description.at("gradient").dump() +
                                     " values, but this problem's has " + std::to_string(n_str) + ".");
        chunk        = std::max<size_t>(1, description.at("batch_size").get<size_t>() / n_gp);
        history_vars = description.at("history").get<vector<std::pair<string, int>>>();
        string history_names, field_names;
        for (const auto &[name, size] : history_vars) {
            n_history += size;
            history_names += (history_names.empty() ? "" : ", ") + name + " (" + std::to_string(size) + ")";
        }

        // Every field, read from the dataset "fields" names for it, goes to the plugin
        // once as a table of a row per voxel or per phase
        const auto model_fields = description.at("fields").get<vector<std::pair<string, int>>>();
        json       datasets     = reader.materialProperties.value("fields", json::object());
        for (size_t f = 0; f < model_fields.size(); ++f) {
            const auto &[name, size] = model_fields[f];
            if (!datasets.contains(name))
                throw std::runtime_error("The " + reader.matmodel + " material takes the field '" + name + "': name its dataset in \"fields\", e.g. \"fields\": {\"" + name + "\": \"<dataset>\"}.");
            const Field  field  = read_field(reader, name, datasets[name].get<string>(), size);
            const size_t n_rows = field.data.size() / size;
            datasets.erase(name);
            if (plugin.set_table(model, f, !field.per_voxel, n_rows, field.data.data(), msg, sizeof(msg)) != 0)
                throw std::runtime_error("Could not hand the " + reader.matmodel + " material its field '" + name + "': " + msg);
            if (!field.per_voxel)
                n_phase_rows = std::min(n_phase_rows, n_rows);
            field_names += (field_names.empty() ? "" : ", ") + name + (field.per_voxel ? " (per voxel)" : " (per phase)");
        }
        if (!datasets.empty())
            throw std::runtime_error("\"fields\" names '" + datasets.begin().key() + "', which is not a field of the " + reader.matmodel + " material.");

        Log::logger().info("# {} material: history {}, fields {}", reader.matmodel, history_names.empty() ? "none" : history_names,
                           field_names.empty() ? "none" : field_names);
    }

    ~PluginModel() override
    {
        plugin.free_model(model);
    }
    PluginModel(const PluginModel &)            = delete;
    PluginModel &operator=(const PluginModel &) = delete;

    // The history of every Gauss point, and the arrays of one plugin call
    void initializeInternalVariables(ptrdiff_t num_elements, int num_gauss_points) override
    {
        previous = trial = converged = vector<double>(num_elements * num_gauss_points * n_history, 0.0);
        gradient_buf.resize(chunk * n_gp * n_str);
        flux_buf.resize(chunk * n_gp * n_str);
        history_old_buf.resize(chunk * n_gp * n_history);
        history_new_buf.resize(chunk * n_gp * n_history);
        voxel_buf.resize(chunk * n_gp);
        phase_buf.resize(chunk * n_gp);
    }

    // A step has converged: its state is what the next step starts from
    void updateInternalVariables() override
    {
        previous.swap(converged);
        converged = trial;
    }

    // The fluxes at all Gauss points of an element, from the last evaluate_batch
    Eigen::Map<const VectorXd> fluxes(ptrdiff_t element_idx) const
    {
        if (!flux_cache)
            throw std::logic_error("PluginModel: an element sweep without Solver::evaluate_batched_stress before it");
        return {flux_cache + size_t(element_idx) * n_gp * n_str, n_gp * n_str};
    }

    // The flux at a Gauss point, from the last evaluate_batch
    void get_sigma(int i, int, ptrdiff_t element_idx) override
    {
        sigma.segment(i, n_str) = fluxes(element_idx).segment(i, n_str);
    }

    // The element's residual straight from its fluxes: the gradient and the loop over
    // the Gauss points of Matmodel::element_residual are not needed
    Matrix<double, howmany * 8, 1> &element_residual(Matrix<double, howmany * 8, 1> &, int, ptrdiff_t element_idx) override
    {
        this->res_e.noalias() = B.transpose() * fluxes(element_idx) * this->v_e / n_gp;
        return this->res_e;
    }

    bool wants_batch() const override
    {
        return true;
    }

    // "linear": true in the material properties; FANS trusts it
    bool is_linear() const override
    {
        return linear;
    }

    // The flux at all Gauss points of `elems`, into flux_all, and their state, into `trial`.
    // With tangent_all, the converged step once more, from the state it started from, for its tangent
    void evaluate_batch(const vector<ptrdiff_t> &elems, const unsigned short *phase, const double *ue_all,
                        double *flux_all, double *tangent_all) override
    {
        const size_t n_dof     = size_t(howmany) * 8; // nodal values per element
        const size_t per_elem  = size_t(n_gp) * n_str;
        const size_t hist_elem = size_t(n_gp) * n_history;
        flux_cache             = flux_all;

        const vector<double> &old = tangent_all ? previous : converged;
        vector<double>        tangent_buf(tangent_all ? chunk * per_elem * n_str : 0);

        for (size_t b = 0; b < elems.size(); b += chunk) {
            const size_t ne = std::min(chunk, elems.size() - b);

            for (size_t k = 0; k < ne; ++k) {
                const ptrdiff_t e = elems[b + k];
                Eigen::Map<VectorXd>(&gradient_buf[k * per_elem], per_elem).noalias() =
                    B * Eigen::Map<const VectorXd>(ue_all + e * n_dof, n_dof) + g0;
                std::fill_n(&voxel_buf[k * n_gp], n_gp, int(e));
                std::fill_n(&phase_buf[k * n_gp], n_gp, phase[e]);
                if (phase[e] >= n_phase_rows)
                    throw std::runtime_error("Phase " + std::to_string(phase[e]) + " has no row in the per-phase fields of its plugin material.");
                std::copy_n(old.data() + e * hist_elem, hist_elem, history_old_buf.data() + k * hist_elem);
                std::copy_n(trial.data() + e * hist_elem, hist_elem, history_new_buf.data() + k * hist_elem); // initial guess
            }

            char err[FANS_PLUGIN_MSGLEN] = {0};
            if (plugin.evaluate(model, ne * n_gp, this->time_old, this->time, gradient_buf.data(), voxel_buf.data(), phase_buf.data(),
                                history_old_buf.data(), flux_buf.data(), history_new_buf.data(), tangent_all ? tangent_buf.data() : nullptr,
                                err, sizeof(err)) != 0)
                throw std::runtime_error(string("Plugin material evaluation failed: ") + err);

            for (size_t k = 0; k < ne; ++k) {
                const ptrdiff_t e = elems[b + k];
                std::copy_n(&flux_buf[k * per_elem], per_elem, flux_all + e * per_elem);
                std::copy_n(history_new_buf.data() + k * hist_elem, hist_elem, trial.data() + e * hist_elem);
                if (tangent_all)
                    std::copy_n(&tangent_buf[k * per_elem * n_str], per_elem * n_str, tangent_all + e * per_elem * n_str);
            }
        }
    }

    // Writes every history variable, e.g. state/internal/Ep as state_internal_Ep
    void postprocess(Solver<howmany, n_str> &solver, Reader &reader, int load_idx, int time_idx) override
    {
        const auto &res     = reader.resultsToWrite;
        const bool  need    = std::find(res.begin(), res.end(), "internal_variables") != res.end();
        const bool  need_gp = std::find(res.begin(), res.end(), "internal_variables_gp") != res.end();
        if (!need && !need_gp)
            return;
        const size_t n_pts = size_t(solver.local_n0 * solver.n_y * solver.n_z) * n_gp;
        const string dir   = string(reader.dataset_name) + "/load" + std::to_string(load_idx) + "/time_step" + std::to_string(time_idx) + "/";

        size_t offset = 0;
        for (auto &[name, size] : history_vars) {
            string field = name;
            std::replace(field.begin(), field.end(), '/', '_');
            vector<double> gp(n_pts * size), elem(n_pts / n_gp * size, 0.0);
            for (size_t p = 0; p < n_pts; ++p)
                for (int c = 0; c < size; ++c) {
                    gp[p * size + c] = converged[p * n_history + offset + c];
                    elem[p / n_gp * size + c] += gp[p * size + c] / n_gp;
                }
            if (need)
                reader.WriteSlab(elem.data(), {size}, (dir + field).c_str());
            if (need_gp)
                reader.WriteSlab(gp.data(), {n_gp, size}, (dir + field + "_gp").c_str());
            offset += size;
        }
    }

    // A plugin cannot know a good reference stiffness, so the input file gives it
    Matrix<double, n_str, n_str> get_reference_stiffness() override
    {
        throw std::runtime_error("Plugin materials need the reference stiffness of the FFT solver in the input file: "
                                 "\"reference_material\", a symmetric positive semi-definite " +
                                 std::to_string(n_str) + "x" + std::to_string(n_str) + " matrix.");
    }

  private:
    Plugin                         plugin;
    FANSPluginModel               *model{nullptr};
    bool                           linear;                                                   // the flux is linear in the gradient, says the input file
    size_t                         chunk;                                                    // elements per plugin call
    vector<std::pair<string, int>> history_vars;                                             // name, doubles per point
    int                            n_history{0};                                             // history doubles per point
    size_t                         n_phase_rows{SIZE_MAX};                                   // phases below have a row in every per-phase field
    vector<double>                 trial, converged, previous;                               // history at every Gauss point
    vector<double>                 gradient_buf, flux_buf, history_old_buf, history_new_buf; // one plugin call
    vector<int>                    voxel_buf, phase_buf;
    const double                  *flux_cache{nullptr};
};

} // namespace

template <int howmany, int n_str>
Matmodel<howmany, n_str> *create_plugin_material(const Reader &reader)
{
    if constexpr (n_str == 3)
        return new PluginModel<ThermalModel, 1, 3>(reader);
    else if constexpr (n_str == 6)
        return new PluginModel<SmallStrainMechModel, 3, 6>(reader);
    else
        return new PluginModel<LargeStrainMechModel, 3, 9>(reader);
}

template Matmodel<1, 3> *create_plugin_material<1, 3>(const Reader &);
template Matmodel<3, 6> *create_plugin_material<3, 6>(const Reader &);
template Matmodel<3, 9> *create_plugin_material<3, 9>(const Reader &);

// Gauss-point stresses of all batching models at `u`, which the next sweep
// over the elements takes
template <int howmany, int n_str>
void Solver<howmany, n_str>::evaluate_batched_stress(double *u)
{
    constexpr size_t dof_elem = howmany * 8;
    const auto      &models   = matmanager->models;

    if (batch_elems.empty()) {
        const ptrdiff_t n_elem = local_n0 * n_y * n_z;
        batch_elems.assign(models.size(), {});
        for (size_t m = 0; m < models.size(); ++m)
            if (models[m]->wants_batch())
                for (ptrdiff_t e = 0; e < n_elem; ++e)
                    if (matmanager->get_info(ms[e]).model == models[m])
                        batch_elems[m].push_back(e);
        batch_ue.resize(size_t(n_elem) * dof_elem);
        batch_gp_stress.resize(size_t(n_elem) * models[0]->n_gp * n_str);
    }

    MPI_Sendrecv(u, n_y * n_z * howmany, MPI_DOUBLE, (world_rank + world_size - 1) % world_size, 0,
                 u + local_n0 * n_y * n_z * howmany, n_y * n_z * howmany, MPI_DOUBLE, (world_rank + 1) % world_size, 0, communicator, MPI_STATUS_IGNORE);

    iterateCubes<0>([&](ptrdiff_t *idx, ptrdiff_t *idxPadding) {
        if (matmanager->get_info(ms[idx[0]]).model->wants_batch())
            for (int i = 0; i < 8; ++i)
                for (int j = 0; j < howmany; ++j)
                    batch_ue[idx[0] * dof_elem + howmany * i + j] = u[howmany * idx[i] + j] - u[howmany * idx[0] + j];
    });

    for (size_t m = 0; m < models.size(); ++m)
        if (models[m]->wants_batch())
            models[m]->evaluate_batch(batch_elems[m], ms, batch_ue.data(), batch_gp_stress.data());
}

template void Solver<1, 3>::evaluate_batched_stress(double *);
template void Solver<3, 6>::evaluate_batched_stress(double *);
template void Solver<3, 9>::evaluate_batched_stress(double *);
