// Material evaluated by a plugin loaded at run time (fans_plugin.h), and the
// solver's batched stress evaluation for it. Kept out of the headers so the
// native models' machine code is unaffected.

#include "general.h"
#include "matmodel.h"
#include "solver.h"
#include "fans_plugin.h"

#include <algorithm>
#include <cctype>
#include <dlfcn.h>

namespace {

struct Plugin {
    FANSPluginModel *(*load)(const char *, const char *, int *, int *, char *, size_t);
    int (*evaluate)(FANSPluginModel *, size_t, const double *, const double *, const double *,
                    double *, double *, char *, size_t);
    void (*free_model)(FANSPluginModel *);
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

    auto sym = [&](auto &f, const char *n) {
        if (!(f = reinterpret_cast<std::remove_reference_t<decltype(f)>>(dlsym(h, n))))
            throw std::runtime_error(lib + " is missing '" + n + "'");
    };
    Plugin p;
    sym(p.load, "fans_plugin_load");
    sym(p.evaluate, "fans_plugin_evaluate");
    sym(p.free_model, "fans_plugin_free");
    return p;
}

class PluginSmallStrainMechModel : public SmallStrainMechModel {
  public:
    PluginSmallStrainMechModel(const Reader &reader)
        : SmallStrainMechModel(reader), plugin(load_plugin(reader.matmodel))
    {
        const string spec   = reader.materialProperties.at("artifact").get<string>();
        const string device = reader.materialProperties.value("device", string("cpu"));

        char msg[FANS_PLUGIN_MSGLEN] = {0};

        model = plugin.load(spec.c_str(), device.c_str(), &n_state, &wants_orientation, msg, sizeof(msg));
        if (!model)
            throw std::runtime_error("Could not load material '" + spec + "': " + msg);
        n_mat = 1;

        batch_points = reader.materialProperties.value("batch_size", device == "cpu" ? size_t{1024} : size_t{65536});

        if (wants_orientation) {
            string path(reader.ms_datasetname);
            path.replace(path.find_last_of('/') + 1, string::npos, "rotation_matrices");
            try {
                H5::H5File  file(reader.ms_filename, H5F_ACC_RDONLY);
                H5::DataSet ds = file.openDataSet(path);
                hsize_t     dims[3];
                ds.getSpace().getSimpleExtentDims(dims);
                grain_rot.resize(9 * dims[0]);
                ds.read(grain_rot.data(), H5::PredType::NATIVE_DOUBLE);
            } catch (const H5::Exception &) {
                throw std::runtime_error("Material '" + spec + "' needs crystal orientations, but " +
                                         reader.ms_filename + " has no " + path + ".");
            }
        }

        string vars;
        state_vars = json::parse(msg).get<vector<std::pair<string, int>>>();
        for (const auto &[name, size] : state_vars)
            vars += (vars.empty() ? "" : ", ") + name + " (" + std::to_string(size) + ")";
        Log::logger().info("# Plugin material '{}' on {}: history {}{}", spec, device, vars.empty() ? "none" : vars,
                           wants_orientation ? ", orientations from rotation_matrices" : "");
    }

    ~PluginSmallStrainMechModel() override
    {
        plugin.free_model(model);
    }
    PluginSmallStrainMechModel(const PluginSmallStrainMechModel &)            = delete;
    PluginSmallStrainMechModel &operator=(const PluginSmallStrainMechModel &) = delete;

    void initializeInternalVariables(ptrdiff_t num_elements, int num_gauss_points) override
    {
        state = state_t = vector<double>(num_elements * num_gauss_points * n_state, 0.0);
    }

    void updateInternalVariables() override
    {
        sig_cache = nullptr;
        state.swap(state_t);
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        if (i != 0)
            return;
        if (!sig_cache)
            throw std::logic_error("PluginSmallStrainMechModel::get_sigma without Solver::evaluate_batched_stress");
        std::copy_n(sig_cache + size_t(element_idx) * n_gp * 6, n_gp * 6, sigma.data());
    }

    bool wants_batch() const override
    {
        return true;
    }

    // A finished_step restarts from the step's initial history (in `state`
    // after the swap), so a rate-dependent flow is not applied twice.
    void evaluate_batch(const vector<ptrdiff_t> &elems, const unsigned short *phase, const double *ue_all,
                        double *sig_all, bool finished_step) override
    {
        constexpr size_t n_dof = 24;
        const size_t     per_e = size_t(n_gp) * 6;
        const size_t     st_e  = size_t(n_gp) * n_state;
        const size_t     chunk = std::max<size_t>(1, batch_points / n_gp);

        in_.resize(chunk * per_e);
        out_.resize(chunk * per_e);
        or_.resize(chunk * n_gp * 9);
        so_.resize(chunk * st_e);
        sn_.resize(chunk * st_e);
        sig_cache = sig_all;

        const vector<double> &history_old = finished_step ? state : state_t;

        for (size_t b = 0; b < elems.size(); b += chunk) {
            const size_t ne = std::min(chunk, elems.size() - b);

            for (size_t k = 0; k < ne; ++k) {
                const ptrdiff_t e = elems[b + k];
                Eigen::Map<VectorXd>(&in_[k * per_e], per_e).noalias() =
                    B * Eigen::Map<const Matrix<double, n_dof, 1>>(ue_all + size_t(e) * n_dof) + g0;
                if (wants_orientation) {
                    const size_t g = phase[e];
                    if (9 * (g + 1) > grain_rot.size())
                        throw std::runtime_error("Grain " + std::to_string(g) + " has no entry in rotation_matrices.");
                    for (int p = 0; p < n_gp; ++p)
                        std::copy_n(&grain_rot[9 * g], 9, &or_[(k * n_gp + p) * 9]);
                }
                std::copy_n(history_old.data() + e * st_e, st_e, so_.data() + k * st_e);
                std::copy_n(state.data() + e * st_e, st_e, sn_.data() + k * st_e); // initial guess: the latest trial
            }

            evaluate(ne * n_gp, in_.data(), or_.data(), so_.data(), out_.data(), sn_.data());

            for (size_t k = 0; k < ne; ++k) {
                const ptrdiff_t e = elems[b + k];
                std::copy_n(out_.data() + k * per_e, per_e, sig_all + e * per_e);
                if (!finished_step)
                    std::copy_n(sn_.data() + k * st_e, st_e, state.data() + e * st_e);
            }
        }
    }

    // Writes every history variable, e.g. state/internal/Ep as state_internal_Ep
    void postprocess(Solver<3, 6> &solver, Reader &reader, int load_idx, int time_idx) override
    {
        const auto &res     = reader.resultsToWrite;
        const bool  need    = std::find(res.begin(), res.end(), "internal_variables") != res.end();
        const bool  need_gp = std::find(res.begin(), res.end(), "internal_variables_gp") != res.end();
        if (!need && !need_gp)
            return;
        const size_t n_pts = size_t(solver.local_n0 * solver.n_y * solver.n_z) * n_gp;
        const string dir   = string(reader.dataset_name) + "/load" + std::to_string(load_idx) + "/time_step" + std::to_string(time_idx) + "/";

        size_t offset = 0;
        for (auto &[name, size] : state_vars) {
            string field = name;
            std::replace(field.begin(), field.end(), '/', '_');
            vector<double> gp(n_pts * size), elem(n_pts / n_gp * size, 0.0);
            for (size_t p = 0; p < n_pts; ++p)
                for (int c = 0; c < size; ++c) {
                    gp[p * size + c] = state_t[p * n_state + offset + c];
                    elem[p / n_gp * size + c] += gp[p * size + c] / n_gp;
                }
            if (need)
                reader.WriteSlab(elem.data(), {size}, (dir + field).c_str());
            if (need_gp)
                reader.WriteSlab(gp.data(), {n_gp, size}, (dir + field + "_gp").c_str());
            offset += size;
        }
    }

    // Initial stiffness
    Matrix<double, 6, 6> get_reference_stiffness() override
    {
        constexpr double     h   = 1e-6;
        Matrix<double, 6, 6> eps = h * Matrix<double, 6, 6>::Identity(), sig;
        vector<double>       s_old(6 * n_state, 0.0), s_new(s_old.size()), rot(6 * 9, 0.0);
        for (int p = 0; p < 6; ++p)
            rot[9 * p] = rot[9 * p + 4] = rot[9 * p + 8] = 1.0;
        evaluate(6, eps.data(), rot.data(), s_old.data(), sig.data(), s_new.data());
        return (sig + sig.transpose()) / (2 * h);
    }

  private:
    void evaluate(size_t n, const double *eps, const double *rot, const double *s_old, double *sig, double *s_new)
    {
        char err[FANS_PLUGIN_MSGLEN] = {0};
        if (plugin.evaluate(model, n, eps, wants_orientation ? rot : nullptr, n_state ? s_old : nullptr, sig,
                            n_state ? s_new : nullptr, err, sizeof(err)) != 0)
            throw std::runtime_error(string("Plugin material evaluation failed: ") + err);
    }

    Plugin                         plugin;
    FANSPluginModel               *model{nullptr};
    int                            n_state{0};
    int                            wants_orientation{0};
    vector<std::pair<string, int>> state_vars; // history variables: name, size
    size_t                         batch_points;
    vector<double>                 grain_rot; // [grain][3][3]
    vector<double>                 state;     // trial history
    vector<double>                 state_t;   // converged history
    vector<double>                 in_, out_, so_, sn_, or_;
    const double                  *sig_cache{nullptr};
};

} // namespace

Matmodel<3, 6> *create_plugin_material(const Reader &reader)
{
    return new PluginSmallStrainMechModel(reader);
}

// Gauss-point stresses of all batching models at `u`, served to the next
// per-element sweep through get_sigma
template <int howmany, int n_str>
void Solver<howmany, n_str>::evaluate_batched_stress(double *u, bool finished_step)
{
    constexpr size_t dof_elem = howmany * 8;
    const auto      &models   = matmanager->models;

    if (batch_elems.empty()) {
        const ptrdiff_t n_elem = local_n0 * n_y * n_z;
        batch_elems.assign(models.size(), {});
        for (ptrdiff_t e = 0; e < n_elem; ++e) {
            Matmodel<howmany, n_str> *model = matmanager->get_info(ms[e]).model;
            if (model->wants_batch())
                batch_elems[std::find(models.begin(), models.end(), model) - models.begin()].push_back(e);
        }
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
            models[m]->evaluate_batch(batch_elems[m], ms, batch_ue.data(), batch_gp_stress.data(), finished_step);
}

template void Solver<1, 3>::evaluate_batched_stress(double *, bool);
template void Solver<3, 6>::evaluate_batched_stress(double *, bool);
template void Solver<3, 9>::evaluate_batched_stress(double *, bool);
