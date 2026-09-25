/*
 * Material whose constitutive response comes from a plugin loaded at run
 * time (see fans_plugin.h), plus the solver's batched evaluation of it.
 *
 * All of it lives in this translation unit rather than in headers on
 * purpose: code added to the headers the solver's hot templates are compiled
 * with changes their LTO inlining, and so the native models' machine code,
 * even in runs that never load a plugin.
 */

#include "general.h"
#include "matmodel.h"
#include "solver.h"
#include "fans_plugin.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <dlfcn.h>

Matmodel<3, 6> *create_plugin_material(const Reader &reader);

namespace {

struct Plugin {
    FANSPluginModel *(*load)(const char *, int *, int *, char *, size_t);
    int (*evaluate)(FANSPluginModel *, size_t, const double *, const double *, const double *,
                    double *, double *, char *, size_t);
    void (*free_model)(FANSPluginModel *);
};

/// "NEML2" -> libfans_neml2.so, found through FANS's own RPATH, i.e. next to
/// libFANS.so in both the build and the install tree.
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
#define FANS_SYM(f, n)                        \
    if (!(p.f = (decltype(p.f)) dlsym(h, n))) \
    throw std::runtime_error(lib + " is missing '" n "'")
    FANS_SYM(load, "fans_plugin_load");
    FANS_SYM(evaluate, "fans_plugin_evaluate");
    FANS_SYM(free_model, "fans_plugin_free");
#undef FANS_SYM
    return p;
}

class PluginSmallStrainMechModel : public SmallStrainMechModel {
  public:
    PluginSmallStrainMechModel(const Reader &reader)
        : SmallStrainMechModel(reader), plugin(load_plugin(reader.matmodel))
    {
        if (!reader.materialProperties.contains("artifact"))
            throw std::runtime_error("Material model '" + reader.matmodel +
                                     "' needs an \"artifact\" naming the compiled material.");
        const string spec = reader.materialProperties["artifact"].get<string>();

        char err[FANS_PLUGIN_ERRLEN] = {0};

        model = plugin.load(spec.c_str(), &n_state, &wants_orientation, err, sizeof(err));
        if (!model)
            throw std::runtime_error("Could not load material '" + spec + "': " + err);
        n_mat = 1; // parameters come from the artifact, so one model per group

        batch_points = reader.materialProperties.value("batch_size", size_t{1024});
        if (batch_points == 0)
            throw std::runtime_error("batch_size must be positive.");

        if (wants_orientation) {
            // One crystal-to-sample rotation per grain, next to the microstructure
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

        Log::logger().info("# Plugin material '{}': {} history variable(s) per Gauss point{}", spec,
                           n_state, wants_orientation ? ", orientations from rotation_matrices" : "");
    }

    ~PluginSmallStrainMechModel() override
    {
        if (model)
            plugin.free_model(model);
    }
    PluginSmallStrainMechModel(const PluginSmallStrainMechModel &)            = delete; // owns `model`
    PluginSmallStrainMechModel &operator=(const PluginSmallStrainMechModel &) = delete;

    void initializeInternalVariables(ptrdiff_t num_elements, int num_gauss_points) override
    {
        if (n_state)
            state = state_t = vector<double>(num_elements * num_gauss_points * n_state, 0.0);
    }

    /// Promote trial history to converged, as the native models do with `_t`.
    /// The swap leaves the history the step started from in `state`.
    void updateInternalVariables() override
    {
        sig_cache = nullptr;
        if (n_state)
            state.swap(state_t);
    }

    /// Every element sweep is prepared by Solver::evaluate_batched_stress, so
    /// this copies that batched call's Gauss-point stresses.
    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        if (i != 0)
            return;
        if (!sig_cache)
            throw std::logic_error("PluginSmallStrainMechModel::get_sigma without Solver::evaluate_batched_stress");
        const double *sig = sig_cache + size_t(element_idx) * n_gp * 6;
        std::copy(sig, sig + size_t(n_gp) * 6, sigma.data());
    }

    bool wants_batch() const override
    {
        return true;
    }

    /**
     * Evaluates `elems` in chunks of batch_points material points. Too small
     * and per-call overhead dominates; too large and the model's temporaries
     * leave cache, and the ranks of a full node compete for memory bandwidth
     * (J2 scales to ~4096 points, crystal plasticity to ~512).
     *
     * A `finished_step` (postprocessing, after the solve committed the
     * history) starts again from the history that step started from, still
     * in `state`, and leaves the committed one alone: re-evaluating from the
     * committed history would apply a rate-dependent model's flow twice.
     */
    void evaluate_batch(const vector<ptrdiff_t> &elems, const unsigned short *phase, const double *ue_all,
                        double *sig_all, bool finished_step) override
    {
        constexpr size_t n_dof = 24;
        const size_t     per_e = size_t(n_gp) * 6;
        const size_t     st_e  = size_t(n_gp) * n_state;
        const size_t     chunk = std::max<size_t>(1, batch_points / n_gp);

        in_.resize(chunk * per_e);
        out_.resize(chunk * per_e);
        or_.resize(wants_orientation ? chunk * n_gp * 9 : 0);
        so_.resize(chunk * st_e);
        sn_.resize(chunk * st_e);
        sig_cache = sig_all;

        const vector<double> &history_old = finished_step ? state : state_t;

        for (size_t b = 0; b < elems.size(); b += chunk) {
            const size_t ne = std::min(chunk, elems.size() - b);

            for (size_t k = 0; k < ne; ++k) {
                const ptrdiff_t e = elems[b + k];
                eps.noalias()     = B * Eigen::Map<const Matrix<double, n_dof, 1>>(ue_all + size_t(e) * n_dof) + g0;
                std::copy(eps.data(), eps.data() + per_e, &in_[k * per_e]);
                if (wants_orientation)
                    fill_orientation(phase[e], &or_[k * n_gp * 9]);
                if (n_state)
                    std::memcpy(&so_[k * st_e], &history_old[size_t(e) * st_e], st_e * sizeof(double));
            }

            char err[FANS_PLUGIN_ERRLEN] = {0};
            if (plugin.evaluate(model, ne * n_gp, in_.data(), wants_orientation ? or_.data() : nullptr,
                                n_state ? so_.data() : nullptr, out_.data(), n_state ? sn_.data() : nullptr,
                                err, sizeof(err)) != 0)
                throw std::runtime_error(string("Plugin material evaluation failed: ") + err);

            for (size_t k = 0; k < ne; ++k) {
                const ptrdiff_t e = elems[b + k];
                std::copy(&out_[k * per_e], &out_[(k + 1) * per_e], sig_all + size_t(e) * per_e);
                if (n_state && !finished_step)
                    std::memcpy(&state[size_t(e) * st_e], &sn_[k * st_e], st_e * sizeof(double));
            }
        }
    }

    /// Forward-difference tangent at zero strain and history (and the
    /// identity orientation), symmetrised because FANS takes a Cholesky
    /// factorisation of it. Linear models are exact; for others
    /// "reference_material" in the input overrides it.
    Matrix<double, 6, 6> get_reference_stiffness() override
    {
        constexpr double     h      = 1e-6;
        Matrix<double, 6, 7> eps_fd = Matrix<double, 6, 7>::Zero(), sig_fd; // column 0: base point
        vector<double>       s_old(7 * n_state, 0.0), s_new(s_old.size()), identity(7 * 9, 0.0);
        char                 err[FANS_PLUGIN_ERRLEN] = {0};
        eps_fd.rightCols<6>().diagonal().setConstant(h);
        for (int p = 0; p < 7; ++p)
            identity[9 * p] = identity[9 * p + 4] = identity[9 * p + 8] = 1.0;
        if (plugin.evaluate(model, 7, eps_fd.data(), wants_orientation ? identity.data() : nullptr,
                            n_state ? s_old.data() : nullptr, sig_fd.data(), n_state ? s_new.data() : nullptr,
                            err, sizeof(err)) != 0)
            throw std::runtime_error(string("Plugin material has no reference stiffness: ") + err +
                                     "\nGive \"reference_material\" in the input instead.");
        const Matrix<double, 6, 6> C = (sig_fd.rightCols<6>().colwise() - sig_fd.col(0)) / h;
        return 0.5 * (C + C.transpose());
    }

  private:
    /// Orientation of grain g, repeated for each Gauss point of an element.
    void fill_orientation(size_t g, double *dst) const
    {
        if (9 * (g + 1) > grain_rot.size())
            throw std::runtime_error("Grain " + std::to_string(g) + " has no entry in rotation_matrices.");
        for (int p = 0; p < n_gp; ++p)
            std::copy_n(&grain_rot[9 * g], 9, dst + 9 * p);
    }

    Plugin           plugin;
    FANSPluginModel *model{nullptr};
    int              n_state{0};
    int              wants_orientation{0};
    size_t           batch_points;             //!< material points per plugin call
    vector<double>   grain_rot;                //!< [grain][3][3] crystal-to-sample rotations
    vector<double>   state;                    //!< trial history
    vector<double>   state_t;                  //!< previous converged history
    vector<double>   in_, out_, so_, sn_, or_; //!< chunk scratch
    const double    *sig_cache{nullptr};       //!< Gauss-point stresses of the last stress recovery
};

} // namespace

Matmodel<3, 6> *create_plugin_material(const Reader &reader)
{
    return new PluginSmallStrainMechModel(reader);
}

/**
 * Evaluates every batching model's Gauss-point stresses at displacements `u`
 * in one batched call per model, so that the following per-element sweep
 * (residual, homogenized stress or postprocessing) is served from the
 * model's cache instead of calling the plugin once per element.
 * `finished_step` asks for the stresses the last solve converged to
 * (postprocessing, after the history was committed).
 */
template <int howmany, int n_str>
void Solver<howmany, n_str>::evaluate_batched_stress(double *u, bool finished_step)
{
    constexpr size_t dof_elem = howmany * 8;
    const auto      &models   = matmanager->models;

    if (batch_elems.empty()) { // first use: each batching model's elements, and the buffers
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

// Instantiated here, not in solver.h, to keep this code out of the solver's
// translation unit. The solver calls it behind a runtime flag, so every
// solver needs it even though only <3,6> can have a batching model.
template void Solver<1, 3>::evaluate_batched_stress(double *, bool);
template void Solver<3, 6>::evaluate_batched_stress(double *, bool);
template void Solver<3, 9>::evaluate_batched_stress(double *, bool);
