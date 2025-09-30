#ifndef MIXED_BC_H
#define MIXED_BC_H

// ============================================================================
//  mixedBCs.h
//  --------------------------------------------------------------------------
//  • Holds MixedBC + LoadCase structs
//  • Provides:   MixedBC::from_json(), finalize()
// ============================================================================

#include <Eigen/Dense>
using namespace Eigen;

// ---------------------------------------------------------------------------
struct MixedBC {
    /*  Index sets (0‑based)  */
    VectorXi idx_E; // strain‑controlled components
    VectorXi idx_F; // stress‑controlled components

    /*  Time paths  */
    MatrixXd F_E_path; // (#steps × |E|)
    MatrixXd P_F_path; // (#steps × |F|)

    /*  Projectors & auxiliary matrix  */
    MatrixXd Q_E, Q_F, M; // M = (Q_Fᵀ C0 Q_F)⁻¹

    // ------------------------------------------------------------
    // build Q_E, Q_F, M
    // ------------------------------------------------------------
    void finalize(const MatrixXd &C0)
    {
        const int n_str = static_cast<int>(C0.rows());

        Q_E = MatrixXd::Zero(n_str, idx_E.size());
        for (int c = 0; c < idx_E.size(); ++c)
            Q_E(idx_E(c), c) = 1.0;

        Q_F = MatrixXd::Zero(n_str, idx_F.size());
        for (int c = 0; c < idx_F.size(); ++c)
            Q_F(idx_F(c), c) = 1.0;

        if (idx_F.size() > 0)
            M = (Q_F.transpose() * C0 * Q_F).inverse();
        else
            M.resize(0, 0); // pure‑strain case
    }

    // ------------------------------------------------------------
    // Factory: parse MixedBC from a JSON object
    // ------------------------------------------------------------
    static MixedBC from_json(const json &jc, int n_str)
    {
        auto toEigen = [](const vector<int> &v) {
            VectorXi e(v.size());
            for (size_t i = 0; i < v.size(); ++i)
                e(static_cast<int>(i)) = v[i];
            return e;
        };

        MixedBC bc;

        if (!jc.contains("strain_indices") || !jc.contains("stress_indices"))
            throw runtime_error("mixed BC: strain_indices or stress_indices missing");

        bc.idx_E = toEigen(jc["strain_indices"].get<vector<int>>());
        bc.idx_F = toEigen(jc["stress_indices"].get<vector<int>>());

        // ---- sanity: disjoint + complementary -----
        vector<char> present(n_str, 0);
        for (int i = 0; i < bc.idx_E.size(); ++i) {
            int k = bc.idx_E(i);
            if (k < 0 || k >= n_str)
                throw runtime_error("strain index out of range");
            present[k] = 1;
        }
        for (int i = 0; i < bc.idx_F.size(); ++i) {
            int k = bc.idx_F(i);
            if (k < 0 || k >= n_str)
                throw runtime_error("stress index out of range");
            if (present[k])
                throw runtime_error("index appears in both strain_indices and stress_indices");
            present[k] = 1;
        }
        for (int k = 0; k < n_str; ++k)
            if (!present[k])
                throw runtime_error("each component must be either strain‑ or stress‑controlled");

        // ---- parse 2‑D arrays (allow empty for |E|==0 etc.) -----
        auto get2D = [](const json &arr) {
            return arr.get<vector<vector<double>>>();
        };

        vector<vector<double>> strain_raw, stress_raw;
        size_t                 n_steps = 0;

        if (bc.idx_E.size() > 0) {
            if (!jc.contains("strain"))
                throw runtime_error("strain array missing");
            strain_raw = get2D(jc["strain"]);
            n_steps    = strain_raw.size();
        }
        if (bc.idx_F.size() > 0) {
            if (!jc.contains("stress"))
                throw runtime_error("stress array missing");
            stress_raw = get2D(jc["stress"]);
            n_steps    = max(n_steps, stress_raw.size());
        }
        if (n_steps == 0)
            throw runtime_error("mixed BC: at least one of strain/stress must have timesteps");

        // default‑fill for missing part (constant 0)
        if (strain_raw.empty())
            strain_raw.resize(n_steps, vector<double>(0));
        if (stress_raw.empty())
            stress_raw.resize(n_steps, vector<double>(0));

        // consistency checks & build Eigen matrices
        bc.F_E_path = MatrixXd::Zero(n_steps, bc.idx_E.size());
        bc.P_F_path = MatrixXd::Zero(n_steps, bc.idx_F.size());

        for (size_t t = 0; t < n_steps; ++t) {
            if (strain_raw[t].size() != static_cast<size_t>(bc.idx_E.size()))
                throw runtime_error("strain row length mismatch");
            for (int c = 0; c < bc.idx_E.size(); ++c)
                bc.F_E_path(static_cast<int>(t), c) = (bc.idx_E.size() ? strain_raw[t][c] : 0.0);

            if (stress_raw[t].size() != static_cast<size_t>(bc.idx_F.size()))
                throw runtime_error("stress row length mismatch");
            for (int c = 0; c < bc.idx_F.size(); ++c)
                bc.P_F_path(static_cast<int>(t), c) = (bc.idx_F.size() ? stress_raw[t][c] : 0.0);
        }
        return bc;
    }
};

// ---------------------------------------------------------------------------
//  LoadCase : covers both legacy and mixed formats (used by Reader)
// ---------------------------------------------------------------------------
struct LoadCase {
    bool                   mixed = false;
    vector<vector<double>> g0_path;     // legacy pure‑strain
    MixedBC                mbc;         // mixed BC data
    size_t                 n_steps = 0; // number of time steps
};

// ---------------------------------------------------------------------------
//  Helper mix‑in with enable/update for Solver (header‑only)
// ---------------------------------------------------------------------------
template <int howmany>
struct MixedBCController {
    bool mixed_active = false;

  protected:
    MixedBC  mbc_local;
    size_t   step_idx = 0;
    VectorXd g0_vec; // current macro strain (size n_str)

    // call from user code after v_u update each iteration
    template <typename SolverType>
    void update(SolverType &solver)
    {
        if (!mixed_active)
            return;

        VectorXd Pbar = solver.get_homogenized_stress();

        VectorXd PF = (mbc_local.idx_F.size() ? mbc_local.P_F_path.row(step_idx).transpose() : VectorXd());

        if (mbc_local.idx_F.size()) {
            VectorXd rhs     = PF - mbc_local.Q_F.transpose() * Pbar;
            VectorXd delta_E = mbc_local.M * rhs;
            g0_vec += mbc_local.Q_F * delta_E;
        }

        vector<double> gvec(g0_vec.data(), g0_vec.data() + g0_vec.size());
        solver.matmodel->setGradient(gvec);
    }

    template <typename SolverType>
    void activate(SolverType &solver, const MixedBC &mbc_in, size_t t)
    {
        mixed_active = true;
        mbc_local    = mbc_in;
        step_idx     = t;
        mbc_local.finalize(solver.matmodel->kapparef_mat);

        int n_str = solver.matmodel->n_str;
        g0_vec    = VectorXd::Zero(n_str);
        if (mbc_local.idx_E.size())
            g0_vec += mbc_local.Q_E * mbc_local.F_E_path.row(t).transpose();

        vector<double> gvec(g0_vec.data(), g0_vec.data() + n_str);
        solver.matmodel->setGradient(gvec);

        // Do one update to set the initial strain
        solver.updateMixedBC();
    }
};

#endif // MIXED_BC_H
