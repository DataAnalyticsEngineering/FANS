#ifndef FINITESTRAINJ2PLASTICITY_H
#define FINITESTRAINJ2PLASTICITY_H

#include "LargeStrainMechModel.h"
#include "solver.h"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>

/**
 * @brief Rate-independent finite-strain J2 plasticity with linear isotropic hardening
 *
 * This Simo-type model uses logarithmic elastic strain and stores the inverse
 * deformation gradient, elastic Finger tensor, and accumulated equivalent
 * plastic strain at every integration point. The spectral radial return
 * updates Kirchhoff stress, which is converted to first Piola stress for FANS.
 *
 * Material parameters:
 * - bulk_modulus
 * - shear_modulus
 * - yield_stress
 * - isotropic_hardening_parameter
 */
class FiniteStrainJ2Plasticity : public LargeStrainMechModel {
  public:
    FiniteStrainJ2Plasticity(const Reader &reader)
        : LargeStrainMechModel(reader)
    {
        vector<double> bulk;
        vector<double> shear;
        vector<double> yield;
        vector<double> hardening;

        try {
            bulk      = reader.materialProperties.at("bulk_modulus").get<vector<double>>();
            shear     = reader.materialProperties.at("shear_modulus").get<vector<double>>();
            yield     = reader.materialProperties.at("yield_stress").get<vector<double>>();
            hardening = reader.materialProperties.at("isotropic_hardening_parameter").get<vector<double>>();
        } catch (const std::exception &e) {
            throw std::runtime_error("Missing or invalid FiniteStrainJ2Plasticity properties: " + std::string(e.what()));
        }

        n_mat = static_cast<int>(bulk.size());
        if (n_mat == 0 || shear.size() != bulk.size() || yield.size() != bulk.size() || hardening.size() != bulk.size())
            throw std::runtime_error("FiniteStrainJ2Plasticity property arrays must have equal, non-zero length.");

        params.reserve(n_mat);
        for (int m = 0; m < n_mat; ++m) {
            const double K      = bulk[m];
            const double mu     = shear[m];
            const double sigma0 = yield[m];
            const double H      = hardening[m];

            if (!std::isfinite(K) || !std::isfinite(mu) || !std::isfinite(sigma0) || !std::isfinite(H) ||
                K <= 0.0 || mu <= 0.0 || sigma0 < 0.0 || H < 0.0)
                throw std::runtime_error("FiniteStrainJ2Plasticity requires finite K > 0, mu > 0, yield_stress >= 0, and hardening >= 0.");

            const double three_mu = 3.0 * mu;
            params.push_back({mu, sigma0, H, 0.5 * K, three_mu, three_mu + H});
        }
    }

    void initializeInternalVariables(ptrdiff_t num_elements, int num_gauss_points) override
    {
        const ptrdiff_t n_states = num_elements * num_gauss_points;

        F_inv.assign(n_states, Matrix3d::Identity());
        F_inv_t.assign(n_states, Matrix3d::Identity());
        be.assign(n_states, Matrix3d::Identity());
        be_t.assign(n_states, Matrix3d::Identity());
        ep   = VectorXd::Zero(n_states);
        ep_t = VectorXd::Zero(n_states);
    }

    void updateInternalVariables() override
    {
        F_inv.swap(F_inv_t);
        be.swap(be_t);
        ep.swap(ep_t);
    }

    Matrix<double, 9, 9> get_reference_stiffness() override
    {
        double lambda = 0.0;
        double mu     = 0.0;
        for (const Params &p : params) {
            lambda += 2.0 * p.half_K - 2.0 * p.mu / 3.0;
            mu += p.mu;
        }
        lambda /= static_cast<double>(n_mat);
        mu /= static_cast<double>(n_mat);

        Matrix<double, 6, 6> C = Matrix<double, 6, 6>::Zero();
        C.topLeftCorner(3, 3).setConstant(lambda);
        C.diagonal().array() += 2.0 * mu;
        return compute_spatial_tangent(Matrix3d::Identity(), Matrix3d::Zero(), C);
    }

    void postprocess(Solver<3, 9> &solver, Reader &reader, int load_idx, int time_idx) override
    {
        const ptrdiff_t n_elem  = solver.local_n0 * solver.n_y * solver.n_z;
        const auto     &results = reader.resultsToWrite;

        const bool need_ep    = std::find(results.begin(), results.end(), "equivalent_plastic_strain") != results.end();
        const bool need_ep_gp = std::find(results.begin(), results.end(), "equivalent_plastic_strain_gp") != results.end();
        const bool need_be    = std::find(results.begin(), results.end(), "elastic_finger") != results.end();
        const bool need_be_gp = std::find(results.begin(), results.end(), "elastic_finger_gp") != results.end();

        if (!(need_ep || need_ep_gp || need_be || need_be_gp))
            return;

        VectorXd ep_elem, ep_gp_data, be_elem, be_gp_data;
        if (need_ep)
            ep_elem.resize(n_elem);
        if (need_ep_gp)
            ep_gp_data.resize(n_elem * n_gp);
        if (need_be)
            be_elem.resize(n_elem * 9);
        if (need_be_gp)
            be_gp_data.resize(n_elem * n_gp * 9);

        for (ptrdiff_t elem = 0; elem < n_elem; ++elem) {
            const ptrdiff_t offset = elem * n_gp;

            if (need_ep)
                ep_elem(elem) = ep_t.segment(offset, n_gp).mean();

            if (need_be) {
                Matrix3d be_mean = Matrix3d::Zero();
                for (int gp = 0; gp < n_gp; ++gp)
                    be_mean += be_t[offset + gp];
                be_mean /= static_cast<double>(n_gp);

                for (int c = 0; c < 9; ++c)
                    be_elem(9 * elem + c) = be_mean(c / 3, c % 3);
            }

            if (need_ep_gp || need_be_gp) {
                for (int gp = 0; gp < n_gp; ++gp) {
                    const ptrdiff_t state = offset + gp;
                    if (need_ep_gp)
                        ep_gp_data(state) = ep_t(state);
                    if (need_be_gp)
                        for (int c = 0; c < 9; ++c)
                            be_gp_data(9 * state + c) = be_t[state](c / 3, c % 3);
                }
            }
        }

        if (need_ep)
            reader.writeSlab("equivalent_plastic_strain", load_idx, time_idx, ep_elem.data(), {1});
        if (need_ep_gp)
            reader.writeSlab("equivalent_plastic_strain_gp", load_idx, time_idx, ep_gp_data.data(), {n_gp});
        if (need_be)
            reader.writeSlab("elastic_finger", load_idx, time_idx, be_elem.data(), {9});
        if (need_be_gp)
            reader.writeSlab("elastic_finger_gp", load_idx, time_idx, be_gp_data.data(), {n_gp, 9});
    }

  protected:
    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        const ptrdiff_t state = element_idx * n_gp + i / 9;
        const Params   &p     = params[mat_index];
        const Matrix3d  F     = extract_F(i);

        if (const double J = F.determinant(); !std::isfinite(J) || J <= 0.0)
            throw std::runtime_error("FiniteStrainJ2Plasticity requires finite det(F) > 0.");

        const Matrix3d Finv = F.inverse();
        const Matrix3d Finc = F * F_inv_t[state];
        Matrix3d       tmp;
        Matrix3d       be_tr;
        tmp.noalias()   = be_t[state] * Finc.transpose();
        be_tr.noalias() = Finc * tmp;
        be_tr           = 0.5 * (be_tr + be_tr.transpose()).eval();

        SelfAdjointEigenSolver<Matrix3d> eig(be_tr);
        if (eig.info() != Eigen::Success)
            throw std::runtime_error("Failed to diagonalize the elastic Finger tensor.");

        const Vector3d beta = eig.eigenvalues();
        if (!beta.allFinite() || (beta.array() <= 0.0).any())
            throw std::runtime_error("Elastic Finger tensor must be symmetric positive definite.");

        const Vector3d  log_beta    = beta.array().log().matrix();
        const double    tr_log_be   = log_beta.sum();
        const double    mean_log_be = tr_log_be / 3.0;
        const Vector3d  dev_log_be  = log_beta.array() - mean_log_be;
        const Vector3d  s_tr        = p.mu * dev_log_be;
        const double    tau_eq      = std::sqrt(1.5 * s_tr.squaredNorm());
        const double    ep_in       = ep_t(state);
        const double    sigma_y     = p.sigma0 + p.H * ep_in;
        const double    phi         = tau_eq - sigma_y;
        const double    tol         = 1e-12 * std::max({1.0, tau_eq, sigma_y});
        const double    tau_vol     = p.half_K * tr_log_be;
        const Matrix3d &Q           = eig.eigenvectors();

        F_inv[state] = Finv;

        if (phi <= tol) {
            const Vector3d tau_principal = Vector3d::Constant(tau_vol) + s_tr;
            Matrix3d       tau;
            tau.noalias() = Q * tau_principal.asDiagonal() * Q.transpose();

            be[state]        = be_tr;
            ep(state)        = ep_in;
            const Matrix3d P = tau * Finv.transpose();
            store_P(i, P);
            return;
        }

        const double   dgamma        = phi / p.D;
        const double   r             = 1.0 - p.three_mu * dgamma / tau_eq;
        const Vector3d log_be_new    = Vector3d::Constant(mean_log_be) + r * dev_log_be;
        const Vector3d tau_principal = Vector3d::Constant(tau_vol) + p.mu * r * dev_log_be;

        Matrix3d be_new;
        Matrix3d tau;
        be_new.noalias() = Q * log_be_new.array().exp().matrix().asDiagonal() * Q.transpose();
        be_new           = 0.5 * (be_new + be_new.transpose()).eval();
        tau.noalias()    = Q * tau_principal.asDiagonal() * Q.transpose();

        be[state]        = be_new;
        ep(state)        = ep_in + dgamma;
        const Matrix3d P = tau * Finv.transpose();
        store_P(i, P);
    }

  private:
    struct Params {
        double mu;
        double sigma0;
        double H;
        double half_K;
        double three_mu;
        double D;
    };

    using Matrix3dVector = std::vector<Matrix3d, Eigen::aligned_allocator<Matrix3d>>;

    vector<Params> params;

    Matrix3dVector F_inv;
    Matrix3dVector F_inv_t;
    Matrix3dVector be;
    Matrix3dVector be_t;
    VectorXd       ep;
    VectorXd       ep_t;
};

#endif // FINITESTRAINJ2PLASTICITY_H
