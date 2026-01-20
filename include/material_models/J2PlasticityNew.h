#ifndef J2PLASTICITYNEW_H
#define J2PLASTICITYNEW_H

#include "matmodel.h"
#include "solver.h"

class J2PlasticityNew : public SmallStrainMechModel {
  public:
    J2PlasticityNew(const Reader &reader)
        : SmallStrainMechModel(reader)
    {
        try {
            bulk_modulus  = reader.materialProperties["bulk_modulus"].get<vector<double>>();
            shear_modulus = reader.materialProperties["shear_modulus"].get<vector<double>>();
            yield_stress  = reader.materialProperties["yield_stress"].get<vector<double>>();
            K             = reader.materialProperties["isotropic_hardening_parameter"].get<vector<double>>();
        } catch (const std::out_of_range &e) {
            throw std::runtime_error("Missing material properties for J2PlasticityNew.");
        }
        n_mat = bulk_modulus.size();

        c23 = sqrt(2.0 / 3.0);

        sig_t.setZero();
        n.setZero();
        eps_elastic.setZero();
    }

    virtual void initializeInternalVariables(ptrdiff_t num_elements, int num_gauss_points) override
    {
        plasticStrain.resize(num_elements, Matrix<double, 6, Dynamic>::Zero(6, num_gauss_points));
        plasticStrain_t.resize(num_elements, Matrix<double, 6, Dynamic>::Zero(6, num_gauss_points));
        q.resize(num_elements, VectorXd::Zero(num_gauss_points));
        q_t.resize(num_elements, VectorXd::Zero(num_gauss_points));
    }

    virtual void updateInternalVariables() override
    {
        plasticStrain_t = plasticStrain;
        q_t             = q;
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        // c23 = np.sqrt(2./3.)
        // (Already initialized in constructor)

        // s = 2.*mat_param["G_mod"]*(e-ep_in)
        double               G     = shear_modulus[mat_index];
        Matrix<double, 6, 1> e     = eps.block<6, 1>(i, 0);
        Matrix<double, 6, 1> ep_in = plasticStrain_t[element_idx].col(i / 6);
        Matrix<double, 6, 1> s     = 2.0 * G * (e - ep_in);

        // lam = (mat_param["K_mod"]-2./3.*mat_param["G_mod"])
        double lam = bulk_modulus[mat_index] - 2.0 / 3.0 * G;

        // s[:3] += lam*e[:3].sum()
        s.head<3>().array() += lam * e.head<3>().sum();

        // sig_t = np.array(s)
        sig_t = s;

        // sig_t[:3] -= s[:3].mean()
        sig_t.head<3>().array() -= s.head<3>().mean();

        // s_t = np.linalg.norm(sig_t)
        double s_t = sig_t.norm();

        // q = float(q_in)
        double q_in  = q_t[element_idx](i / 6);
        double q_var = q_in;

        // dq = 0.
        double dq = 0.0;

        // sy, dsy = YieldStress(q, dq, dt, mat_param)
        // YieldStress_LinHard: return mat_param["s_y"] + q*mat_param["h_iso"], mat_param["h_iso"]
        double sy  = yield_stress[mat_index] + q_var * K[mat_index];
        double dsy = K[mat_index];

        // phi = s_t - c23*sy
        double phi = s_t - c23 * sy;

        // dgamma = np.max((0., phi/(2.*mat_param["G_mod"]+2./3.*dsy)))
        double dgamma = std::max(0.0, phi / (2.0 * G + 2.0 / 3.0 * dsy));

        // ep = ep_in + dgamma*sig_t/s_t
        Matrix<double, 6, 1> ep;
        if (s_t > 1e-12) {
            ep = ep_in + dgamma * sig_t / s_t;
        } else {
            ep = ep_in;
        }

        // q += c23*dgamma
        q_var += c23 * dgamma;

        // s = s - dgamma * 2. * mat_param["G_mod"] * sig_t/s_t
        if (s_t > 1e-12) {
            s = s - dgamma * 2.0 * G * sig_t / s_t;
        }

        // Update internal variables
        plasticStrain[element_idx].col(i / 6) = ep;
        q[element_idx](i / 6)                 = q_var;

        // return s, ep, q
        sigma.block<6, 1>(i, 0) = s;
    }

    Matrix<double, 6, 6> get_reference_stiffness() const override
    {
        const double Kbar      = std::accumulate(bulk_modulus.begin(), bulk_modulus.end(), 0.0) / static_cast<double>(n_mat);
        const double Gbar      = std::accumulate(shear_modulus.begin(), shear_modulus.end(), 0.0) / static_cast<double>(n_mat);
        const double lambdabar = Kbar - 2.0 * Gbar / 3.0;

        Matrix<double, 6, 6> kappa_ref = Matrix<double, 6, 6>::Zero();
        kappa_ref.topLeftCorner(3, 3).setConstant(lambdabar);
        kappa_ref.diagonal().array() += 2.0 * Gbar;
        return kappa_ref;
    }

    void postprocess(Solver<3, 6> &solver, Reader &reader, int load_idx, int time_idx) override;

  protected:
    // Material properties
    vector<double> bulk_modulus;  // K_mod
    vector<double> shear_modulus; // G_mod (variable G in code)
    vector<double> yield_stress;  // s_y
    vector<double> K;             // h_iso (isotropic hardening modulus)

    // Internal variables
    vector<Matrix<double, 6, Dynamic>> plasticStrain;   // eps_p
    vector<Matrix<double, 6, Dynamic>> plasticStrain_t; // Previous time step
    vector<VectorXd>                   q;               // q (isotropic hardening variable)
    vector<VectorXd>                   q_t;             // Previous time step

    // Preallocated temporary variables
    Matrix<double, 6, 1> sig_t;
    Matrix<double, 6, 1> n;
    Matrix<double, 6, 1> eps_elastic;
    double               c23;
};

// Derived class: Linear Isotropic Hardening
class J2PlasticityNew_LinearIsotropicHardening : public J2PlasticityNew {
  public:
    J2PlasticityNew_LinearIsotropicHardening(const Reader &reader)
        : J2PlasticityNew(reader)
    {
    }
};

// Postprocess implementation
inline void J2PlasticityNew::postprocess(Solver<3, 6> &solver, Reader &reader, int load_idx, int time_idx)
{
    int n_str = 6;
    int n_gp  = plasticStrain_t[0].cols();

    auto &results                = reader.resultsToWrite;
    bool  need_plastic_strain    = std::find(results.begin(), results.end(), "plastic_strain") != results.end();
    bool  need_plastic_strain_gp = std::find(results.begin(), results.end(), "plastic_strain_gp") != results.end();
    bool  need_iso_hard          = std::find(results.begin(), results.end(), "isotropic_hardening_variable") != results.end();
    bool  need_iso_hard_gp       = std::find(results.begin(), results.end(), "isotropic_hardening_variable_gp") != results.end();

    VectorXd *plastic_strain_elem    = need_plastic_strain ? new VectorXd(solver.local_n0 * solver.n_y * solver.n_z * n_str) : nullptr;
    VectorXd *plastic_strain_gp_data = need_plastic_strain_gp ? new VectorXd(solver.local_n0 * solver.n_y * solver.n_z * n_gp * n_str) : nullptr;
    VectorXd *iso_hard_elem          = need_iso_hard ? new VectorXd(solver.local_n0 * solver.n_y * solver.n_z) : nullptr;
    VectorXd *iso_hard_gp_data       = need_iso_hard_gp ? new VectorXd(solver.local_n0 * solver.n_y * solver.n_z * n_gp) : nullptr;

    for (ptrdiff_t elem_idx = 0; elem_idx < solver.local_n0 * solver.n_y * solver.n_z; ++elem_idx) {
        if (need_plastic_strain) {
            (*plastic_strain_elem).segment(n_str * elem_idx, n_str) = plasticStrain_t[elem_idx].rowwise().mean();
        }
        if (need_iso_hard) {
            (*iso_hard_elem)(elem_idx) = q_t[elem_idx].mean();
        }

        if (need_plastic_strain_gp || need_iso_hard_gp) {
            for (int gp = 0; gp < n_gp; ++gp) {
                if (need_plastic_strain_gp) {
                    (*plastic_strain_gp_data).segment(n_str * n_gp * elem_idx + n_str * gp, n_str) = plasticStrain_t[elem_idx].col(gp);
                }
                if (need_iso_hard_gp) {
                    (*iso_hard_gp_data)(n_gp *elem_idx + gp) = q_t[elem_idx](gp);
                }
            }
        }
    }

    if (need_plastic_strain)
        reader.writeSlab("plastic_strain", load_idx, time_idx, plastic_strain_elem->data(), {n_str});
    if (need_plastic_strain_gp)
        reader.writeSlab("plastic_strain_gp", load_idx, time_idx, plastic_strain_gp_data->data(), {n_gp, n_str});
    if (need_iso_hard)
        reader.writeSlab("isotropic_hardening_variable", load_idx, time_idx, iso_hard_elem->data(), {1});
    if (need_iso_hard_gp)
        reader.writeSlab("isotropic_hardening_variable_gp", load_idx, time_idx, iso_hard_gp_data->data(), {n_gp});

    if (plastic_strain_elem)
        delete plastic_strain_elem;
    if (plastic_strain_gp_data)
        delete plastic_strain_gp_data;
    if (iso_hard_elem)
        delete iso_hard_elem;
    if (iso_hard_gp_data)
        delete iso_hard_gp_data;
}

#endif // J2PLASTICITYNEW_H
