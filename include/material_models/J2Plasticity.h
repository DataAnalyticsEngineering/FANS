#ifndef J2PLASTICITY_H
#define J2PLASTICITY_H

#include "matmodel.h"
#include "solver.h"

class J2Plasticity : public SmallStrainMechModel {
  public:
    J2Plasticity(const Reader &reader)
        : SmallStrainMechModel(reader)
    {
        try {
            bulk_modulus  = reader.materialProperties["bulk_modulus"].get<vector<double>>();
            shear_modulus = reader.materialProperties["shear_modulus"].get<vector<double>>();
            yield_stress  = reader.materialProperties["yield_stress"].get<vector<double>>();                  // Initial yield stress
            K             = reader.materialProperties["isotropic_hardening_parameter"].get<vector<double>>(); // Isotropic hardening parameter
            H             = reader.materialProperties["kinematic_hardening_parameter"].get<vector<double>>(); // Kinematic hardening parameter
            eta           = reader.materialProperties["viscosity"].get<vector<double>>();                     // Viscosity parameter
            dt            = reader.materialProperties["time_step"].get<double>();                             // Time step
        } catch (const std::out_of_range &e) {
            throw std::runtime_error("Missing material properties for the requested material model.");
        }
        n_mat = bulk_modulus.size();

        const double Kbar      = std::accumulate(bulk_modulus.begin(), bulk_modulus.end(), 0.0) / static_cast<double>(n_mat);
        const double Gbar      = std::accumulate(shear_modulus.begin(), shear_modulus.end(), 0.0) / static_cast<double>(n_mat);
        const double lambdabar = Kbar - 2.0 * Gbar / 3.0;
        kapparef_mat.setZero();
        kapparef_mat.topLeftCorner(3, 3).setConstant(lambdabar);
        kapparef_mat.diagonal().array() += 2.0 * Gbar;

        // Allocate the member matrices/vectors for performance optimization
        sqrt_two_over_three = sqrt(2.0 / 3.0);
        sigma_trial_n1.setZero();
        dev.setZero();
        qbar_trial_n1.setZero();
        n.setZero();
        eps_elastic.setZero();
        dev_minus_qbar.setZero();
    }

    /**
     * @brief Initializes internal variables for the J2 plasticity model.
     *
     * This function sets up the internal variables required for the J2 plasticity model.
     * It initializes the plastic strain and other internal variables for the given number
     * of elements and Gauss points.
     *
     * @param num_elements The number of elements in the model.
     * @param num_gauss_points The number of Gauss points per element.
     *
     * @note Variables with the suffix '_t' represent values from the previous time step.
     */
    virtual void initializeInternalVariables(ptrdiff_t num_elements, int num_gauss_points) override
    {
        // Initialize plastic strain and other internal variables
        plasticStrain.resize(num_elements, Matrix<double, 6, Dynamic>::Zero(6, num_gauss_points));
        plasticStrain_t.resize(num_elements, Matrix<double, 6, Dynamic>::Zero(6, num_gauss_points));
        psi.resize(num_elements, VectorXd::Zero(num_gauss_points));
        psi_t.resize(num_elements, VectorXd::Zero(num_gauss_points));
        psi_bar.resize(num_elements, Matrix<double, 6, Dynamic>::Zero(6, num_gauss_points));
        psi_bar_t.resize(num_elements, Matrix<double, 6, Dynamic>::Zero(6, num_gauss_points));
    }

    virtual void updateInternalVariables() override
    {
        plasticStrain_t = plasticStrain;
        psi_t           = psi;
        psi_bar_t       = psi_bar;
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        // Elastic Predictor
        eps_elastic = eps.block<6, 1>(i, 0) - plasticStrain_t[element_idx].col(i / 6);
        treps       = eps_elastic.head<3>().sum();

        // Compute trial stress
        sigma_trial_n1.head<3>().setConstant((bulk_modulus[mat_index] - 2.0 * shear_modulus[mat_index] / 3.0) * treps);
        sigma_trial_n1.head<3>() += 2 * shear_modulus[mat_index] * eps_elastic.head<3>();
        sigma_trial_n1.tail<3>() = 2 * shear_modulus[mat_index] * eps_elastic.tail<3>();

        // Deviatoric stress and norm
        dev = sigma_trial_n1;
        dev.head<3>().array() -= sigma_trial_n1.head<3>().mean();

        // Compute trial q and q_bar
        q_trial_n1    = compute_q_trial(psi_t[element_idx](i / 6), mat_index);
        qbar_trial_n1 = -(2.0 / 3.0) * H[mat_index] * psi_bar_t[element_idx].col(i / 6);

        // Calculate the trial yield function
        dev_minus_qbar      = dev - qbar_trial_n1;
        norm_dev_minus_qbar = dev_minus_qbar.norm();

        // Avoid division by zero
        if (norm_dev_minus_qbar < 1e-12) {
            n.setZero();
        } else {
            n = dev_minus_qbar / norm_dev_minus_qbar;
        }

        f_trial = norm_dev_minus_qbar - sqrt_two_over_three * (yield_stress[mat_index] - q_trial_n1);

        // Compute plastic multiplier
        gamma_n1 = (f_trial < 0) ? 0 : compute_gamma(f_trial, mat_index, i, element_idx);

        // Update stress and internal variables
        sigma_trial_n1 -= gamma_n1 * 2 * shear_modulus[mat_index] * n;
        plasticStrain[element_idx].col(i / 6) = plasticStrain_t[element_idx].col(i / 6) + gamma_n1 * n;
        psi[element_idx](i / 6) += gamma_n1 * sqrt_two_over_three;
        psi_bar[element_idx].col(i / 6) -= gamma_n1 * n;

        // Assign final stress
        sigma.block<6, 1>(i, 0) = sigma_trial_n1;
    }

    // Virtual methods for derived classes to implement different behaviors
    virtual double compute_q_trial(double psi_val, int mat_index)                             = 0;
    virtual double compute_gamma(double f_trial, int mat_index, int i, ptrdiff_t element_idx) = 0;

    void postprocess(Solver<3, 6> &solver, Reader &reader, const char *resultsFileName, int load_idx, int time_idx) override;

  protected:
    // Material properties
    vector<double> bulk_modulus;
    vector<double> shear_modulus;
    vector<double> yield_stress;
    vector<double> K;   // Isotropic hardening parameter
    vector<double> H;   // Kinematic hardening parameter
    vector<double> eta; // Viscosity parameter
    double         dt;  // Time step

    // Internal variables
    vector<Matrix<double, 6, Dynamic>> plasticStrain;
    vector<Matrix<double, 6, Dynamic>> plasticStrain_t;
    vector<VectorXd>                   psi;
    vector<VectorXd>                   psi_t;
    vector<Matrix<double, 6, Dynamic>> psi_bar;
    vector<Matrix<double, 6, Dynamic>> psi_bar_t;

    // Preallocated member variables for reuse
    Matrix<double, 6, 1> sigma_trial_n1;
    Matrix<double, 6, 1> dev;
    Matrix<double, 6, 1> qbar_trial_n1;
    Matrix<double, 6, 1> n;
    Matrix<double, 6, 1> eps_elastic;
    double               treps;
    double               q_trial_n1;
    double               f_trial;
    Matrix<double, 6, 1> dev_minus_qbar;
    double               norm_dev_minus_qbar;
    double               gamma_n1;
    double               sqrt_two_over_three;
};

// Derived Class Linear Isotropic Hardening
class J2ViscoPlastic_LinearIsotropicHardening : public J2Plasticity {
  public:
    J2ViscoPlastic_LinearIsotropicHardening(const Reader &reader)
        : J2Plasticity(reader)
    {
    }

    double compute_q_trial(double psi_val, int mat_index) override
    {
        return -K[mat_index] * psi_val;
    }

    double compute_gamma(double f_trial, int mat_index, int i, ptrdiff_t element_idx) override
    {
        return f_trial / (2 * shear_modulus[mat_index] + (2.0 / 3.0) * (K[mat_index] + H[mat_index]) + eta[mat_index] / dt);
    }
};

// Derived Class Non-Linear (Exponential law) Isotropic Hardening
class J2ViscoPlastic_NonLinearIsotropicHardening : public J2Plasticity {
  public:
    J2ViscoPlastic_NonLinearIsotropicHardening(const Reader &reader)
        : J2Plasticity(reader)
    {
        try {
            sigma_inf = reader.materialProperties["saturation_stress"].get<vector<double>>();   // Saturation stress
            delta     = reader.materialProperties["saturation_exponent"].get<vector<double>>(); // Saturation exponent
        } catch (const std::out_of_range &e) {
            throw std::runtime_error("Missing material properties for the requested material model.");
        }

        // Precompute constants for optimization
        denominator.resize(n_mat);
        sigma_diff.resize(n_mat);
        for (size_t i = 0; i < n_mat; ++i) {
            denominator[i] = 2 * shear_modulus[i] + (2.0 / 3.0) * (K[i] + H[i]) + eta[i] / dt;
            sigma_diff[i]  = sqrt_two_over_three * (sigma_inf[i] - yield_stress[i]);
        }
    }

    double compute_q_trial(double psi_val, int mat_index) override
    {
        return -K[mat_index] * psi_val - (sigma_inf[mat_index] - yield_stress[mat_index]) * (1 - exp(-delta[mat_index] * psi_val));
    }

    double compute_gamma(double f_trial, int mat_index, int i, ptrdiff_t element_idx) override
    {
        gamma_n1  = 0;
        gamma_inc = 1;
        NR_iter   = 0;
        while (gamma_inc > NR_tol && NR_iter < NR_max_iter) {
            g = f_trial - gamma_n1 * denominator[mat_index] -
                sigma_diff[mat_index] *
                    (-exp(-delta[mat_index] * (psi_t[element_idx](i / 6) + sqrt_two_over_three * gamma_n1)) + exp(-delta[mat_index] * psi_t[element_idx](i / 6)));
            dg = -denominator[mat_index] -
                 (2 / 3) * (sigma_inf[mat_index] - yield_stress[mat_index]) * delta[mat_index] * exp(-delta[mat_index] * (psi_t[element_idx](i / 6) + sqrt_two_over_three * gamma_n1));
            gamma_inc = -g / dg;
            gamma_n1 += gamma_inc;
            NR_iter++;
        }
        return gamma_n1;
    }

  protected:
    // Material properties
    vector<double> sigma_inf; // Saturation stress
    vector<double> delta;     // Saturation exponent
  private:
    // Newton-Raphson parameters
    double NR_tol      = 1e-10;
    int    NR_max_iter = 10;
    int    NR_iter;

    // Preallocated member variables for reuse
    double gamma_inc;
    double g;
    double dg;

    // Precomputed constants
    vector<double> denominator; // 2 * mu + H * (2/3) + eta/dt
    vector<double> sigma_diff;  // sqrt(2/3) * (sigma_inf - yield_stress)
};

inline void J2Plasticity::postprocess(Solver<3, 6> &solver, Reader &reader, const char *resultsFileName, int load_idx, int time_idx)
{
    int      n_str                             = 6; // The plastic strain and stress vectors have 6 components each
    VectorXd mean_plastic_strain               = VectorXd::Zero(solver.local_n0 * solver.n_y * solver.n_z * n_str);
    VectorXd mean_isotropic_hardening_variable = VectorXd::Zero(solver.local_n0 * solver.n_y * solver.n_z);
    VectorXd mean_kinematic_hardening_variable = VectorXd::Zero(solver.local_n0 * solver.n_y * solver.n_z * n_str);

    // Compute the mean values for each element
    for (ptrdiff_t elem_idx = 0; elem_idx < solver.local_n0 * solver.n_y * solver.n_z; ++elem_idx) {
        mean_plastic_strain.segment(n_str * elem_idx, n_str)               = plasticStrain_t[elem_idx].rowwise().mean();
        mean_isotropic_hardening_variable(elem_idx)                        = psi_t[elem_idx].mean();
        mean_kinematic_hardening_variable.segment(n_str * elem_idx, n_str) = psi_bar_t[elem_idx].rowwise().mean();
    }

    reader.writeSlab("plastic_strain", resultsFileName, solver.dataset_name, load_idx, time_idx, mean_plastic_strain.data(), n_str);
    reader.writeSlab("isotropic_hardening_variable", resultsFileName, solver.dataset_name, load_idx, time_idx, mean_isotropic_hardening_variable.data(), 1);
    reader.writeSlab("kinematic_hardening_variable", resultsFileName, solver.dataset_name, load_idx, time_idx, mean_kinematic_hardening_variable.data(), n_str);
}

#endif // J2PLASTICITY_H
