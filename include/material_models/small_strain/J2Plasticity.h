#ifndef J2PLASTICITY_H
#define J2PLASTICITY_H

#include "matmodel.h"
#include "solver.h"

class J2Plasticity : public SmallStrainMechModel {
  public:
    struct YieldStressResponse {
        double value;
        double derivative;
    };

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
        } catch (const std::exception &e) {
            throw std::runtime_error("Missing or invalid material properties for J2Plasticity.");
        }
        n_mat = bulk_modulus.size();
        if (n_mat == 0 || shear_modulus.size() != n_mat || yield_stress.size() != n_mat ||
            K.size() != n_mat || H.size() != n_mat || eta.size() != n_mat)
            throw std::runtime_error("Inconsistent J2Plasticity material-property sizes.");

        if (dt <= 0.0)
            throw std::runtime_error("J2Plasticity requires time_step > 0.");

        lambda.resize(n_mat);
        for (int m = 0; m < n_mat; ++m)
            lambda[m] = bulk_modulus[m] - two_thirds * shear_modulus[m];
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
    void initializeInternalVariables(ptrdiff_t num_elements, int num_gauss_points) override
    {
        const ptrdiff_t num_states = num_elements * num_gauss_points;
        // Initialize plastic strain and other internal variables
        plasticStrain.setZero(6, num_states);
        plasticStrain_t.setZero(6, num_states);
        q.setZero(num_states);
        q_t.setZero(num_states);
        alpha.setZero(6, num_states);
        alpha_t.setZero(6, num_states);
    }

    void updateInternalVariables() override
    {
        plasticStrain.swap(plasticStrain_t);
        q.swap(q_t);
        alpha.swap(alpha_t);
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        const int                  gp        = i / 6;
        const ptrdiff_t            state_idx = element_idx * n_gp + gp;
        const double               G         = shear_modulus[mat_index];
        const Matrix<double, 6, 1> e         = eps.segment<6>(i);
        const Matrix<double, 6, 1> ep_in     = plasticStrain_t.col(state_idx);
        const Matrix<double, 6, 1> alpha_in  = alpha_t.col(state_idx);
        const double               q_in      = q_t(state_idx);

        // Elastic trial stress
        const Matrix<double, 6, 1> ee_t = e - ep_in;
        Matrix<double, 6, 1>       s    = 2.0 * G * ee_t;
        s.head<3>().array() += lambda[mat_index] * ee_t.head<3>().sum();

        // Shifted trial deviatoric stress
        Matrix<double, 6, 1> xi_t = s - alpha_in;
        xi_t.head<3>().array() -= s.head<3>().mean();

        const YieldStressResponse yield       = compute_yield_stress(q_in, mat_index);
        const double              tolerance   = 1e-12 * std::max(1.0, yield.value);
        const double              yield_limit = c23 * yield.value + tolerance;
        const double              xi_norm_sq  = xi_t.squaredNorm();

        // Elastic step
        if (xi_norm_sq <= yield_limit * yield_limit) {
            plasticStrain.col(state_idx) = ep_in;
            q(state_idx)                 = q_in;
            alpha.col(state_idx)         = alpha_in;
            sigma.segment<6>(i)          = s;
            return;
        }

        // Plastic step
        const double               xi_norm = std::sqrt(xi_norm_sq);
        const double               phi     = xi_norm - c23 * yield.value; // Trial Von-Mises yield function: phi = ||xi_trial|| - sqrt(2/3) sigma_y
        const Matrix<double, 6, 1> n       = xi_t / xi_norm;
        const double               dgamma  = compute_gamma(phi, q_in, yield, mat_index);

        s -= 2.0 * G * dgamma * n;

        plasticStrain.col(state_idx) = ep_in + dgamma * n;
        q(state_idx)                 = q_in + c23 * dgamma;
        alpha.col(state_idx)         = alpha_in + two_thirds * H[mat_index] * dgamma * n;

        sigma.segment<6>(i) = s;
    }

    virtual YieldStressResponse compute_yield_stress(double q, int mat_index) const = 0;

    double compute_gamma(double phi, double q_in, const YieldStressResponse &yield_response_in, int mat_index) const
    {
        const double D                   = 2.0 * shear_modulus[mat_index] + two_thirds * H[mat_index] + eta[mat_index] / dt;
        const double initial_denominator = D + two_thirds * yield_response_in.derivative;
        const double residual_tolerance  = 1e-12 * std::max(1.0, std::max(std::abs(phi), std::abs(yield_response_in.value)));
        const int    max_iterations      = 50;

        double dgamma = phi / initial_denominator;
        int    iter   = 0;

        while (iter < max_iterations) {
            const double              q_new              = q_in + c23 * dgamma;
            const YieldStressResponse yield_response_new = compute_yield_stress(q_new, mat_index);
            const double              residual           = phi - D * dgamma - c23 * (yield_response_new.value - yield_response_in.value);

            if (std::abs(residual) <= residual_tolerance)
                return dgamma;

            const double derivative = -D - two_thirds * yield_response_new.derivative;
            dgamma -= residual / derivative;
            ++iter;
        }

        throw std::runtime_error("J2Plasticity return mapping did not converge.");
    }

    Matrix<double, 6, 6> get_reference_stiffness() override
    {
        const double Kbar      = std::accumulate(bulk_modulus.begin(), bulk_modulus.end(), 0.0) / static_cast<double>(n_mat);
        const double Gbar      = std::accumulate(shear_modulus.begin(), shear_modulus.end(), 0.0) / static_cast<double>(n_mat);
        const double lambdabar = Kbar - two_thirds * Gbar;

        Matrix<double, 6, 6> kappa_ref = Matrix<double, 6, 6>::Zero();
        kappa_ref.topLeftCorner(3, 3).setConstant(lambdabar);
        kappa_ref.diagonal().array() += 2.0 * Gbar;
        return kappa_ref;
    }

    void postprocess(Solver<3, 6> &solver, Reader &reader, int load_idx, int time_idx) override;

  protected:
    // Material properties
    vector<double> bulk_modulus;
    vector<double> shear_modulus;
    vector<double> yield_stress;
    vector<double> K;      // Isotropic hardening parameter
    vector<double> H;      // Kinematic hardening parameter
    vector<double> eta;    // Viscosity parameter
    vector<double> lambda; // First Lame parameter
    double         dt;     // Time step

    // Internal variables
    Matrix<double, 6, Dynamic> plasticStrain;
    Matrix<double, 6, Dynamic> plasticStrain_t;
    VectorXd                   q;
    VectorXd                   q_t;
    Matrix<double, 6, Dynamic> alpha;
    Matrix<double, 6, Dynamic> alpha_t;

    static constexpr double c23        = 0.81649658092772603273;
    static constexpr double two_thirds = 2.0 / 3.0;
};

// Derived Class Linear Isotropic Hardening
class J2ViscoPlastic_LinearIsotropicHardening : public J2Plasticity {
  public:
    J2ViscoPlastic_LinearIsotropicHardening(const Reader &reader)
        : J2Plasticity(reader)
    {
    }

    YieldStressResponse compute_yield_stress(double q, int mat_index) const override
    {
        return {yield_stress[mat_index] + K[mat_index] * q, K[mat_index]};
    }
};

// Derived Class Non-Linear (Exponential law) Isotropic Hardening
class J2ViscoPlastic_NonLinearIsotropicHardening : public J2Plasticity {
  public:
    J2ViscoPlastic_NonLinearIsotropicHardening(const Reader &reader)
        : J2Plasticity(reader)
    {
        try {
            sigma_inf = reader.materialProperties["saturation_stress"].get<vector<double>>();
            delta     = reader.materialProperties["saturation_exponent"].get<vector<double>>();
        } catch (const std::exception &e) {
            throw std::runtime_error("Missing or invalid nonlinear J2Plasticity properties.");
        }

        if (sigma_inf.size() != static_cast<size_t>(n_mat) ||
            delta.size() != static_cast<size_t>(n_mat))
            throw std::runtime_error("Nonlinear J2Plasticity property vectors have invalid size.");
    }

    YieldStressResponse compute_yield_stress(double q, int mat_index) const override
    {
        const double A           = sigma_inf[mat_index] - yield_stress[mat_index];
        const double exponential = std::exp(-delta[mat_index] * q);

        return {yield_stress[mat_index] + K[mat_index] * q + A * (1.0 - exponential),
                K[mat_index] + A * delta[mat_index] * exponential};
    }

  protected:
    // Material properties
    vector<double> sigma_inf; // Saturation stress
    vector<double> delta;     // Saturation exponent
};

inline void J2Plasticity::postprocess(Solver<3, 6> &solver, Reader &reader, int load_idx, int time_idx)
{
    const int       n_str  = 6;
    const ptrdiff_t n_elem = solver.local_n0 * solver.n_y * solver.n_z;

    // Check what user requested
    auto      &results                = reader.resultsToWrite;
    const bool need_plastic_strain    = std::find(results.begin(), results.end(), "plastic_strain") != results.end();
    const bool need_plastic_strain_gp = std::find(results.begin(), results.end(), "plastic_strain_gp") != results.end();
    const bool need_iso_hard          = std::find(results.begin(), results.end(), "isotropic_hardening_variable") != results.end();
    const bool need_iso_hard_gp       = std::find(results.begin(), results.end(), "isotropic_hardening_variable_gp") != results.end();
    const bool need_kin_hard          = std::find(results.begin(), results.end(), "kinematic_hardening_variable") != results.end();
    const bool need_kin_hard_gp       = std::find(results.begin(), results.end(), "kinematic_hardening_variable_gp") != results.end();

    VectorXd plastic_strain_elem, plastic_strain_gp_data;
    VectorXd iso_hard_elem, iso_hard_gp_data;
    VectorXd kin_hard_elem, kin_hard_gp_data;

    if (need_plastic_strain)
        plastic_strain_elem.resize(n_elem * n_str);
    if (need_plastic_strain_gp)
        plastic_strain_gp_data.resize(n_elem * n_gp * n_str);
    if (need_iso_hard)
        iso_hard_elem.resize(n_elem);
    if (need_iso_hard_gp)
        iso_hard_gp_data.resize(n_elem * n_gp);
    if (need_kin_hard)
        kin_hard_elem.resize(n_elem * n_str);
    if (need_kin_hard_gp)
        kin_hard_gp_data.resize(n_elem * n_gp * n_str);

    for (ptrdiff_t elem_idx = 0; elem_idx < n_elem; ++elem_idx) {
        const ptrdiff_t state_offset = elem_idx * n_gp;
        // Element averages
        if (need_plastic_strain)
            plastic_strain_elem.segment(n_str * elem_idx, n_str) = plasticStrain_t.middleCols(state_offset, n_gp).rowwise().mean();

        if (need_iso_hard)
            iso_hard_elem(elem_idx) = q_t.segment(state_offset, n_gp).mean();

        if (need_kin_hard)
            kin_hard_elem.segment(n_str * elem_idx, n_str) = alpha_t.middleCols(state_offset, n_gp).rowwise().mean();

        // All Gauss point data
        for (int gp = 0; gp < n_gp; ++gp) {
            const ptrdiff_t state_idx = state_offset + gp;

            if (need_plastic_strain_gp)
                plastic_strain_gp_data.segment(n_str * n_gp * elem_idx + n_str * gp, n_str) = plasticStrain_t.col(state_idx);

            if (need_iso_hard_gp)
                iso_hard_gp_data(n_gp * elem_idx + gp) = q_t(state_idx);

            if (need_kin_hard_gp)
                kin_hard_gp_data.segment(n_str * n_gp * elem_idx + n_str * gp, n_str) = alpha_t.col(state_idx);
        }
    }

    // Write only what was requested
    if (need_plastic_strain)
        reader.writeSlab("plastic_strain", load_idx, time_idx, plastic_strain_elem.data(), {n_str});
    if (need_plastic_strain_gp)
        reader.writeSlab("plastic_strain_gp", load_idx, time_idx, plastic_strain_gp_data.data(), {n_gp, n_str});
    if (need_iso_hard)
        reader.writeSlab("isotropic_hardening_variable", load_idx, time_idx, iso_hard_elem.data(), {1});
    if (need_iso_hard_gp)
        reader.writeSlab("isotropic_hardening_variable_gp", load_idx, time_idx, iso_hard_gp_data.data(), {n_gp});
    if (need_kin_hard)
        reader.writeSlab("kinematic_hardening_variable", load_idx, time_idx, kin_hard_elem.data(), {n_str});
    if (need_kin_hard_gp)
        reader.writeSlab("kinematic_hardening_variable_gp", load_idx, time_idx, kin_hard_gp_data.data(), {n_gp, n_str});
}

#endif // J2PLASTICITY_H
