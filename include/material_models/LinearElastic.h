#ifndef LINEARELASTIC_H
#define LINEARELASTIC_H

#include "matmodel.h"
#include <Eigen/StdVector> // For Eigen's aligned_allocator

class LinearElasticIsotropic : public MechModel, public LinearModel<3> {
  public:
    LinearElasticIsotropic(vector<double> l_e, map<string, vector<double>> materialProperties)
        : MechModel(l_e)
    {
        try {
            bulk_modulus = materialProperties["bulk_modulus"];
            mu           = materialProperties["shear_modulus"];
        } catch (const std::exception &e) {
            throw std::runtime_error("Missing material properties for the requested material model.");
        }

        n_mat = bulk_modulus.size();
        lambda.resize(n_mat);
        mu.resize(n_mat);

        for (int i = 0; i < n_mat; ++i) {
            lambda[i] = bulk_modulus[i] - (2.0 / 3.0) * mu[i];
        }

        double lambda_ref = (*max_element(lambda.begin(), lambda.end()) +
                             *min_element(lambda.begin(), lambda.end())) /
                            2;
        double mu_ref = (*max_element(mu.begin(), mu.end()) +
                         *min_element(mu.begin(), mu.end())) /
                        2;

        kapparef_mat = Matrix<double, 6, 6>::Zero();
        kapparef_mat.topLeftCorner(3, 3).setConstant(lambda_ref);
        kapparef_mat += 2 * mu_ref * Matrix<double, 6, 6>::Identity();

        phase_stiffness = new Matrix<double, 24, 24>[n_mat];
        Matrix<double, 6, 6> phase_kappa;

        for (int i = 0; i < n_mat; i++) {
            phase_kappa.setZero();
            phase_stiffness[i] = Matrix<double, 24, 24>::Zero();

            phase_kappa.topLeftCorner(3, 3).setConstant(lambda[i]);
            phase_kappa += 2 * mu[i] * Matrix<double, 6, 6>::Identity();

            for (int p = 0; p < 8; ++p) {
                phase_stiffness[i] += B_int[p].transpose() * phase_kappa * B_int[p] * v_e * 0.1250;
            }
        }
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        double buf1     = lambda[mat_index] * (eps(i, 0) + eps(i + 1, 0) + eps(i + 2, 0));
        double buf2     = 2 * mu[mat_index];
        sigma(i + 0, 0) = buf1 + buf2 * eps(i + 0, 0);
        sigma(i + 1, 0) = buf1 + buf2 * eps(i + 1, 0);
        sigma(i + 2, 0) = buf1 + buf2 * eps(i + 2, 0);
        sigma(i + 3, 0) = buf2 * eps(i + 3, 0);
        sigma(i + 4, 0) = buf2 * eps(i + 4, 0);
        sigma(i + 5, 0) = buf2 * eps(i + 5, 0);
    }

  private:
    vector<double> bulk_modulus;
    vector<double> lambda;
    vector<double> mu;
};

class LinearElasticTriclinic : public MechModel, public LinearModel<3> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // Ensure proper alignment for Eigen structures

    LinearElasticTriclinic(const vector<double> &l_e, const map<string, vector<double>> &materialProperties)
        : MechModel(l_e)
    {
        vector<string> C_keys = {
            "C_11", "C_12", "C_13", "C_14", "C_15", "C_16",
            "C_22", "C_23", "C_24", "C_25", "C_26",
            "C_33", "C_34", "C_35", "C_36",
            "C_44", "C_45", "C_46",
            "C_55", "C_56",
            "C_66"};

        n_mat                = materialProperties.at("C_11").size();
        size_t num_constants = C_keys.size();

        // Initialize matrix to hold all constants
        MatrixXd C_constants(num_constants, n_mat);

        // Load material constants into matrix
        for (size_t k = 0; k < num_constants; ++k) {
            const auto &values = materialProperties.at(C_keys[k]);
            if (values.size() != n_mat) {
                throw std::runtime_error("Inconsistent size for material property: " + C_keys[k]);
            }
            C_constants.row(k) = Eigen::Map<const RowVectorXd>(values.data(), values.size());
        }

        // Assemble stiffness matrices for each material
        C_mats.resize(n_mat);
        kapparef_mat = Matrix<double, 6, 6>::Zero();

        for (size_t i = 0; i < n_mat; ++i) {
            Matrix<double, 6, 6> C_i = Matrix<double, 6, 6>::Zero();
            int                  k   = 0; // Index for C_constants

            // Assign constants to the upper triangular part
            for (int row = 0; row < 6; ++row) {
                for (int col = row; col < 6; ++col) {
                    C_i(row, col) = C_constants(k++, i);
                }
            }

            // Symmetrize the matrix
            C_i = C_i.selfadjointView<Eigen::Upper>();

            C_mats[i] = C_i;
            kapparef_mat += C_i;
        }

        kapparef_mat /= n_mat;

        // Compute phase stiffness matrices
        phase_stiffness = new Matrix<double, 24, 24>[n_mat];
        for (size_t i = 0; i < n_mat; ++i) {
            phase_stiffness[i] = Matrix<double, 24, 24>::Zero();
            for (int p = 0; p < 8; ++p) {
                phase_stiffness[i] += B_int[p].transpose() * C_mats[i] * B_int[p] * v_e * 0.1250;
            }
        }
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        sigma.segment<6>(i) = C_mats[mat_index] * eps.segment<6>(i);
    }

  private:
    std::vector<Matrix<double, 6, 6>, Eigen::aligned_allocator<Matrix<double, 6, 6>>> C_mats;
};

#endif // LINEARELASTIC_H
