#ifndef LINEARTHERMAL_H
#define LINEARTHERMAL_H

#include "matmodel.h"
#include <Eigen/StdVector> // For Eigen's aligned_allocator

class LinearThermalIsotropic : public ThermalModel, public LinearModel<1, 3> {
  public:
    LinearThermalIsotropic(const Reader &reader)
        : ThermalModel(reader)
    {
        try {
            conductivity = reader.materialProperties["conductivity"].get<vector<double>>();
        } catch (const std::exception &e) {
            throw std::runtime_error("Missing material properties for the requested material model.");
        }
        n_mat = conductivity.size();

        double kappa_ref = (*max_element(conductivity.begin(), conductivity.end()) +
                            *min_element(conductivity.begin(), conductivity.end())) /
                           2;
        kapparef_mat = kappa_ref * Matrix3d::Identity();

        Matrix3d phase_kappa;
        phase_stiffness = new Matrix<double, 8, 8>[n_mat];

        for (size_t i = 0; i < n_mat; ++i) {
            phase_stiffness[i] = Matrix<double, 8, 8>::Zero();
            phase_kappa        = conductivity[i] * Matrix3d::Identity();

            for (int p = 0; p < n_gp; ++p) {
                phase_stiffness[i] += B_int[p].transpose() * phase_kappa * B_int[p] * v_e / n_gp;
            }
        }
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        sigma(i + 0, 0) = conductivity[mat_index] * eps(i + 0, 0);
        sigma(i + 1, 0) = conductivity[mat_index] * eps(i + 1, 0);
        sigma(i + 2, 0) = conductivity[mat_index] * eps(i + 2, 0);
    }

  private:
    vector<double> conductivity;
};

class LinearThermalTriclinic : public ThermalModel, public LinearModel<1, 3> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // Ensure proper alignment for Eigen structures

    LinearThermalTriclinic(const Reader &reader)
        : ThermalModel(reader)
    {
        vector<string> K_keys = {
            "K_11", "K_12", "K_13",
            "K_22", "K_23",
            "K_33"};

        try {
            n_mat                = reader.materialProperties.at("K_11").get<vector<double>>().size();
            size_t num_constants = K_keys.size();

            // Initialize matrix to hold all constants
            K_constants.resize(num_constants, n_mat);

            // Load material constants into matrix
            for (size_t k = 0; k < num_constants; ++k) {
                const auto &values = reader.materialProperties.at(K_keys[k]).get<vector<double>>();
                if (values.size() != n_mat) {
                    throw std::runtime_error("Inconsistent size for material property: " + K_keys[k]);
                }
                K_constants.row(k) = Eigen::Map<const RowVectorXd>(values.data(), values.size());
            }
        } catch (const std::exception &e) {
            throw std::runtime_error("Missing or inconsistent material properties for the requested material model.");
        }

        // Assemble conductivity matrices for each material
        K_mats.resize(n_mat);
        kapparef_mat = Matrix3d::Zero();

        for (size_t i = 0; i < n_mat; ++i) {
            Matrix3d K_i;
            K_i << K_constants(0, i), K_constants(1, i), K_constants(2, i),
                K_constants(1, i), K_constants(3, i), K_constants(4, i),
                K_constants(2, i), K_constants(4, i), K_constants(5, i);

            K_mats[i] = K_i;
            kapparef_mat += K_i;
        }

        kapparef_mat /= n_mat;

        // Compute phase stiffness matrices
        phase_stiffness = new Matrix<double, 8, 8>[n_mat];
        for (size_t i = 0; i < n_mat; ++i) {
            phase_stiffness[i] = Matrix<double, 8, 8>::Zero();
            for (int p = 0; p < n_gp; ++p) {
                phase_stiffness[i] += B_int[p].transpose() * K_mats[i] * B_int[p] * v_e / n_gp;
            }
        }
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        sigma.segment<3>(i) = K_mats[mat_index] * eps.segment<3>(i);
    }

  private:
    std::vector<Matrix3d, Eigen::aligned_allocator<Matrix3d>> K_mats;
    MatrixXd                                                  K_constants;
};

#endif // LINEARTHERMAL_H
