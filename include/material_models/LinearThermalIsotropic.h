#ifndef LINEARTHERMALISOTROPIC_H
#define LINEARTHERMALISOTROPIC_H

#include "matmodel.h"

class LinearThermalIsotropic : public ThermalModel, public LinearModel<1> {
  public:
    LinearThermalIsotropic(vector<double> l_e, json materialProperties)
        : ThermalModel(l_e)
    {
        try {
            conductivity = materialProperties["conductivity"].get<vector<double>>();
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

            for (int p = 0; p < 8; ++p) {
                phase_stiffness[i] += B_int[p].transpose() * phase_kappa * B_int[p] * v_e * 0.1250;
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

#endif // LINEARTHERMALISOTROPIC_H
