#ifndef SAINTVENANTKIRCHHOFF_H
#define SAINTVENANTKIRCHHOFF_H

#include "LargeStrainMechModel.h"

/**
 * @brief Saint-Venant Kirchhoff hyperelastic material model
 *
 * Constitutive law: S = lambda tr(E) I + 2mu E
 * where E = 0.5 * (F^T F - I) is the Green-Lagrange strain tensor
 */
class SaintVenantKirchhoff : public LargeStrainMechModel {
  public:
    SaintVenantKirchhoff(const Reader &reader)
        : LargeStrainMechModel(reader)
    {
        try {
            bulk_modulus  = reader.materialProperties["bulk_modulus"].get<vector<double>>();
            shear_modulus = reader.materialProperties["shear_modulus"].get<vector<double>>();
        } catch (json::exception &e) {
            throw std::runtime_error("Error reading SaintVenantKirchhoff material properties: " + string(e.what()));
        }
        n_mat = bulk_modulus.size();
        lambda.resize(n_mat);
        mu.resize(n_mat);

        for (int i = 0; i < n_mat; ++i) {
            mu[i]     = shear_modulus[i];
            lambda[i] = bulk_modulus[i] - (2.0 / 3.0) * mu[i];
        }
    }

    Matrix3d compute_S(const Matrix3d &F, int mat_index, ptrdiff_t element_idx, int i) override
    {
        Matrix3d E   = compute_E_from_F(F);
        double   trE = E.trace();
        return lambda[mat_index] * trE * Matrix3d::Identity() + 2.0 * mu[mat_index] * E;
    }

    Matrix<double, 6, 6> compute_material_tangent(const Matrix3d &F, int mat_index, ptrdiff_t element_idx, int i) override
    {
        Matrix<double, 6, 6> C = Matrix<double, 6, 6>::Zero();
        C.topLeftCorner(3, 3).setConstant(lambda[mat_index]);
        C += 2.0 * mu[mat_index] * Matrix<double, 6, 6>::Identity();
        return C;
    }

    Matrix<double, 9, 9> get_reference_stiffness() override
    {
        // Compute reference tangent at F=I
        Matrix<double, 9, 9> kapparef   = Matrix<double, 9, 9>::Zero();
        const Matrix3d       F_identity = Matrix3d::Identity();
        for (int mat_idx = 0; mat_idx < n_mat; ++mat_idx) {
            const Matrix3d             S = compute_S(F_identity, mat_idx, 0, 0);
            const Matrix<double, 6, 6> C = compute_material_tangent(F_identity, mat_idx, 0, 0);
            kapparef += compute_spatial_tangent(F_identity, S, C);
        }
        return kapparef / static_cast<double>(n_mat);
    }

  private:
    vector<double> lambda;
    vector<double> mu;
    vector<double> bulk_modulus;
    vector<double> shear_modulus;
};

#endif // SAINTVENANTKIRCHHOFF_H
