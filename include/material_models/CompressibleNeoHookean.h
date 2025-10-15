#ifndef COMPRESSIBLENEOHOOKEAN_H
#define COMPRESSIBLENEOHOOKEAN_H

#include "LargeStrainMechModel.h"

/**
 * @brief Compressible Neo-Hookean hyperelastic material model
 *
 * Constitutive law: S = lambda log(J) C^{-1} + mu (I - C^{-1})
 * where:
 *   C = F^T F is the right Cauchy-Green tensor
 *   J = det(F) is the Jacobian
 */
class CompressibleNeoHookean : public LargeStrainMechModel {
  public:
    CompressibleNeoHookean(const Reader &reader)
        : LargeStrainMechModel(reader)
    {
        try {
            bulk_modulus  = reader.materialProperties["bulk_modulus"].get<vector<double>>();
            shear_modulus = reader.materialProperties["shear_modulus"].get<vector<double>>();
        } catch (json::exception &e) {
            throw std::runtime_error("Error reading CompressibleNeoHookean material properties: " + string(e.what()));
        }
        n_mat = bulk_modulus.size();
        lambda.resize(n_mat);
        mu.resize(n_mat);

        for (int i = 0; i < n_mat; ++i) {
            mu[i]     = shear_modulus[i];
            lambda[i] = bulk_modulus[i] - (2.0 / 3.0) * mu[i];
        }

        // Compute reference tangent at F=I
        kapparef_mat = Matrix<double, 9, 9>::Zero();
        for (int mat_idx = 0; mat_idx < n_mat; ++mat_idx) {
            Matrix3d             F_identity = Matrix3d::Identity();
            Matrix3d             S          = compute_S(F_identity, mat_idx, 0);
            Matrix<double, 6, 6> C_mandel   = compute_material_tangent(F_identity, mat_idx);
            Matrix<double, 9, 9> A          = compute_spatial_tangent(F_identity, S, C_mandel);
            kapparef_mat += A;
        }
        kapparef_mat /= static_cast<double>(n_mat);
    }

    Matrix3d compute_S(const Matrix3d &F, int mat_index, ptrdiff_t element_idx) override
    {
        Matrix3d C = compute_C(F);
        double   J = F.determinant();

        if (J <= 0.0) {
            throw std::runtime_error("Negative Jacobian determinant in CompressibleNeoHookean!");
        }

        double   logJ  = log(J);
        Matrix3d C_inv = C.inverse();

        return lambda[mat_index] * logJ * C_inv + mu[mat_index] * (Matrix3d::Identity() - C_inv);
    }

    Matrix<double, 6, 6> compute_material_tangent(const Matrix3d &F, int mat_index) override
    {
        Matrix3d C = compute_C(F);
        double   J = F.determinant();

        if (J <= 0.0) {
            throw std::runtime_error("Negative Jacobian determinant in CompressibleNeoHookean!");
        }

        double   logJ  = log(J);
        Matrix3d C_inv = C.inverse();

        Matrix<double, 6, 1> C_inv_mandel;
        C_inv_mandel << C_inv(0, 0), C_inv(1, 1), C_inv(2, 2),
            sqrt(2.0) * C_inv(0, 1), sqrt(2.0) * C_inv(0, 2), sqrt(2.0) * C_inv(1, 2);

        // PP1 = C_inv X C_inv (dyadic product in Mandel notation)
        Matrix<double, 6, 6> PP1 = C_inv_mandel * C_inv_mandel.transpose();

        // PP2 = symmetric product of C_inv in Mandel notation
        // PP2_IJKL = C_inv_IK C_inv_JL + C_inv_IL C_inv_JK
        Matrix<double, 6, 6> PP2                = Matrix<double, 6, 6>::Zero();
        int                  mandel_to_ij[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}};
        double               idx_fac[6]         = {1.0, 1.0, 1.0, sqrt(2.0), sqrt(2.0), sqrt(2.0)};

        for (int i = 0; i < 6; ++i) {
            int I = mandel_to_ij[i][0];
            int J = mandel_to_ij[i][1];
            for (int j = 0; j < 6; ++j) {
                int K = mandel_to_ij[j][0];
                int L = mandel_to_ij[j][1];

                PP2(i, j) = C_inv(I, K) * C_inv(J, L) + C_inv(I, L) * C_inv(J, K);
                PP2(i, j) *= idx_fac[i] * idx_fac[j];
            }
        }

        // Material tangent: C = lambda PP1 + (mu - lambda log(J)) PP2
        Matrix<double, 6, 6> C_tangent = lambda[mat_index] * PP1 +
                                         (mu[mat_index] - lambda[mat_index] * logJ) * PP2;

        return C_tangent;
    }

  private:
    vector<double> lambda;
    vector<double> mu;
    vector<double> bulk_modulus;
    vector<double> shear_modulus;
};

#endif // COMPRESSIBLENEOHOOKEAN_H
