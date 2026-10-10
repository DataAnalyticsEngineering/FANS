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
    }

    Matrix3d compute_S(const Matrix3d &F, int mat_index, ptrdiff_t element_idx, int i) override
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

    bool has_tangent() const override
    {
        return true;
    }

    // A = dP/dF of P = mu F + (lambda log(J) - mu) G with G = F^{-T}:
    //   A_iJkL = mu d_ik d_JL + lambda G_iJ G_kL + (mu - lambda log(J)) G_iL G_kJ
    void get_tangent(int i_eps, int mat_index, ptrdiff_t element_idx, Tangent A) override
    {
        const Matrix3d F = extract_F(i_eps);
        const Matrix3d G = F.inverse().transpose();
        const double   c = mu[mat_index] - lambda[mat_index] * log(F.determinant());

        for (int i = 0; i < 3; ++i)
            for (int J = 0; J < 3; ++J)
                for (int k = 0; k < 3; ++k)
                    for (int L = 0; L < 3; ++L)
                        A(3 * i + J, 3 * k + L) = mu[mat_index] * (i == k) * (J == L) + lambda[mat_index] * G(i, J) * G(k, L) + c * G(i, L) * G(k, J);
    }

    // At F = I the tangent is that of linear isotropic elasticity
    Matrix<double, 9, 9> get_reference_stiffness() override
    {
        Matrix<double, 9, 9> kapparef = Matrix<double, 9, 9>::Zero();
        for (int mat_idx = 0; mat_idx < n_mat; ++mat_idx)
            kapparef += isotropic_reference_stiffness(lambda[mat_idx], mu[mat_idx]);
        return kapparef / static_cast<double>(n_mat);
    }

  private:
    vector<double> lambda;
    vector<double> mu;
    vector<double> bulk_modulus;
    vector<double> shear_modulus;
};

#endif // COMPRESSIBLENEOHOOKEAN_H
