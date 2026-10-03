#ifndef LARGESTRAINMECHMODEL_H
#define LARGESTRAINMECHMODEL_H

#include "matmodel.h"

/**
 * @brief Base class for large strain mechanical material models
 * This class provides the kinematic framework for finite deformation mechanics.
 * Derived classes may implement the constitutive response through either compute_S() or get_sigma().
 */
class LargeStrainMechModel : public Matmodel<3, 9> {
  public:
    LargeStrainMechModel(const Reader &reader)
        : Matmodel(reader)
    {
        Construct_B();
    }

  protected:
    // ============================================================================
    // Kinematic helper functions
    // ============================================================================

    /**
     * @brief Extract deformation gradient from eps vector at Gauss point
     * @param i Index in eps array (should be n_str * gauss_point)
     * @return 3×3 deformation gradient matrix F
     */
    inline Matrix3d extract_F(int i) const
    {
        Matrix3d F;
        F << eps(i), eps(i + 1), eps(i + 2),
            eps(i + 3), eps(i + 4), eps(i + 5),
            eps(i + 6), eps(i + 7), eps(i + 8);
        return F;
    }
    /**
     * @brief Store 1st Piola-Kirchhoff stress in sigma array
     * @param i Index in sigma array (should be n_str * gauss_point)
     * @param P 3×3 1st Piola-Kirchhoff stress tensor
     */
    inline void store_P(int i, const Matrix3d &P)
    {
        sigma(i)     = P(0, 0);
        sigma(i + 1) = P(0, 1);
        sigma(i + 2) = P(0, 2);
        sigma(i + 3) = P(1, 0);
        sigma(i + 4) = P(1, 1);
        sigma(i + 5) = P(1, 2);
        sigma(i + 6) = P(2, 0);
        sigma(i + 7) = P(2, 1);
        sigma(i + 8) = P(2, 2);
    }

    inline Matrix3d compute_C(const Matrix3d &F) const
    {
        return F.transpose() * F; // C = F^T F
    }
    inline Matrix3d compute_E(const Matrix3d &C) const
    {
        return 0.5 * (C - Matrix3d::Identity()); // E = 0.5 * (C - I)
    }
    inline Matrix3d compute_E_from_F(const Matrix3d &F) const
    {
        return compute_E(compute_C(F)); // E = 0.5 * (F^T F - I)
    }
    inline Matrix3d push_forward(const Matrix3d &F, const Matrix3d &S) const
    {
        return F * S; // P = F S
    }

    // ============================================================================
    // Virtual functions for derived material models to implement
    // ============================================================================

    /**
     * @brief Convenience constitutive interface for models formulated in terms of S
     * @param F Deformation gradient
     * @param mat_index Material-property index local to this material model
     * @param element_idx Local element index
     * @param i Offset of this Gauss point in the flattened strain vector
     */
    virtual Matrix3d compute_S(const Matrix3d &F, int mat_index, ptrdiff_t element_idx, int i)
    {
        throw std::logic_error("Large-strain material model must override compute_S() or get_sigma().");
    }

    /**
     * @brief Reference stiffness dP/dF at F = I of an isotropic material with Lame
     * constants lambda and mu, as a 9x9 matrix on row-major F and P:
     *   A_iJkL = lambda d_iJ d_kL + mu (d_ik d_JL + d_iL d_Jk)
     */
    static Matrix<double, 9, 9> isotropic_reference_stiffness(double lambda, double mu)
    {
        Matrix<double, 9, 9> A = Matrix<double, 9, 9>::Zero();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                A(3 * i + i, 3 * j + j) += lambda; // lambda d_iJ d_kL
                A(3 * i + j, 3 * i + j) += mu;     // mu d_ik d_JL
                A(3 * i + j, 3 * j + i) += mu;     // mu d_iL d_Jk
            }
        return A;
    }

    // ============================================================================
    // Standard interface implementation
    // ============================================================================

    /**
     * @brief Compute B matrix for deformation gradient (9×24)
     * Returns B_F matrix that relates nodal displacements to F components:
     */
    Matrix<double, 9, 24> Compute_B(const double x, const double y, const double z) override;

    /**
     * @brief Default stress update for models formulated in terms of S
     * This function:
     * 1. Extracts F from eps
     * 2. Calls compute_S(F) - material model implements this
     * 3. Computes P = F S
     * 4. Stores P in sigma
     */
    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        Matrix3d F = extract_F(i);                            // Extract deformation gradient
        Matrix3d S = compute_S(F, mat_index, element_idx, i); // Material model computes 2nd PK stress from F
        Matrix3d P = push_forward(F, S);                      // Push forward to 1st PK stress
        store_P(i, P);                                        // Store in sigma array
    }
};

inline Matrix<double, 9, 24> LargeStrainMechModel::Compute_B(const double x, const double y, const double z)
{
    Matrix<double, 9, 24> B_F = Matrix<double, 9, 24>::Zero();
    Matrix<double, 3, 8>  dN  = Matmodel<3, 9>::Compute_basic_B(x, y, z);

    for (int i = 0; i < 3; ++i) {     // Displacement component
        for (int J = 0; J < 3; ++J) { // Material coordinate derivative
            int row = 3 * i + J;      // Row-major: F_iJ
            for (int node = 0; node < 8; ++node) {
                int col       = 3 * node + i; // Column: DOF for u_i at node
                B_F(row, col) = dN(J, node);
            }
        }
    }

    return B_F;
}

#endif // LARGESTRAINMECHMODEL_H
