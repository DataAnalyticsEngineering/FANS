#ifndef LARGESTRAINMECHMODEL_H
#define LARGESTRAINMECHMODEL_H

#include "matmodel.h"

/**
 * @brief Base class for large strain mechanical material models
 * This class provides the kinematic framework for finite deformation mechanics.
 * Derived classes only need to implement the constitutive response (S from E).
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
     * @brief Compute 2nd Piola-Kirchhoff stress and material tangent C = dS/dE from deformation gradient
     *
     * This is the core constitutive function that derived classes must implement.
     * It should compute S based on the material's constitutive law.
     * @param F Deformation gradient (3×3)
     * @param mat_index Material index
     * @param element_idx Element index
     * @return 2nd Piola-Kirchhoff stress S (symmetric 3×3)
     * @return the 6×6 material tangent C = dS/dE in Mandel notation: Ordering: (11, 22, 33, sqrt(2)·12, sqrt(2)·13, sqrt(2)·23)
     */
    virtual Matrix3d             compute_S(const Matrix3d &F, int mat_index, ptrdiff_t element_idx) = 0;
    virtual Matrix<double, 6, 6> compute_material_tangent(const Matrix3d &F, int mat_index)         = 0;

    /**
     * @brief Convert material tangent C (dS/dE) to spatial tangent A (dP/dF)
     *
     * Computes the spatial tangent dP/dF from the material tangent dS/dE.
     *
     * The full expression is:
     *   dP_iJ/dF_kL = δ_ik S_LJ + F_iM * dS_MJ/dE_PQ * dE_PQ/dF_kL
     * where:
     *   dE_PQ/dF_kL = 0.5 * (F_kP δ_QL + F_kQ δ_PL)
     *
     * @param F Deformation gradient (3×3)
     * @param S 2nd Piola-Kirchhoff stress (3×3 symmetric)
     * @param C_mandel Material tangent dS/dE in Mandel notation (6×6)
     * @return Spatial tangent A = dP/dF in full notation (9×9)
     */
    inline Matrix<double, 9, 9> compute_spatial_tangent(const Matrix3d             &F,
                                                        const Matrix3d             &S,
                                                        const Matrix<double, 6, 6> &C_mandel) const
    {
        Matrix<double, 9, 9> A = Matrix<double, 9, 9>::Zero();

        // Mandel to tensor index mapping
        int mandel_to_ij[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}};

        // Compute A_iJkL = dP_iJ/dF_kL
        for (int i = 0; i < 3; ++i) {
            for (int J = 0; J < 3; ++J) {
                int row = 3 * i + J; // Index for P_iJ

                for (int k = 0; k < 3; ++k) {
                    for (int L = 0; L < 3; ++L) {
                        int col = 3 * k + L; // Index for F_kL

                        // First term: δ_ik S_LJ
                        if (i == k) {
                            A(row, col) += S(L, J);
                        }

                        // Second term: F_iM * C_MJPQ * dE_PQ/dF_kL
                        for (int M = 0; M < 3; ++M) {
                            // Find MJ in Mandel notation
                            int MJ_mandel = -1;
                            for (int idx = 0; idx < 6; ++idx) {
                                if ((mandel_to_ij[idx][0] == M && mandel_to_ij[idx][1] == J) ||
                                    (mandel_to_ij[idx][0] == J && mandel_to_ij[idx][1] == M)) {
                                    MJ_mandel = idx;
                                    break;
                                }
                            }
                            if (MJ_mandel < 0)
                                continue;

                            for (int P = 0; P < 3; ++P) {
                                for (int Q = P; Q < 3; ++Q) { // Q >= P due to symmetry
                                    // Find PQ in Mandel notation
                                    int PQ_mandel = -1;
                                    for (int idx = 0; idx < 6; ++idx) {
                                        if ((mandel_to_ij[idx][0] == P && mandel_to_ij[idx][1] == Q) ||
                                            (mandel_to_ij[idx][0] == Q && mandel_to_ij[idx][1] == P)) {
                                            PQ_mandel = idx;
                                            break;
                                        }
                                    }
                                    if (PQ_mandel < 0)
                                        continue;

                                    // Get C_MJPQ and account for Mandel factors
                                    double C_val = C_mandel(MJ_mandel, PQ_mandel);
                                    if (MJ_mandel >= 3)
                                        C_val /= sqrt(2.0);
                                    if (PQ_mandel >= 3)
                                        C_val /= sqrt(2.0);

                                    // Compute dE_PQ/dF_kL = 0.5 * (F_kP δ_QL + F_kQ δ_PL)
                                    double dE_dF = 0.0;
                                    if (Q == L)
                                        dE_dF += 0.5 * F(k, P);
                                    if (P == L)
                                        dE_dF += 0.5 * F(k, Q);

                                    A(row, col) += F(i, M) * C_val * dE_dF;
                                }
                            }
                        }
                    }
                }
            }
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
     * @brief Compute stress at Gauss point (implements base class interface)
     * This function:
     * 1. Extracts F from eps
     * 2. Calls compute_S(F) - material model implements this
     * 3. Computes P = F S
     * 4. Stores P in sigma
     */
    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override final
    {
        Matrix3d F = extract_F(i);                         // Extract deformation gradient
        Matrix3d S = compute_S(F, mat_index, element_idx); // Material model computes 2nd PK stress from F
        Matrix3d P = push_forward(F, S);                   // Push forward to 1st PK stress
        store_P(i, P);                                     // Store in sigma array
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
