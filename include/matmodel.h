
#ifndef MATMODEL_H
#define MATMODEL_H

#include "general.h"

template <int howmany, int n_str>
class Solver;

template <int howmany, int n_str>
class Matmodel {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html

        static constexpr int num_str = n_str; // length of strain and stress

    int    verbosity; //!< output verbosity
    int    n_mat;     //!< Number of Materials
    int    n_gp;      //!< Number of Gauss points (computed from FE_type)
    string FE_type;   //!< Finite element type: "HEX8", "HEX8R", "BBAR"

    double *strain; //!< Gradient
    double *stress; //!< Flux

    Matmodel(const Reader &reader);

    Matrix<double, howmany * 8, howmany * 8> Compute_Reference_ElementStiffness(const Matrix<double, n_str, n_str> &kapparef_mat);
    Matrix<double, howmany * 8, 1>          &element_residual(Matrix<double, howmany * 8, 1> &ue, int mat_index, ptrdiff_t element_idx);
    void                                     getStrainStress(double *strain, double *stress, Matrix<double, howmany * 8, 1> &ue, int mat_index, ptrdiff_t element_idx);
    void                                     setGradient(vector<double> _g0);

    virtual void postprocess(Solver<howmany, n_str> &solver, Reader &reader, int load_idx, int time_idx) {};

    virtual void initializeInternalVariables(ptrdiff_t num_elements, int num_gauss_points) {}
    virtual void updateInternalVariables() {}

    vector<double>                       macroscale_loading;
    virtual Matrix<double, n_str, n_str> get_reference_stiffness() const = 0;

    virtual ~Matmodel() = default;

  protected:
    double l_e_x;
    double l_e_y;
    double l_e_z;
    double v_e;

    vector<Matrix<double, n_str, howmany * 8>> B_int; //!< precomputed B matrix at all integration points (size = n_gp)
    MatrixXd                                   B;     //!< Dynamic B matrix (n_str * n_gp x howmany * 8)

    VectorXd                       eps;   //!< Dynamic strain vector (size = n_str * n_gp)
    VectorXd                       g0;    //!< Macro-scale loading (size = n_str * n_gp)
    VectorXd                       sigma; //!< Dynamic stress vector (size = n_str * n_gp)
    Matrix<double, howmany * 8, 1> res_e;

    Matrix<double, 3, 8>                       Compute_basic_B(const double x, const double y, const double z) const;
    virtual Matrix<double, n_str, howmany * 8> Compute_B(const double x, const double y, const double z) = 0;
    void                                       Construct_B();

    virtual void get_sigma(int i, int mat_index, ptrdiff_t element_idx) = 0;
};

template <int howmany, int n_str>
Matmodel<howmany, n_str>::Matmodel(const Reader &reader)
    : FE_type(reader.FE_type)
{
    // Set number of Gauss points based on FE type
    if (FE_type == "HEX8" || FE_type == "BBAR") {
        n_gp = 8; // Full integration
    } else if (FE_type == "HEX8R") {
        n_gp = 1; // Reduced integration
    } else {
        throw std::runtime_error("Unknown FE_type: '" + FE_type + "'. Supported types: HEX8, HEX8R, BBAR");
    }

    l_e_x = reader.l_e[0];
    l_e_y = reader.l_e[1];
    l_e_z = reader.l_e[2];

    v_e = l_e_x * l_e_y * l_e_z;

    // Resize dynamic matrices based on number of Gauss points
    B_int.resize(n_gp);
    B.resize(n_str * n_gp, howmany * 8);
    eps.resize(n_str * n_gp);
    g0.resize(n_str * n_gp);
    sigma.resize(n_str * n_gp);
}

template <int howmany, int n_str>
void Matmodel<howmany, n_str>::Construct_B()
{
    if (FE_type == "HEX8" || FE_type == "BBAR") {
        // Full integration with 8 Gauss points
        const double xi_p     = 0.5 + sqrt(3.) / 6.;
        const double xi_m     = 0.5 - sqrt(3.) / 6.;
        const double xi[8][3] = {{xi_m, xi_m, xi_m}, {xi_p, xi_m, xi_m}, {xi_m, xi_p, xi_m}, {xi_p, xi_p, xi_m}, {xi_m, xi_m, xi_p}, {xi_p, xi_m, xi_p}, {xi_m, xi_p, xi_p}, {xi_p, xi_p, xi_p}};

        if (FE_type == "BBAR") {
            // BBAR: Selective reduced integration to prevent volumetric locking
            auto B_vol = Compute_B(0.5, 0.5, 0.5);

            for (int p = 0; p < 8; p++) {
                auto B_full = Compute_B(xi[p][0], xi[p][1], xi[p][2]);

                if constexpr (n_str > 3) {
                    for (int col = 0; col < howmany * 8; ++col) {
                        double vol_full = (B_full(0, col) + B_full(1, col) + B_full(2, col)) / 3.0;
                        double vol_bar  = (B_vol(0, col) + B_vol(1, col) + B_vol(2, col)) / 3.0;

                        // B-bar = B_dev + B_vol_bar
                        B_int[p](0, col) = B_full(0, col) - vol_full + vol_bar;
                        B_int[p](1, col) = B_full(1, col) - vol_full + vol_bar;
                        B_int[p](2, col) = B_full(2, col) - vol_full + vol_bar;

                        // Shear strains (deviatoric only)
                        for (int row = 3; row < n_str; ++row) {
                            B_int[p](row, col) = B_full(row, col);
                        }
                    }
                } else {
                    B_int[p] = B_full;
                }

                B.block(n_str * p, 0, n_str, howmany * 8) = B_int[p];
            }
        } else {
            // HEX8: Standard full integration (8 Gauss points)
            for (int p = 0; p < 8; p++) {
                B_int[p]                                  = Compute_B(xi[p][0], xi[p][1], xi[p][2]);
                B.block(n_str * p, 0, n_str, howmany * 8) = B_int[p];
            }
        }
    } else if (FE_type == "HEX8R") {
        // HEX8R: Reduced integration (1 Gauss point)
        B_int[0]                          = Compute_B(0.5, 0.5, 0.5);
        B.block(0, 0, n_str, howmany * 8) = B_int[0];
    } else {
        throw std::runtime_error("Unsupported FE_type in Construct_B: " + FE_type);
    }
}

template <int howmany, int n_str>
Matrix<double, 3, 8> Matmodel<howmany, n_str>::Compute_basic_B(const double x, const double y, const double z) const
{
    Matrix<double, 3, 8> out;
    out(0, 0) = -(1. - y) * (1. - z) / l_e_x;
    out(0, 1) = (1. - y) * (1. - z) / l_e_x;
    out(0, 2) = -y * (1. - z) / l_e_x;
    out(0, 3) = y * (1. - z) / l_e_x;
    out(0, 4) = -(1. - y) * z / l_e_x;
    out(0, 5) = (1. - y) * z / l_e_x;
    out(0, 6) = -y * z / l_e_x;
    out(0, 7) = y * z / l_e_x;

    out(1, 0) = -(1. - x) * (1. - z) / l_e_y;
    out(1, 1) = -x * (1. - z) / l_e_y;
    out(1, 2) = (1. - x) * (1. - z) / l_e_y;
    out(1, 3) = x * (1. - z) / l_e_y;
    out(1, 4) = -(1. - x) * z / l_e_y;
    out(1, 5) = -x * z / l_e_y;
    out(1, 6) = (1. - x) * z / l_e_y;
    out(1, 7) = x * z / l_e_y;

    out(2, 0) = -(1. - x) * (1. - y) / l_e_z;
    out(2, 1) = -x * (1. - y) / l_e_z;
    out(2, 2) = -(1. - x) * y / l_e_z;
    out(2, 3) = -x * y / l_e_z;
    out(2, 4) = (1. - x) * (1. - y) / l_e_z;
    out(2, 5) = x * (1. - y) / l_e_z;
    out(2, 6) = (1. - x) * y / l_e_z;
    out(2, 7) = x * y / l_e_z;
    return out;
}

template <int howmany, int n_str>
Matrix<double, howmany * 8, 1> &Matmodel<howmany, n_str>::element_residual(Matrix<double, howmany * 8, 1> &ue, int mat_index, ptrdiff_t element_idx)
{

    eps.noalias() = B * ue + g0;

    for (int i = 0; i < n_gp; ++i) {
        get_sigma(n_str * i, mat_index, element_idx);
    }
    res_e.noalias() = B.transpose() * sigma * v_e / n_gp;
    return res_e;
}
template <int howmany, int n_str>
void Matmodel<howmany, n_str>::getStrainStress(double *strain, double *stress, Matrix<double, howmany * 8, 1> &ue, int mat_index, ptrdiff_t element_idx)
{
    eps.noalias() = B * ue + g0;
    sigma.setZero();
    for (int i = 0; i < n_gp; ++i) {
        get_sigma(n_str * i, mat_index, element_idx);
    }

    Matrix<double, n_str, 1> avg_strain = Matrix<double, n_str, 1>::Zero();
    Matrix<double, n_str, 1> avg_stress = Matrix<double, n_str, 1>::Zero();

    for (int i = 0; i < n_gp; ++i) {
        for (int j = 0; j < n_str; ++j) {
            avg_strain(j) += eps(i * n_str + j) / n_gp;
            avg_stress(j) += sigma(i * n_str + j) / n_gp;
        }
    }

    for (int i = 0; i < n_str; ++i) {
        strain[i] = avg_strain(i);
        stress[i] = avg_stress(i);
    }
}

template <int howmany, int n_str>
void Matmodel<howmany, n_str>::setGradient(vector<double> _g0)
{
    macroscale_loading = _g0;
    for (int i = 0; i < n_str; i++) {
        for (int j = 0; j < n_gp; ++j) {
            g0(n_str * j + i, 0) = _g0[i];
        }
    }
}
template <int howmany, int n_str>
Matrix<double, howmany * 8, howmany * 8> Matmodel<howmany, n_str>::Compute_Reference_ElementStiffness(const Matrix<double, n_str, n_str> &kapparef_mat)
{
    Matrix<double, howmany * 8, howmany * 8> Reference_ElementStiffness = Matrix<double, howmany * 8, howmany * 8>::Zero();
    Matrix<double, howmany * 8, howmany * 8> tmp                        = Matrix<double, howmany * 8, howmany * 8>::Zero();

    for (int p = 0; p < n_gp; ++p) {
        tmp += B_int[p].transpose() * kapparef_mat * B_int[p] * v_e / n_gp;
    }
    // before: 8 groups of howmany      after: howmany groups of 8
    for (int i = 0; i < howmany * 8; ++i) {
        for (int j = 0; j < howmany * 8; ++j) {
            Reference_ElementStiffness((i % howmany) * 8 + i / howmany, (j % howmany) * 8 + j / howmany) = tmp(i, j);
        }
    }
    return Reference_ElementStiffness;
}

class ThermalModel : public Matmodel<1, 3> {
  public:
    ThermalModel(const Reader &reader)
        : Matmodel(reader)
    {
        Construct_B();
    };

  protected:
    Matrix<double, 3, 8> Compute_B(const double x, const double y, const double z);
};

inline Matrix<double, 3, 8> ThermalModel::Compute_B(const double x, const double y, const double z)
{
    return Matmodel<1, 3>::Compute_basic_B(x, y, z);
}

class SmallStrainMechModel : public Matmodel<3, 6> {
  public:
    SmallStrainMechModel(const Reader &reader)
        : Matmodel(reader)
    {
        Construct_B();
    };

  protected:
    Matrix<double, 6, 24> Compute_B(const double x, const double y, const double z);
};

inline Matrix<double, 6, 24> SmallStrainMechModel::Compute_B(const double x, const double y, const double z)
{
    Matrix<double, 6, 24> out       = Matrix<double, 6, 24>::Zero();
    Matrix<double, 3, 8>  B_tmp     = Matmodel<3, 6>::Compute_basic_B(x, y, z);
    const double          sqrt_half = 7.071067811865476e-01;

    for (int q = 0; q < 8; q++) {
        out(0, 3 * q + 0) = B_tmp(0, q);
        out(1, 3 * q + 1) = B_tmp(1, q);
        out(2, 3 * q + 2) = B_tmp(2, q);

        out(3, 3 * q + 0) = sqrt_half * B_tmp(1, q);
        out(4, 3 * q + 0) = sqrt_half * B_tmp(2, q);
        out(5, 3 * q + 1) = sqrt_half * B_tmp(2, q);

        out(3, 3 * q + 1) = sqrt_half * B_tmp(0, q);
        out(4, 3 * q + 2) = sqrt_half * B_tmp(0, q);
        out(5, 3 * q + 2) = sqrt_half * B_tmp(1, q);
    }
    return out;
}

template <int howmany, int n_str>
class LinearModel {
  public:
    Matrix<double, howmany * 8, howmany * 8> *phase_stiffness;
};

#endif // MATMODEL_H
