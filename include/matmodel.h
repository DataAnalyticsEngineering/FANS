
#ifndef MATMODEL_H
#define MATMODEL_H

#include "general.h"

constexpr int get_n_str(int h)
{
    switch (h) {
    case 1:
        return 3;
    case 3:
        return 6;
    }
}
template <int howmany>
class Solver;

template <int howmany>
class Matmodel {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html

        const static int n_str = get_n_str(howmany); // length of strain and stress

    int verbosity; //!< output verbosity
    int n_mat;     // Number of Materials

    double *strain; //!< Gradient
    double *stress; //!< Flux

    Matmodel(vector<double> l_e);

    Matrix<double, howmany * 8, howmany * 8> Compute_Reference_ElementStiffness();
    Matrix<double, howmany * 8, 1>          &element_residual(Matrix<double, howmany * 8, 1> &ue, int mat_index, ptrdiff_t element_idx);
    void                                     getStrainStress(double *strain, double *stress, Matrix<double, howmany * 8, 1> &ue, int mat_index, ptrdiff_t element_idx);
    void                                     setGradient(vector<double> _g0);

    virtual void postprocess(Solver<howmany> &solver, Reader &reader, const char *resultsFileName, int load_idx, int time_idx) {}

    virtual void initializeInternalVariables(ptrdiff_t num_elements, int num_gauss_points) {}
    virtual void updateInternalVariables() {}

    vector<double> macroscale_loading;

  protected:
    double l_e_x;
    double l_e_y;
    double l_e_z;
    double v_e;

    Matrix<double, n_str, n_str> kapparef_mat; // Reference conductivity matrix

    Matrix<double, n_str, howmany * 8>     B_el_mean; //!< precomputed mean B matrix over the element
    Matrix<double, n_str, howmany * 8>     B_int[8];  //!< precomputed B matrix at all integration points
    Matrix<double, 8 * n_str, howmany * 8> B;

    Matrix<double, n_str * 8, 1>   eps;
    Matrix<double, n_str * 8, 1>   g0; //!< Macro-scale Gradient
    Matrix<double, n_str * 8, 1>   sigma;
    Matrix<double, howmany * 8, 1> res_e;

    Matrix<double, 3, 8>                       Compute_basic_B(const double x, const double y, const double z) const;
    virtual Matrix<double, n_str, howmany * 8> Compute_B(const double x, const double y, const double z) = 0;
    void                                       Construct_B();

    virtual void get_sigma(int i, int mat_index, ptrdiff_t element_idx) = 0;
};

template <int howmany>
Matmodel<howmany>::Matmodel(vector<double> l_e)
{
    l_e_x = l_e[0];
    l_e_y = l_e[1];
    l_e_z = l_e[2];

    v_e = l_e_x * l_e_y * l_e_z;
}

template <int howmany>
void Matmodel<howmany>::Construct_B()
{
    const double xi_p     = 0.5 + sqrt(3.) / 6.;
    const double xi_m     = 0.5 - sqrt(3.) / 6.;
    const double xi[8][3] = {{xi_m, xi_m, xi_m}, {xi_p, xi_m, xi_m}, {xi_m, xi_p, xi_m}, {xi_p, xi_p, xi_m}, {xi_m, xi_m, xi_p}, {xi_p, xi_m, xi_p}, {xi_m, xi_p, xi_p}, {xi_p, xi_p, xi_p}};

    B_el_mean = Compute_B(0.5, 0.5, 0.5);

    // fetch B at the integration sites
    for (int p = 0; p < 8; p++) {
        B_int[p]                                  = Compute_B(xi[p][0], xi[p][1], xi[p][2]);
        B.block(n_str * p, 0, n_str, howmany * 8) = B_int[p];
    }
}

template <int howmany>
Matrix<double, 3, 8> Matmodel<howmany>::Compute_basic_B(const double x, const double y, const double z) const
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

template <int howmany>
Matrix<double, howmany * 8, 1> &Matmodel<howmany>::element_residual(Matrix<double, howmany * 8, 1> &ue, int mat_index, ptrdiff_t element_idx)
{

    eps.noalias() = B * ue + g0;

    for (int i = 0; i < 8; ++i) {
        get_sigma(n_str * i, mat_index, element_idx);
    }
    res_e.noalias() = B.transpose() * sigma * v_e * 0.125;
    return res_e;
}
template <int howmany>
void Matmodel<howmany>::getStrainStress(double *strain, double *stress, Matrix<double, howmany * 8, 1> &ue, int mat_index, ptrdiff_t element_idx)
{
    eps.noalias() = B * ue + g0;
    sigma.setZero();
    for (int i = 0; i < 8; ++i) {
        get_sigma(n_str * i, mat_index, element_idx);
    }

    Matrix<double, n_str, 1> avg_strain = Matrix<double, n_str, 1>::Zero();
    Matrix<double, n_str, 1> avg_stress = Matrix<double, n_str, 1>::Zero();

    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < n_str; ++j) {
            avg_strain(j) += eps(i * n_str + j) * 0.125;
            avg_stress(j) += sigma(i * n_str + j) * 0.125;
        }
    }

    for (int i = 0; i < n_str; ++i) {
        strain[i] = avg_strain(i);
        stress[i] = avg_stress(i);
    }
}

template <int howmany>
void Matmodel<howmany>::setGradient(vector<double> _g0)
{
    macroscale_loading = _g0;
    for (int i = 0; i < n_str; i++) {
        for (int j = 0; j < 8; ++j) {
            g0(n_str * j + i, 0) = _g0[i];
        }
    }
}
template <int howmany>
Matrix<double, howmany * 8, howmany * 8> Matmodel<howmany>::Compute_Reference_ElementStiffness()
{
    Matrix<double, howmany * 8, howmany * 8> Reference_ElementStiffness = Matrix<double, howmany * 8, howmany * 8>::Zero();
    Matrix<double, howmany * 8, howmany * 8> tmp                        = Matrix<double, howmany * 8, howmany * 8>::Zero();

    for (int p = 0; p < 8; ++p) {
        tmp += B_int[p].transpose() * kapparef_mat * B_int[p] * v_e * 0.1250;
    }
    // before: 8 groups of howmany      after: howmany groups of 8
    for (int i = 0; i < howmany * 8; ++i) {
        for (int j = 0; j < howmany * 8; ++j) {
            Reference_ElementStiffness((i % howmany) * 8 + i / howmany, (j % howmany) * 8 + j / howmany) = tmp(i, j);
        }
    }
    return Reference_ElementStiffness;
}

class ThermalModel : public Matmodel<1> {
  public:
    ThermalModel(vector<double> l_e)
        : Matmodel(l_e)
    {
        Construct_B();
    };

  protected:
    Matrix<double, 3, 8> Compute_B(const double x, const double y, const double z);
};

inline Matrix<double, 3, 8> ThermalModel::Compute_B(const double x, const double y, const double z)
{
    return Matmodel<1>::Compute_basic_B(x, y, z);
}

class MechModel : public Matmodel<3> {
  public:
    MechModel(vector<double> l_e)
        : Matmodel(l_e)
    {
        Construct_B();
    };

  protected:
    Matrix<double, 6, 24> Compute_B(const double x, const double y, const double z);
};

inline Matrix<double, 6, 24> MechModel::Compute_B(const double x, const double y, const double z)
{
    Matrix<double, 6, 24> out       = Matrix<double, 6, 24>::Zero();
    Matrix<double, 3, 8>  B_tmp     = Matmodel<3>::Compute_basic_B(x, y, z);
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

template <int howmany>
class LinearModel {
  public:
    Matrix<double, howmany * 8, howmany * 8> *phase_stiffness;
};

#endif // MATMODEL_H
