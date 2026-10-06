// The homogenized tangent of a converged step from the materials' own tangents.
//
// With the tangent C = d flux / d gradient of that step frozen at every Gauss
// point, the linearised microstructure is a linear material, flux = C gradient.
// FANS homogenizes it as any linear material: a solve per unit macroscopic
// gradient, whose mean flux is a column of the homogenized tangent.

#include "solver.h"
#include "LargeStrainMechModel.h"

namespace {

template <int n_str>
using Kinematics = std::conditional_t<n_str == 3, ThermalModel, std::conditional_t<n_str == 6, SmallStrainMechModel, LargeStrainMechModel>>;

// flux = C gradient, with C [element][n_gp][n_str][n_str]
template <int n_str>
class FrozenTangent : public Kinematics<n_str> {
  public:
    FrozenTangent(const Reader &reader, size_t n_elem)
        : Kinematics<n_str>(reader), C(n_elem * this->n_gp * n_str * n_str) {}

    vector<double> C;

    void get_sigma(int i, int, ptrdiff_t element_idx) override
    {
        const double *C_gp = &C[(element_idx * this->n_gp * n_str + i) * n_str];
        this->sigma.template segment<n_str>(i).noalias() =
            Map<const Matrix<double, n_str, n_str, RowMajor>>(C_gp) * this->eps.template segment<n_str>(i);
    }

    Matrix<double, n_str, n_str> get_reference_stiffness() override
    {
        throw std::logic_error("FrozenTangent: the solver already has its reference material");
    }
};

} // namespace

// False if a material has no tangent of its own: then get_homogenized_tangent perturbs
template <int howmany, int n_str>
bool Solver<howmany, n_str>::consistent_homogenized_tangent()
{
    if (!matmanager->any_batched)
        return false;

    // The materials' tangents at the converged step, whose nodal values batch_ue still holds
    const size_t                    n_elem = local_n0 * n_y * n_z;
    auto                           *frozen = new FrozenTangent<n_str>(reader, n_elem);
    MaterialManager<howmany, n_str> linearised(frozen, reader.n_mat); // owns `frozen`
    for (size_t m = 0; m < matmanager->models.size(); ++m)
        if (!matmanager->models[m]->evaluate_tangent(batch_elems[m], ms, batch_ue.data(), frozen->C.data()))
            return false;

    Log::logger().info("# Homogenized tangent from the materials' tangents: {} linear solves", n_str);
    MaterialManager<howmany, n_str> *materials = matmanager;
    const vector<double>             u(v_u, v_u + n_elem * howmany);
    matmanager = &linearised;
    disableMixedBC();

    homogenized_tangent.resize(n_str, n_str);
    for (int j = 0; j < n_str; ++j) {
        vector<double> unit(n_str, 0.0);
        unit[j] = 1.0;
        matmanager->set_gradient(unit);
        std::fill_n(v_u, u.size(), 0.0);
        solve();
        homogenized_tangent.col(j) = get_homogenized_stress();
    }

    // Back to the materials and the converged step
    matmanager = materials;
    std::copy(u.begin(), u.end(), v_u);
    return true;
}

template bool Solver<1, 3>::consistent_homogenized_tangent();
template bool Solver<3, 6>::consistent_homogenized_tangent();
template bool Solver<3, 9>::consistent_homogenized_tangent();
