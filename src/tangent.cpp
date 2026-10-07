// The homogenized tangent of a converged step.
//
// If every material has a tangent C = d flux / d gradient, the tangents of the
// step, frozen at every Gauss point, are a linear material: flux = C gradient.
// Its homogenized stiffness is the consistent homogenized tangent, and FANS
// finds it as for any linear material: a solve per unit macroscopic gradient,
// whose mean flux is a column.
// Otherwise the materials themselves are perturbed, by a finite difference. That
// is not for materials with internal variables: every perturbed solve advances them.

#include "solver.h"
#include "LargeStrainMechModel.h"

namespace {

// The model base of the problem: thermal, small or large strain
template <int n_str>
using Kinematics = std::conditional_t<n_str == 3, ThermalModel, std::conditional_t<n_str == 6, SmallStrainMechModel, LargeStrainMechModel>>;

// The linear material flux = C gradient, C being [element][n_gp][n_str][n_str]
template <int n_str>
class FrozenTangent : public Kinematics<n_str> {
  public:
    FrozenTangent(const Reader &reader, size_t n_elem)
        : Kinematics<n_str>(reader), C(n_elem * this->n_gp * n_str * n_str) {}

    vector<double> C;

    void get_sigma(int i, int, ptrdiff_t element_idx) override
    {
        const size_t                                      gp = element_idx * this->n_gp + i / n_str;
        Map<const Matrix<double, n_str, n_str, RowMajor>> C_gp(&C[gp * n_str * n_str]);
        this->sigma.template segment<n_str>(i).noalias() = C_gp * this->eps.template segment<n_str>(i);
    }

    // Not needed: the solver has its reference material already
    Matrix<double, n_str, n_str> get_reference_stiffness() override
    {
        return Matrix<double, n_str, n_str>::Zero();
    }
};

// The tangent of every material at the converged step, into C
template <int howmany, int n_str>
void collect_tangents(Solver<howmany, n_str> &solver, vector<double> &C)
{
    const auto  &materials = *solver.matmanager;
    const size_t per_elem  = C.size() / (solver.local_n0 * solver.n_y * solver.n_z);

    // The native materials element by element, as Solver::get_homogenized_stress goes through them
    Matrix<double, howmany * 8, 1> ue;
    solver.template iterateCubes<0>([&](ptrdiff_t *idx, ptrdiff_t *) {
        const MaterialInfo<howmany, n_str> &info = materials.get_info(solver.ms[idx[0]]);
        if (info.model->wants_batch())
            return;
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < howmany; ++j)
                ue(howmany * i + j) = solver.v_u[howmany * idx[i] + j];
        info.model->getTangent(&C[idx[0] * per_elem], ue, info.local_mat_id, idx[0]);
    });

    // The batched materials in batches
    for (size_t m = 0; m < materials.models.size(); ++m)
        if (materials.models[m]->wants_batch())
            materials.models[m]->evaluate_batch(solver.batch_elems[m], solver.ms, solver.batch_ue.data(), solver.batch_gp_stress.data(), C.data());
}

} // namespace

template <int howmany, int n_str>
MatrixXd Solver<howmany, n_str>::get_homogenized_tangent(double pert_param)
{
    // The step: its stress (evaluating it also readies the nodal values the tangents are
    // collected from), and what to come back to
    const VectorXd                   stress    = get_homogenized_stress();
    MaterialManager<howmany, n_str> *materials = matmanager;
    const vector<double>             g0        = matmanager->models[0]->macroscale_loading;
    const vector<double>             u(v_u, v_u + local_n0 * n_y * n_z * howmany);

    // Does every material have a tangent? Linear models with element stiffnesses are
    // left as they are: perturbing them is exact, and faster.
    bool consistent = !matmanager->all_stiffness;
    for (int phase = 0; phase < matmanager->get_num_phases(); ++phase)
        if (!matmanager->get_info(phase).is_linear && !matmanager->get_info(phase).model->has_tangent())
            consistent = false;

    // If so, the frozen tangents take the place of the materials
    std::unique_ptr<MaterialManager<howmany, n_str>> linearised;
    if (consistent) {
        auto *frozen = new FrozenTangent<n_str>(reader, local_n0 * n_y * n_z);
        collect_tangents(*this, frozen->C);
        linearised = std::make_unique<MaterialManager<howmany, n_str>>(frozen, matmanager->get_num_phases()); // owns `frozen`
        matmanager = linearised.get();
        Log::logger().info("# Homogenized tangent from the materials' tangents");
    } else if (!matmanager->all_linear) {
        Log::logger().warn("# Homogenized tangent by perturbation: For materials with internal variables, this is WRONG!");
    }

    const bool linear = matmanager->all_linear;
    disableMixedBC();

    homogenized_tangent.resize(n_str, n_str);
    for (int j = 0; j < n_str; ++j) {
        vector<double> gradient(n_str, 0.0);
        if (linear) {
            gradient[j] = 1.0;
            std::fill_n(v_u, u.size(), 0.0);
        } else {
            gradient = g0;
            gradient[j] += pert_param;
        }
        matmanager->set_gradient(gradient);
        solve();

        if (linear)
            homogenized_tangent.col(j) = get_homogenized_stress();
        else
            homogenized_tangent.col(j) = (get_homogenized_stress() - stress) / pert_param;
    }
    homogenized_tangent = 0.5 * (homogenized_tangent + homogenized_tangent.transpose()).eval();

    // Back to the step
    matmanager = materials;
    std::copy(u.begin(), u.end(), v_u);
    return homogenized_tangent;
}

template MatrixXd Solver<1, 3>::get_homogenized_tangent(double);
template MatrixXd Solver<3, 6>::get_homogenized_tangent(double);
template MatrixXd Solver<3, 9>::get_homogenized_tangent(double);
