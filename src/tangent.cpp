// The homogenized tangent of a converged step.
//
// If every material has a tangent C = d flux / d gradient, the tangents of the
// step, frozen at every Gauss point, are a linear material: flux = C gradient.
// Its homogenized stiffness is the consistent homogenized tangent, and FANS
// finds it as for any linear material: a solve per unit macroscopic gradient,
// whose mean flux is a column.
// Otherwise the materials themselves are perturbed, by a finite difference.

#include "solver.h"
#include "LargeStrainMechModel.h"
#include "material_models/small_strain/J2Plasticity.h"

namespace {

// The B matrices of the problem
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
    const auto     &materials = *solver.matmanager;
    const ptrdiff_t n_slab    = solver.n_y * solver.n_z * howmany; // nodal values of a slab of the grid
    const size_t    per_elem  = C.size() / (solver.local_n0 * solver.n_y * solver.n_z);
    double         *u         = solver.v_u;

    // The native materials element by element, as Solver::get_homogenized_stress goes through them
    MPI_Sendrecv(u, n_slab, MPI_DOUBLE, (solver.world_rank + solver.world_size - 1) % solver.world_size, 0,
                 u + solver.local_n0 * n_slab, n_slab, MPI_DOUBLE, (solver.world_rank + 1) % solver.world_size, 0, solver.communicator, MPI_STATUS_IGNORE);
    Matrix<double, howmany * 8, 1> ue;
    solver.template iterateCubes<0>([&](ptrdiff_t *idx, ptrdiff_t *) {
        const MaterialInfo<howmany, n_str> &info = materials.get_info(solver.ms[idx[0]]);
        if (info.model->wants_batch())
            return;
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < howmany; ++j)
                ue(howmany * i + j) = u[howmany * idx[i] + j];
        info.model->getTangent(&C[idx[0] * per_elem], ue, info.local_mat_id, idx[0]);
    });

    // The batched materials in batches, from the nodal values of the step, which batch_ue still holds
    for (size_t m = 0; m < materials.models.size(); ++m)
        if (materials.models[m]->wants_batch())
            materials.models[m]->evaluate_batch(solver.batch_elems[m], solver.ms, solver.batch_ue.data(), solver.batch_gp_stress.data(), C.data());
}

} // namespace

template <int howmany, int n_str>
MatrixXd Solver<howmany, n_str>::get_homogenized_tangent(double pert_param)
{
    // The step, to come back to
    MaterialManager<howmany, n_str> *materials = matmanager;
    const vector<double>             g0        = matmanager->models[0]->macroscale_loading;
    const vector<double>             u(v_u, v_u + local_n0 * n_y * n_z * howmany);

    // Does every material have a tangent? Linear models with element stiffnesses are
    // left as they are: perturbing them is exact, and faster.
    bool consistent = !matmanager->all_stiffness;
    for (int phase = 0; phase < reader.n_mat; ++phase) {
        const MaterialInfo<howmany, n_str> &info = matmanager->get_info(phase);
        consistent                               = consistent && (info.is_linear || info.model->has_tangent());
    }

    // If so, the frozen tangents take the place of the materials
    std::unique_ptr<MaterialManager<howmany, n_str>> linearised;
    if (consistent) {
        auto *frozen = new FrozenTangent<n_str>(reader, local_n0 * n_y * n_z);
        collect_tangents(*this, frozen->C);
        linearised = std::make_unique<MaterialManager<howmany, n_str>>(frozen, reader.n_mat); // owns `frozen`
        matmanager = linearised.get();
    } else {
        for (auto *model : matmanager->models) // perturbing would change its history
            if (dynamic_cast<J2Plasticity *>(model) != nullptr)
                throw std::runtime_error("Homogenized tangent computation not implemented for J2Plasticity models.");
    }
    Log::logger().info("# Homogenized tangent {}: {} solves", consistent ? "from the materials' tangents" : "by perturbation", n_str);

    const bool linear = matmanager->all_linear;
    VectorXd   stress;
    if (!linear)
        stress = get_homogenized_stress();
    disableMixedBC();

    homogenized_tangent.resize(n_str, n_str);
    for (int j = 0; j < n_str; ++j) {
        vector<double> gradient(n_str, 0.0);
        if (linear) { // a unit gradient, from rest
            gradient[j] = 1.0;
            std::fill_n(v_u, u.size(), 0.0);
        } else { // a small step from the gradient of the step
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
    if (!consistent)
        homogenized_tangent = 0.5 * (homogenized_tangent + homogenized_tangent.transpose()).eval();

    // Back to the step
    matmanager = materials;
    matmanager->set_gradient(g0);
    std::copy(u.begin(), u.end(), v_u);
    return homogenized_tangent;
}

template MatrixXd Solver<1, 3>::get_homogenized_tangent(double);
template MatrixXd Solver<3, 6>::get_homogenized_tangent(double);
template MatrixXd Solver<3, 9>::get_homogenized_tangent(double);
