#ifndef PSEUDOPLASTICLINEARHARDENING_H
#define PSEUDOPLASTICLINEARHARDENING_H

#include "matmodel.h"
#include "solver.h"

class PseudoPlasticLinearHardening : public MechModel {
  public:
    PseudoPlasticLinearHardening(vector<double> l_e, map<string, vector<double>> materialProperties)
        : MechModel(l_e)
    {
        bulk_modulus        = materialProperties["bulk_modulus"];
        shear_modulus       = materialProperties["shear_modulus"];
        yield_stress        = materialProperties["yield_stress"];
        hardening_parameter = materialProperties["hardening_parameter"];

        a     = 2. / 3;
        b     = sqrt(a);
        n_mat = bulk_modulus.size();

        eps_crit.resize(n_mat);
        E_s.resize(n_mat);
        Matrix<double, 6, 6> *Ce      = new Matrix<double, 6, 6>[n_mat];
        Matrix<double, 6, 6>  topLeft = Matrix<double, 6, 6>::Zero();
        topLeft.topLeftCorner(3, 3).setConstant(1);

        for (int i = 0; i < n_mat; ++i) {
            eps_crit[i] = sqrt(2. / 3.) * yield_stress[i] / (2. * shear_modulus[i]);
            E_s[i]      = (3. * shear_modulus[i]) / (3. * shear_modulus[i] + hardening_parameter[i]);

            Ce[i] = 3 * bulk_modulus[i] * topLeft +
                    2 * shear_modulus[i] * (-1. / 3 * topLeft + Matrix<double, 6, 6>::Identity());
        }
        kapparef_mat = 0.5 * (Ce[0] + Ce[1]); // Note: only works for two materials
    }

    void initializeInternalVariables(ptrdiff_t num_elements, int num_gauss_points) override
    {
        plastic_flag.resize(num_elements, VectorXi::Zero(num_gauss_points));
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        treps             = eps.block<3, 1>(i, 0).sum();
        dev_eps.head<3>() = eps.block<3, 1>(i, 0) - (1.0 / 3.0) * treps * Matrix<double, 3, 1>::Ones();
        dev_eps.tail<3>() = eps.block<3, 1>(i + 3, 0);

        norm_dev_eps = dev_eps.norm();
        buf1         = bulk_modulus[mat_index] * treps;

        if (norm_dev_eps <= eps_crit[mat_index]) {
            buf2                                 = 2.0 * shear_modulus[mat_index];
            plastic_flag[element_idx](i / n_str) = mat_index;
        } else {
            buf2 = (b * yield_stress[mat_index] + a * E_s[mat_index] * hardening_parameter[mat_index] *
                                                      (norm_dev_eps - eps_crit[mat_index])) /
                   norm_dev_eps;
            plastic_flag[element_idx](i / n_str) = this->n_mat + mat_index;
        }
        sigma.block<3, 1>(i, 0).setConstant(buf1);
        sigma.block<3, 1>(i, 0) += buf2 * dev_eps.head<3>();
        sigma.block<3, 1>(i + 3, 0) = buf2 * dev_eps.tail<3>();
    }

    void postprocess(Solver<3> &solver, Reader &reader, const char *resultsFileName, int load_idx, int time_idx) override;

  private:
    vector<double>   bulk_modulus;
    vector<double>   shear_modulus;
    vector<double>   yield_stress;
    vector<double>   hardening_parameter;
    vector<double>   eps_crit;
    vector<double>   E_s;
    vector<VectorXi> plastic_flag;

    Matrix<double, 6, 1> dev_eps;
    double               treps, norm_dev_eps, buf1, buf2;
    double               a;
    double               b;
};

void PseudoPlasticLinearHardening::postprocess(Solver<3> &solver, Reader &reader, const char *resultsFileName, int load_idx, int time_idx)
{

    VectorXf element_plastic_flag = VectorXf::Zero(solver.local_n0 * solver.n_y * solver.n_z);
    for (ptrdiff_t elem_idx = 0; elem_idx < solver.local_n0 * solver.n_y * solver.n_z; ++elem_idx) {
        element_plastic_flag(elem_idx) = plastic_flag[elem_idx].cast<float>().mean();
    }

    if (find(reader.resultsToWrite.begin(), reader.resultsToWrite.end(), "plastic_flag") != reader.resultsToWrite.end()) {
        for (int i = 0; i < solver.world_size; ++i) {
            if (i == solver.world_rank) {
                char name[5096];
                sprintf(name, "%s/load%i/time_step%i/plastic_flag", reader.ms_datasetname, load_idx, time_idx);
                reader.WriteSlab<float>(element_plastic_flag.data(), 1, resultsFileName, name);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}

#endif // PSEUDOPLASTICLINEARHARDENING_H
