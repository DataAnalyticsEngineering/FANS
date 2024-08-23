#ifndef VONMISESPLASTICLINEARISOTROPICHARDENING_H
#define VONMISESPLASTICLINEARISOTROPICHARDENING_H

#include "matmodel.h"
#include "solver.h"

class VonMisesPlasticLinearIsotropicHardening : public MechModel {
public:
    VonMisesPlasticLinearIsotropicHardening(vector<double> l_e, map<string, vector<double>> materialProperties)
        : MechModel(l_e) {
        bulk_modulus = materialProperties["bulk_modulus"];
        shear_modulus = materialProperties["shear_modulus"];
        yield_stress = materialProperties["yield_stress"];
        hardening_parameter = materialProperties["hardening_parameter"];

        n_mat = bulk_modulus.size();

        two_G_plus_H.resize(n_mat);
        sqrt_two_over_three = sqrt(2.0 / 3.0);

        Matrix<double, 6, 6>* Ce = new Matrix<double, 6, 6>[n_mat];
        Matrix<double, 6, 6> topLeft = Matrix<double, 6, 6>::Zero();
        topLeft.topLeftCorner(3, 3).setConstant(1);

        for (int i = 0; i < n_mat; ++i) {
            Ce[i] = 3 * bulk_modulus[i] * topLeft +
                    2 * shear_modulus[i] * (-1.0 / 3.0 * topLeft + Matrix<double, 6, 6>::Identity());
            two_G_plus_H[i] = 2 * shear_modulus[i] + (2.0 / 3.0) * hardening_parameter[i];
        }
        kapparef_mat = 0.5 * (Ce[0] + Ce[1]); // Note: only works for two materials
    }

    void initializeInternalVariables(ptrdiff_t num_elements, int num_gauss_points) override {
        plastic_strain.resize(num_elements, Matrix<double, 6, Dynamic>::Zero(6, num_gauss_points));
        hardening_variable.resize(num_elements, VectorXd::Zero(num_gauss_points));

        plastic_strain_t.resize(num_elements, Matrix<double, 6, Dynamic>::Zero(6, num_gauss_points));
        hardening_variable_t.resize(num_elements, VectorXd::Zero(num_gauss_points));
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override {

        // Elastic Predictor
        eps_elastic = eps.template block<6, 1>(i, 0) - plastic_strain_t[element_idx].col(i / n_str);
        treps = eps_elastic.head<3>().sum();

        // Compute trial stress
        stress_trial.head<3>().setConstant(bulk_modulus[mat_index] * treps);
        stress_trial.head<3>() += 2 * shear_modulus[mat_index] * eps_elastic.head<3>();
        stress_trial.tail<3>() = 2 * shear_modulus[mat_index] * eps_elastic.tail<3>();

        // Deviatoric stress and von Mises equivalent stress
        sigma_hydrostatic = stress_trial.head<3>().mean();
        dev_stress = stress_trial;
        dev_stress.head<3>().array() -= sigma_hydrostatic;
        norm_dev_stress = dev_stress.head<3>().norm();

        // Yield function
        yield_function = norm_dev_stress - sqrt_two_over_three * (yield_stress[mat_index] + hardening_parameter[mat_index] * hardening_variable_t[element_idx](i / n_str));

        if (yield_function > 0) {
            // Plastic correction
            delta_lambda = yield_function / two_G_plus_H[mat_index];

            // Update plastic strain and hardening variable
            plastic_strain_increment = delta_lambda * (dev_stress / norm_dev_stress);
            plastic_strain[element_idx].col(i / n_str) = plastic_strain_t[element_idx].col(i / n_str) + plastic_strain_increment;

            // Update hardening variable
            hardening_variable[element_idx](i / n_str) = hardening_variable_t[element_idx](i / n_str) + sqrt_two_over_three * delta_lambda;

            // Update stress
            stress_trial -= 2 * shear_modulus[mat_index] * plastic_strain_increment;
        }

        // Final stress
        sigma.template block<6, 1>(i, 0) = stress_trial;
    }

    void postprocess(Solver<3>& solver, Reader& reader, const char* resultsFileName, int load_idx, int time_idx) override;

private:
    vector<double> bulk_modulus;
    vector<double> shear_modulus;
    vector<double> yield_stress;
    vector<double> hardening_parameter;
    vector<double> two_G_plus_H;
    double sqrt_two_over_three;

    // Internal variables
    vector<Matrix<double, 6, Dynamic>> plastic_strain; // Stores plastic strain at each Gauss point of each element
    vector<VectorXd> hardening_variable; // Stores hardening variable at each Gauss point of each element

    vector<Matrix<double, 6, Dynamic>> plastic_strain_t; // Stores plastic strain at each Gauss point of each element
    vector<VectorXd> hardening_variable_t; // Stores hardening variable at each Gauss point of each element

    Matrix<double, 6, 1> eps_elastic;
    Matrix<double, 6, 1> stress_trial;
    Matrix<double, 6, 1> dev_stress;
    Matrix<double, 6, 1> plastic_strain_increment;
    double norm_dev_stress;
    double yield_function;
    double treps;
    double sigma_hydrostatic;
    double delta_lambda;
};

void VonMisesPlasticLinearIsotropicHardening::postprocess(Solver<3>& solver, Reader& reader, const char* resultsFileName, int load_idx, int time_idx) {
    plastic_strain_t = plastic_strain;
    hardening_variable_t = hardening_variable;

    int n_str = 6; // The plastic strain and stress vectors have 6 components each
    VectorXd mean_plastic_strain = VectorXd::Zero(solver.local_n0 * solver.n_y * solver.n_z * n_str);
    VectorXd mean_hardening_variable = VectorXd::Zero(solver.local_n0 * solver.n_y * solver.n_z);

    // Compute the mean values for each element
    for (ptrdiff_t elem_idx = 0; elem_idx < solver.local_n0 * solver.n_y * solver.n_z; ++elem_idx) {
        mean_plastic_strain.segment(n_str * elem_idx, n_str) = plastic_strain[elem_idx].rowwise().mean();
        mean_hardening_variable(elem_idx) = hardening_variable[elem_idx].mean();
    }

    if (find(reader.resultsToWrite.begin(), reader.resultsToWrite.end(), "plastic_strain") != reader.resultsToWrite.end()) {
        for (int i = 0; i < solver.world_size; ++i) {
            if (i == solver.world_rank) {
                char name[5096];
                sprintf(name, "%s/load%i/time_step%i/plastic_strain", reader.ms_datasetname, load_idx, time_idx);
                reader.WriteSlab<double>(mean_plastic_strain.data(), n_str, resultsFileName, name);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    if (find(reader.resultsToWrite.begin(), reader.resultsToWrite.end(), "hardening_variable") != reader.resultsToWrite.end()) {
        for (int i = 0; i < solver.world_size; ++i) {
            if (i == solver.world_rank) {
                char name[5096];
                sprintf(name, "%s/load%i/time_step%i/hardening_variable", reader.ms_datasetname, load_idx, time_idx);
                reader.WriteSlab<double>(mean_hardening_variable.data(), 1, resultsFileName, name);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

}

#endif // VONMISESPLASTICLINEARISOTROPICHARDENING_H
