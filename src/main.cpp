#include "general.h"
#include "matmodel.h"
#include "setup.h"
#include "solver.h"

// Version
#include "version.h"

template <int howmany, int n_str>
void runSolver(Reader &reader, char input_fn[])
{
    /*-------------------------------------------------------------------------------------------------------------------------------*/
    // Compute elastic homogenized tangent
    Reader reader_linelastic;
    reader_linelastic.ReadInputFile(input_fn);
    reader_linelastic.ReadMS(howmany);
    // Override material models to be linear elastic
    for (auto &mg : reader_linelastic.inputJson["materials"]) {
        mg["matmodel"] = "LinearElasticIsotropic";
    }
    // Create a temporary material manager and solver for the linear elastic case
    MaterialManager<howmany, n_str> *matmanager_lin = createMaterialManager<howmany, n_str>(reader_linelastic);
    Solver<howmany, n_str>          *solver_lin     = createSolver(reader_linelastic, matmanager_lin);
    // Comptue homogenized tangent
    auto C_elastic = solver_lin->get_homogenized_tangent(0.0);
    if (reader.world_rank == 0) {
        printf("║   Homogenized stiffness tensor (C_elastic): \n");
        cout << "# Homogenized tangent: " << endl
             << setprecision(12) << C_elastic << endl
             << endl;
    }
    // Force write C_elastic to results file for reference - location - ms_dataset_name_results/load0/time_step0/homogenized_tangent
    hsize_t dims[2] = {static_cast<hsize_t>(n_str), static_cast<hsize_t>(n_str)};
    char    name[5096];
    snprintf(name, sizeof(name), "%s/load%i/time_step%i/%s", reader.dataset_name, 0, 0, "homogenized_tangent");
    reader.WriteData(C_elastic.data(), name, dims, 2);
    delete solver_lin;
    delete matmanager_lin;
    /*-------------------------------------------------------------------------------------------------------------------------------*/

    reader.ReadMS(howmany);

    for (size_t load_path_idx = 0; load_path_idx < reader.load_cases.size(); ++load_path_idx) {
        MaterialManager<howmany, n_str> *matmanager = createMaterialManager<howmany, n_str>(reader);
        Solver<howmany, n_str>          *solver     = createSolver(reader, matmanager);

        if (reader.world_rank == 0) {
            printf("\n╔════════════════════════════════════════════════════════════ Load case %zu/%zu: %zu time steps ════════════════════════════════════════════════════════════╗\n",
                   load_path_idx + 1, reader.load_cases.size(), reader.load_cases[load_path_idx].n_steps);
        }

        for (size_t time_step_idx = 0; time_step_idx < reader.load_cases[load_path_idx].n_steps; ++time_step_idx) {
            if (reader.world_rank == 0) {
                printf("║   ▶ Time step %zu/%zu (load case %zu/%zu) ◀ \n",
                       time_step_idx + 1, reader.load_cases[load_path_idx].n_steps,
                       load_path_idx + 1, reader.load_cases.size());
            }
            if (reader.load_cases[load_path_idx].mixed) {
                solver->enableMixedBC(reader.load_cases[load_path_idx].mbc, time_step_idx);
            } else {
                const auto &g0 = reader.load_cases[load_path_idx].g0_path[time_step_idx];
                matmanager->set_gradient(g0);
            }
            solver->solve();
            solver->postprocess(reader, load_path_idx, time_step_idx);

            /*-------------------------------------------------------------------------------------------------------------------------------*/
            { // Break the time loop early after meeting the yield criterion (norm(eps - C^-1 : sigma) > threshold)
                VectorXd eps_p = solver->homogenized_strain - C_elastic.ldlt().solve(solver->homogenized_stress);
                if (reader.world_rank == 0) {
                    printf("║    Yield criterion check: ||eps_p|| = %e \n", eps_p.norm());
                }
                // Force write eps_bar_p to results file
                hsize_t dims[1] = {static_cast<hsize_t>(n_str)};
                char    name[5096];
                snprintf(name, sizeof(name), "%s/load%i/time_step%i/%s", reader.dataset_name, load_path_idx, time_step_idx, "eps_bar_p");
                reader.WriteData(eps_p.data(), name, dims, 1);
                // This threshold is somewhat arbitrary and should be tuned. As i briefly suggested to AO, it shoudl be a superset of all future regions of interest.
                if (eps_p.norm() > 0.002) {
                    if (reader.world_rank == 0) {
                        printf("║    Yield criterion met at time step %zu/%zu (load case %zu/%zu), breaking the time loop. \n",
                               time_step_idx + 1, reader.load_cases[load_path_idx].n_steps,
                               load_path_idx + 1, reader.load_cases.size());
                    }
                    break;
                }
            }
            /*-------------------------------------------------------------------------------------------------------------------------------*/

            if (reader.extrapolate_displacement)
                solver->extrapolateDisplacement();
        }
        delete solver;
        delete matmanager;
        if (reader.world_rank == 0) {
            printf("╚══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝\n");
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc > 1 && string(argv[1]) == "--version") {
        cout << "FANS version " << PROJECT_VERSION << endl;
        return 0;
    }

    if (argc != 3) {
        fprintf(stderr, "USAGE: %s [input file basename] [output file basename]\n", argv[0]);
        return 10;
    }

    MPI_Init(NULL, NULL);
    fftw_mpi_init();

    Reader reader{MPI_COMM_WORLD};
    reader.ReadInputFile(argv[1]);
    reader.OpenResultsFile(argv[2]);

    if (reader.problemType == "thermal") {
        runSolver<1, 3>(reader, argv[1]);
    } else if (reader.problemType == "mechanical" && reader.strain_type == "small") {
        runSolver<3, 6>(reader, argv[1]);
    } else if (reader.problemType == "mechanical" && reader.strain_type == "large") {
        runSolver<3, 9>(reader, argv[1]);
    } else {
        throw std::invalid_argument(reader.problemType + " is not a valid problem type");
    }
    reader.CloseResultsFile();

    MPI_Finalize();
    return 0;
}
