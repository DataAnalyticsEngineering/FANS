#include "general.h"
#include "matmodel.h"
#include "setup.h"
#include "solver.h"

// Version
#include "logging.h"
#include "version.h"

template <int howmany, int n_str>
void runSolver(Reader &reader)
{
    reader.ReadMS(howmany);

    for (size_t load_path_idx = 0; load_path_idx < reader.load_cases.size(); ++load_path_idx) {
        MaterialManager<howmany, n_str> *matmanager = createMaterialManager<howmany, n_str>(reader);
        Solver<howmany, n_str>          *solver     = createSolver(reader, matmanager);

        Log::general->info() << "\n╔════════════════════════════════════════════════════════════ Load case "
                             << load_path_idx + 1 << "/" << reader.load_cases.size() << ": " << reader.load_cases[load_path_idx].n_steps
                             << " time steps ════════════════════════════════════════════════════════════╗\n";

        for (size_t time_step_idx = 0; time_step_idx < reader.load_cases[load_path_idx].n_steps; ++time_step_idx) {
            Log::general->info() << "║   ▶ Time step " << time_step_idx + 1 << "/" << reader.load_cases[load_path_idx].n_steps
                                 << " (load case " << load_path_idx + 1 << "/" << reader.load_cases.size() << ") ◀ \n";
            if (reader.load_cases[load_path_idx].mixed) {
                solver->enableMixedBC(reader.load_cases[load_path_idx].mbc, time_step_idx);
            } else {
                const auto &g0 = reader.load_cases[load_path_idx].g0_path[time_step_idx];
                matmanager->set_gradient(g0);
            }
            solver->solve();
            solver->postprocess(reader, load_path_idx, time_step_idx);
        }
        delete solver;
        delete matmanager;
        Log::general->info() << "╚══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝\n";
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

    Reader reader;
    reader.ReadInputFile(argv[1]);
    reader.OpenResultsFile(argv[2]);

    if (reader.problemType == "thermal") {
        runSolver<1, 3>(reader);
    } else if (reader.problemType == "mechanical" && reader.strain_type == "small") {
        runSolver<3, 6>(reader);
    } else if (reader.problemType == "mechanical" && reader.strain_type == "large") {
        runSolver<3, 9>(reader);
    } else {
        throw std::invalid_argument(reader.problemType + " is not a valid problem type");
    }
    reader.CloseResultsFile();

    Log::finalize();
    MPI_Finalize();
    return 0;
}
