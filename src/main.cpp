#include "general.h"
#include "matmodel.h"
#include "setup.h"
#include "solver.h"

// Version
#include "version.h"

template <int howmany>
void runSolver(Reader &reader, const char *output_file_basename)
{
    reader.ReadMS(howmany);

    for (size_t load_path_idx = 0; load_path_idx < reader.load_cases.size(); ++load_path_idx) {
        Matmodel<howmany> *matmodel = createMatmodel<howmany>(reader);
        Solver<howmany>   *solver   = createSolver(reader, matmodel);

        for (size_t time_step_idx = 0; time_step_idx < reader.load_cases[load_path_idx].n_steps; ++time_step_idx) {
            if (reader.load_cases[load_path_idx].mixed) {
                solver->enableMixedBC(reader.load_cases[load_path_idx].mbc, time_step_idx);
            } else {
                const auto &g0 = reader.load_cases[load_path_idx].g0_path[time_step_idx];
                matmodel->setGradient(g0);
            }
            solver->solve();
            solver->postprocess(reader, output_file_basename, load_path_idx, time_step_idx);
        }
        delete solver;
        delete matmodel;
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

    if (reader.problemType == "thermal") {
        runSolver<1>(reader, argv[2]);
    } else if (reader.problemType == "mechanical") {
        runSolver<3>(reader, argv[2]);
    } else {
        throw std::invalid_argument(reader.problemType + " is not a valid problem type");
    }

    MPI_Finalize();
    return 0;
}
