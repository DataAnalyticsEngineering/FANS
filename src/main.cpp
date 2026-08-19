#include "general.h"
#include "matmodel.h"
#include "setup.h"
#include "solver.h"

// Version
#include "version.h"

template <int howmany, int n_str>
void runSolver(Reader &reader, char input_fn[])
{
    reader.ReadMS(howmany);

    for (size_t load_path_idx = 0; load_path_idx < reader.load_cases.size(); ++load_path_idx) {
        MaterialManager<howmany, n_str> *matmanager = createMaterialManager<howmany, n_str>(reader);
        Solver<howmany, n_str>          *solver     = createSolver(reader, matmanager);

        Log::logger().info("\n╔════════════════════════════════════════════════════════════ Load case {}/{}: {} time steps ════════════════════════════════════════════════════════════╗",
                           load_path_idx + 1, reader.load_cases.size(), reader.load_cases[load_path_idx].n_steps);

        for (size_t time_step_idx = 0; time_step_idx < reader.load_cases[load_path_idx].n_steps; ++time_step_idx) {
            Log::logger().info("║   ▶ Time step {}/{} (load case {}/{}) ◀ ",
                               time_step_idx + 1, reader.load_cases[load_path_idx].n_steps,
                               load_path_idx + 1, reader.load_cases.size());
            if (reader.load_cases[load_path_idx].mixed) {
                solver->enableMixedBC(reader.load_cases[load_path_idx].mbc, time_step_idx);
            } else {
                const auto &g0 = reader.load_cases[load_path_idx].g0_path[time_step_idx];
                matmanager->set_gradient(g0);
            }
            solver->solve();
            solver->postprocess(reader, load_path_idx, time_step_idx);
            if (reader.extrapolate_displacement)
                solver->extrapolateDisplacement();
        }
        delete solver;
        delete matmanager;
        Log::logger().info("╚══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝");
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(nullptr, nullptr);
    Log::init(Log::parse_options(argc, argv));

    if (argc > 1 && string(argv[1]) == "--version") {
        Log::logger().info("FANS version {}", PROJECT_VERSION);
        MPI_Finalize();
        return 0;
    }

    if (argc < 3) {
        Log::logger().error("USAGE: {} [input file] [output file] [--log-level LEVEL]", argv[0]);
        MPI_Finalize();
        return 10;
    }

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
