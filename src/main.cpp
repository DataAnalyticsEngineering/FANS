#include "general.h"
#include "matmodel.h"
#include "setup.h"
#include "solver.h"

template <int howmany>
void runSolver(Reader &reader, const char *output_file_basename)
{
    reader.ReadMS(howmany);

    for (size_t load_path_idx = 0; load_path_idx < reader.g0.size(); ++load_path_idx) {
        Matmodel<howmany> *matmodel = createMatmodel<howmany>(reader);
        Solver<howmany>   *solver   = createSolver(reader, matmodel);

        const auto &load_path = reader.g0[load_path_idx];
        for (size_t time_step_idx = 0; time_step_idx < load_path.size(); ++time_step_idx) {
            const auto &g0 = load_path[time_step_idx];

            if (g0.size() != matmodel->n_str) {
                throw std::invalid_argument("Invalid length of loading g0");
            }

            matmodel->setGradient(g0);
            solver->solve();
            solver->postprocess(reader, output_file_basename, load_path_idx, time_step_idx);
        }

        delete solver;
        delete matmodel;
    }
}

int main(int argc, char *argv[])
{
    if (argc != 3) {
        fprintf(stderr, "USAGE: %s [input file basename] [output file basename]\n", argv[0]);
        return 10;
    }

    MPI_Init(NULL, NULL);
    fftw_mpi_init();

    // std::string filePath = "/home/torben/Documents/git/FANS/test/input_files/sphere_ThermalLinear.json";
    const char* actual_path = "/home/torben/Documents/git/FANS/test/input_files/sphere_ThermalLinear.json";  // Replace with your actual path
  int path_length = strlen(actual_path) + 1;  // Add 1 for null terminator

  char* temp_path = new char[path_length];  // Allocate memory for the path
  strcpy(temp_path, actual_path);
    printf("filePath: %s\n", temp_path);
    Reader reader;
    reader.ReadInputFile(temp_path);

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
