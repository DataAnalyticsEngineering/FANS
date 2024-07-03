#include "general.h"
#include "matmodel.h"
#include "solver.h"
#include "solverFP.h"
#include "solverCG.h"

int main( int argc, char ** argv ) {
    if( argc != 2 ) {
        fprintf( stderr, "USAGE: %s [input file basename] \n", argv[0] );
        return 10;
    }

    MPI_Init(NULL, NULL);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // initialize fftw mpi
    //
	fftw_mpi_init();

    Reader reader;

    // only necessary for sequential reading
    // for (int i = 0; i < world_size; i++){
    // 	if(i == world_rank){
    //     	reader.ReadInputFile(argv[1]);
    //   	}
    //   	MPI_Barrier(MPI_COMM_WORLD);
    // }
    reader.ReadInputFile(argv[1]);
    MPI_Barrier(MPI_COMM_WORLD);

    if(reader.problemType == "thermal"){

        for (int i = 0; i < world_size; i++){
    	    if(i == world_rank){
            	reader.ReadMS(1);
      	    }
      	    MPI_Barrier(MPI_COMM_WORLD);
        }
        reader.ComputeVolumeFractions();

        Matmodel<1>* matmodel;
        Solver<1>* solver;

        if(reader.matmodel == "ThermalLinear"){
            matmodel = new ThermalLinear(reader.l_e, reader.materialProperties);
        }else{
            throw invalid_argument(reader.matmodel + " is not a valid matmodel");
        }

        if(reader.method == "fp"){
            solver = new SolverFP<1>(reader, matmodel);
        }else if(reader.method == "cg"){
            solver = new SolverCG<1>(reader, matmodel);
        }else{
            throw invalid_argument(reader.method + " is not a valid method");
        }


        vector<double> g0_all = reader.g0;
        uint n_loads = g0_all.size() / matmodel->n_str;
        if (g0_all.size() % matmodel->n_str != 0) throw invalid_argument("Invalid length of loading g0");

        vector<double> g0(matmodel->n_str);
        for(int i_load = 0; i_load < n_loads; i_load++){
            for (int i = 0; i < matmodel->n_str; ++i) {
                g0[i] = g0_all[i_load*matmodel->n_str + i];
            }

            matmodel->setGradient(g0);
            solver->solve();
	    	solver->postprocess(reader, "results.h5", i_load);
        }
    } else if(reader.problemType == "mechanical"){

        for (int i = 0; i < world_size; i++){
    	    if(i == world_rank){
            	reader.ReadMS(3);
      	    }
      	    MPI_Barrier(MPI_COMM_WORLD);
        }
        reader.ComputeVolumeFractions();

        Matmodel<3>* matmodel;
        Solver<3>* solver;

        if(reader.matmodel == "MechLinear"){
            matmodel = new MechLinear(reader.l_e, reader.materialProperties);
        }else if(reader.matmodel == "HyperElastic"){
            matmodel = new HyperElastic(reader.l_e, reader.materialProperties);
        }else{
            throw invalid_argument(reader.matmodel + " is not a valid matmodel");
        }

        if(reader.method == "fp"){
            solver = new SolverFP<3>(reader, matmodel);
        }else if(reader.method == "cg"){
            solver = new SolverCG<3>(reader, matmodel);
        }else{
            throw invalid_argument(reader.method + " is not a valid method");
        }

        vector<double> g0_all = reader.g0;
        uint n_loads = g0_all.size() / matmodel->n_str;
        if (g0_all.size() % matmodel->n_str != 0) throw invalid_argument("Invalid length of loading g0");

        vector<double> g0(matmodel->n_str);
        for(int i_load = 0; i_load < n_loads; i_load++){
            for (int i = 0; i < matmodel->n_str; ++i) {
                g0[i] = g0_all[i_load*matmodel->n_str + i];
            }
            matmodel->setGradient(g0);
            solver->solve();
	    	solver->postprocess(reader, "results.h5", i_load);
        }
    }

    MPI_Finalize();
    return 0;
}