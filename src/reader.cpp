#include "general.h"
#include "reader.h"

#include "H5Cpp.h"
#include "hdf5.h"
#include "stdlib.h"
#include "fftw3-mpi.h"
#include "mpi.h"

#include "H5FDmpi.h"
#include "H5FDmpio.h"

using namespace std;

#include <json.hpp>
using nlohmann::json;
using namespace nlohmann;


void Reader::ComputeVolumeFractions(){

    if (world_rank == 0)
        printf("# Volume fractions\n");

    long vol_frac[n_mat];
    double v_frac[n_mat];
    for(int i=0; i < n_mat; i++){
        vol_frac[i] = 0;
    }
    for(size_t i=0; i < local_n0 * dims[1] * dims[2]; i++) {
        vol_frac[ms[i]]++;
    }
    for(int i=0; i < n_mat; i++){
        long vf;
        MPI_Allreduce(&(vol_frac[i]), &vf, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        v_frac[i] = double(vf) / double(dims[0] * dims[1] * dims[2]);
        if (world_rank == 0)
            printf("# material %4i    vol. frac. %10.4f%%  \n", i, 100.*v_frac[i]);
    }
    // for(auto it = materialProperties.begin(); it != materialProperties.end(); it++){
    //     vector<double> prop = materialProperties[it->first];
    //     double vol_av = 0.;
    //     for(int i=0; i < n_mat; i++){
    //         vol_av += v_frac[i] * prop[i];
    //     }
    //     cout << "# " << it->first << ":    ";  //printf only works for c-strings
    //     printf("minimum : %f, maximum: %f, mean: %f\n", 
    //         *min_element(prop.begin(), prop.end()), *max_element(prop.begin(), prop.end()), vol_av);
    // }
}

string Reader::ReadFileLocations(char fn[]){
    try{
        ifstream i(fn);
        json j;
        i >> j;
        output_path = j["output_path"].get<string>();
        return j["input_file"].get<string>();
    }catch(const std::exception& e){
        fprintf(stderr, "ERROR trying to read input file '%s' for FANS\n", fn );
        exit(10);
    }
}

void Reader :: ReadInputFile(char fn[]){
    try{
    
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    ifstream i(fn);
    json j;
    i >> j;

    strcpy(ms_filename, j["ms_filename"].get<string>().c_str());
    strcpy(ms_datasetname, j["ms_datasetname"].get<string>().c_str());

    L = j["ms_L"].get<vector<double>>();

    TOL = j["TOL"].get<double>();
    n_it = j["n_it"].get<int>();
    if (j.contains("g0")){
        g0 = j["g0"].get<vector<double>>();
    }

    problemType = j["problem_type"].get<string>();
    matmodel = j["matmodel"].get<string>();
    method = j["method"].get<string>();
    
    json j_mat = j["material_properties"];
    for (auto it = j_mat.begin(); it != j_mat.end(); ++it){
        
        materialProperties[it.key()] = it.value().get<vector<double>>();
        n_mat = materialProperties[it.key()].size();

        if (world_rank == 0){
        cout << "# " << it.key() << ":\t ";
        for(double d : materialProperties[it.key()]){
            printf("   %10.5f", d);
        }
        printf("\n");
        }
    }
    if (world_rank == 0){
        printf("# microstructure: \t '%s'\n", ms_filename);
        printf("# FANS Tolerance: \t %10.5e\n# Max iterations: \t %6i\n", TOL, n_it);
    }

    // Read the results_to_write field
    resultsToWrite = j["results"].get<vector<string>>();

    }catch(const std::exception& e){
        fprintf(stderr, "ERROR trying to read input file '%s' for FANS\n", fn );
        exit(10);
    }
}

//    printf("# Macro-scale Gradient - (");
//    for (const auto& number : g0) {
//        printf("%10.5f ", number);}
//    printf(")\n");

void Reader::safe_create_group( hid_t file, const char * const name )
{
    // no leading '/' --> exit
    const char DELIMITER = '/';
    if( name[0] != DELIMITER )
        return;
    
    // copy name to buffer
    char buffer[4096];
    strcpy(buffer,name);
    char * str = buffer;
    str=strchr(str + 1,DELIMITER);
    while( str != NULL )
    {
        // while another / character is found 
        long int l=str-buffer; // length of substring
        buffer[l] = '\0'; // temporary 'end of string'

        // safely create the group if needed
        hid_t group;

        /* Save old error handler */
        herr_t (*old_func)(hid_t, void*);
        void *old_client_data;
        H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
        /* Turn off error handling */
        H5Eset_auto(H5E_DEFAULT, NULL, NULL);

        group = H5Gopen(file, buffer, H5P_DEFAULT);
        /* Restore previous error handler */
        H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

        if(group < 0){
            group = H5Gcreate(file, buffer, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        }
        H5Gclose(group);

        buffer[l] = DELIMITER; // restore original string

        str=strchr(str + 1, DELIMITER); // find next delimiter
    }
}


void Reader :: ReadMS(int hm){
	
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hid_t       data_type;
    hsize_t     _dims[3];                 /* dataset dimensions */
    hsize_t	    count[3];	          /* hyperslab selection parameters */
    hsize_t	    offset[3];
    hid_t	    plist_id;                 /* property list identifier */
    herr_t	    status;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Info info  = MPI_INFO_NULL;

    //Set up file access property list with parallel I/O access
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);        // "set File Access Property List"

    //Open the file collectively and release property list identifier.
    file_id = H5Fopen(ms_filename, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    //Create property list for collective dataset write.
    //plist_id = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);   // "set Data Transfer Property List" (x means transfer)
    plist_id = H5P_DEFAULT;


    // read physical dimensions: Not needed. Added ms_L to json file
//    double* L = FANS_malloc<double>(3);
//    char name_L[5096];
//    sprintf(name_L,"%s_L", ms_datasetname);
//    dset_id = H5Dopen2(file_id, name_L, plist_id);
//    data_type = H5Dget_type(dset_id);
//    status = H5Dread(dset_id, data_type, H5S_ALL, H5S_ALL, plist_id, L);

    dset_id = H5Dopen2(file_id, ms_datasetname, plist_id);

    hid_t dspace = H5Dget_space(dset_id);
    int rank = H5Sget_simple_extent_dims(dspace, _dims, NULL);
    data_type = H5T_NATIVE_UCHAR;   //could also use H5Dget_type(dset_id)

    dims.resize(3);
    dims[0] = _dims[0];
    dims[1] = _dims[1];
    dims[2] = _dims[2];

    l_e.resize(3);
    l_e[0] = L[0] / double(dims[0]);
    l_e[1] = L[1] / double(dims[1]);
    l_e[2] = L[2] / double(dims[2]);

    if (world_rank == 0){
        printf("# grid size set to [%i x %i x %i] --> %i voxels \nMicrostructure length: [%3.6f x %3.6f x %3.6f]\n", dims[0], dims[1], dims[2], dims[0]*dims[1]*dims[2], L[0], L[1], L[2] );
        // if(dims[0] % 2 != 0)	fprintf(stderr, "[ FANS3D_Grid ] WARNING: n_x is not a multiple of 2\n");
        // if(dims[1] % 2 != 0)	fprintf(stderr, "[ FANS3D_Grid ] WARNING: n_y is not a multiple of 2\n");
        // if(dims[2] % 2 != 0)	fprintf(stderr, "[ FANS3D_Grid ] WARNING: n_z is not a multiple of 2\n");
        printf("Voxel length: [%1.8f, %1.8f, %1.8f]\n", l_e[0], l_e[1],l_e[2]);
    }

    const ptrdiff_t n[3]  = {dims[0], dims[1], dims[2] / 2 + 1};
    ptrdiff_t block0 = FFTW_MPI_DEFAULT_BLOCK;
    ptrdiff_t block1 = FFTW_MPI_DEFAULT_BLOCK;


    // see https://fftw.org/doc/Basic-and-advanced-distribution-interfaces.html
    // and https://www.fftw.org/fftw3_doc/Transposed-distributions.html
    // on there it is recommended to use one of fftw's allocation functions "to ensure optimal alignment"

    /* there is no documentation for this method, so here is the signature from "fftw3-mpi.h"
    FFTW_EXTERN ptrdiff_t XM(local_size_many_transposed)	\
     (int rnk, const ptrdiff_t *n, ptrdiff_t howmany,		\
      ptrdiff_t block0, ptrdiff_t block1, MPI_Comm comm,	\
      ptrdiff_t *local_n0, ptrdiff_t *local_0_start,		\
      ptrdiff_t *local_n1, ptrdiff_t *local_1_start);		\
    */

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n, hm, block0, block1, MPI_COMM_WORLD, &local_n0, &local_0_start, &local_n1, &local_1_start);
    

    //Each process defines a dataset in memory which reads a hyperslab from the file
    count[0] = local_n0;
    count[1] = dims[1];
    count[2] = dims[2];
    memspace = H5Screate_simple(rank, count, NULL);

    //Select hyperslab in the file.
    offset[0] = local_0_start;
    offset[1] = 0;
    offset[2] = 0;
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    ms = FANS_malloc<unsigned char>(count[0] * count[1] * count[2]);
    status = H5Dread(dset_id, data_type, memspace, filespace, plist_id, this->ms);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Fclose(file_id);
    //H5Tclose(data_type);
}






// The code above is based on this example: Hyperslab_by_row.c


// /*  
//  *  This example writes data to the HDF5 file by rows.
//  *  Number of processes is assumed to be 1 or multiples of 2 (up to 8)
//  */
 
// #include "hdf5.h"
// #include "stdlib.h"

// #define H5FILE_NAME     "SDS_row.h5"
// #define DATASETNAME 	"IntArray" 
// #define NX     8                      /* dataset dimensions */
// #define NY     5 
// #define RANK   2

// int
// main (int argc, char **argv)
// {
//     /*
//      * HDF5 APIs definitions
//      */ 	
//     hid_t       file_id, dset_id;         /* file and dataset identifiers */
//     hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
//     hsize_t     dimsf[2];                 /* dataset dimensions */
//     int         *data;                    /* pointer to data buffer to write */
//     hsize_t	count[2];	          /* hyperslab selection parameters */
//     hsize_t	offset[2];
//     hid_t	plist_id;                 /* property list identifier */
//     int         i;
//     herr_t	status;

//     /*
//      * MPI variables
//      */
//     int mpi_size, mpi_rank;
//     MPI_Comm comm  = MPI_COMM_WORLD;
//     MPI_Info info  = MPI_INFO_NULL;

//     /*
//      * Initialize MPI
//      */
//     MPI_Init(&argc, &argv);
//     MPI_Comm_size(comm, &mpi_size);
//     MPI_Comm_rank(comm, &mpi_rank);  
 
//     /* 
//      * Set up file access property list with parallel I/O access
//      */
//      plist_id = H5Pcreate(H5P_FILE_ACCESS);
//      H5Pset_fapl_mpio(plist_id, comm, info);

//     /*
//      * Create a new file collectively and release property list identifier.
//      */
//     file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
//     H5Pclose(plist_id);
   

//     /*
//      * Create the dataspace for the dataset.
//      */
//     dimsf[0] = NX;
//     dimsf[1] = NY;
//     filespace = H5Screate_simple(RANK, dimsf, NULL); 

//     /*
//      * Create the dataset with default properties and close filespace.
//      */
//     dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_INT, filespace,
// 			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     H5Sclose(filespace);

//     /* 
//      * Each process defines dataset in memory and writes it to the hyperslab
//      * in the file.
//      */
//     count[0] = dimsf[0]/mpi_size;
//     count[1] = dimsf[1];
//     offset[0] = mpi_rank * count[0];
//     offset[1] = 0;
//     memspace = H5Screate_simple(RANK, count, NULL);

//     /*
//      * Select hyperslab in the file.
//      */
//     filespace = H5Dget_space(dset_id);
//     H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

//     /*
//      * Initialize data buffer 
//      */
//     data = (int *) malloc(sizeof(int)*count[0]*count[1]);
//     for (i=0; i < count[0]*count[1]; i++) {
//         data[i] = mpi_rank + 10;
//     }

//     /*
//      * Create property list for collective dataset write.
//      */
//     plist_id = H5Pcreate(H5P_DATASET_XFER);
//     H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
//     status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
// 		      plist_id, data);
//     free(data);

//     /*
//      * Close/release resources.
//      */
//     H5Dclose(dset_id);
//     H5Sclose(filespace);
//     H5Sclose(memspace);
//     H5Pclose(plist_id);
//     H5Fclose(file_id);
 
//     MPI_Finalize();

//     return 0;
// }     