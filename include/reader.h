
#ifndef READER_H
#define READER_H

#include <vector>
#include <string>
#include <map>

using namespace std;

class Reader{
    public:
        // contents of input file:
        char ms_filename[4096];                                                  // Name of Micro-structure hdf5 file
        char ms_datasetname[4096];                                                  // Absolute path of Micro-structure in hdf5 file
        int n_mat;
        map<string, vector<double>> materialProperties;
        double TOL;
        int n_it;
        vector<double> g0;
        string problemType;
        string matmodel;
        string method;

        vector<string> resultsToWrite;

        // contents of microstructure file:
        vector<int> dims;
        vector<double> l_e;
        vector<double> L;
        unsigned char *ms;                                                            // Micro-structure Binary

        int world_rank;
        int world_size;

        ptrdiff_t alloc_local;
        ptrdiff_t local_n0;
        ptrdiff_t local_0_start;        // this is the x-value of the start point, not the index in the array
        ptrdiff_t local_n1;
        ptrdiff_t local_1_start;

        string output_path;

        // void Setup(ptrdiff_t howmany);
        string ReadFileLocations(char fn[]);
        void ReadInputFile(char fn[]);
        void ReadMS(int hm);
        void ComputeVolumeFractions();
        // void ReadHDF5(char file_name[], char dset_name[]);
        void safe_create_group( hid_t file, const char * const name );
        
        template<typename T>
        void WriteSlab(T *data, int _howmany, const char* file_name, const char* dset_name);

        template<typename T>
        void WriteData(T *data, const char* file_name, const char* dset_name, hsize_t* dims, int rank);        
};

template<typename T>
void Reader::WriteData(T *data, const char* file_name, const char* dset_name, hsize_t* dims, int rank)
{
    hid_t data_type;
    if (std::is_same<T, double>::value) {
        data_type = H5T_NATIVE_DOUBLE;
    } else if (std::is_same<T, unsigned char>::value) {
        data_type = H5T_NATIVE_UCHAR;
    } else {
        throw std::invalid_argument("Conversion of this data type to H5 data type not yet implemented");
    }

    hid_t	    plist_id;                
    hid_t       file_id;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    //TODO: refactor this into a general error handling method
    /* Save old error handler */
    herr_t (*old_func)(hid_t, void*);
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
    /* Turn off error handling */
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    file_id = H5Fopen(file_name, H5F_ACC_RDWR, plist_id);
    /* Restore previous error handler */
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

    if(file_id < 0){
        file_id = H5Fcreate(file_name, H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
    }
    H5Pclose(plist_id);

    // Ensure all groups in the path are created
    safe_create_group(file_id, dset_name);

    // Check if the dataset already exists
    if (H5Lexists(file_id, dset_name, H5P_DEFAULT) > 0) {
        // If it exists, delete it
        H5Ldelete(file_id, dset_name, H5P_DEFAULT);
    }

    // Create the data space for the dataset
    hid_t dataspace_id = H5Screate_simple(rank, dims, NULL);
    if (dataspace_id < 0) {
        H5Fclose(file_id);
        throw std::runtime_error("Error creating dataspace");
    }

    // Create the dataset with default properties
    hid_t dataset_id = H5Dcreate2(file_id, dset_name, data_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
        throw std::runtime_error("Error creating dataset");
    }

    // Write the data to the dataset
    if (H5Dwrite(dataset_id, data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0) {
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
        throw std::runtime_error("Error writing data to dataset");
    }

    // Close the dataset and the file
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Fclose(file_id);
}

// this function has to be here because of the template
template<typename T>
void Reader::WriteSlab(T *data, int _howmany, const char* file_name, const char* dset_name)
{

    hid_t data_type;
    if(std::is_same<T, double>()){
        data_type = H5T_NATIVE_DOUBLE;
    } else if(std::is_same<T, unsigned char>()){
        data_type = H5T_NATIVE_UCHAR;
    } else {
        throw std::invalid_argument("Conversion of this data type to H5 data type not yet implemented");
    }
	
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hid_t	    plist_id;                 /* property list identifier */
    herr_t	    status;
    int         rank = 4;
    hsize_t     dimsf[rank];                 /* dataset dimensions */
    hsize_t	    count[rank];	          /* hyperslab selection parameters */
    hsize_t	    offset[rank]; 
 
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    // H5Pset_fapl_mpio(plist_id, comm, info);

    //TODO: refactor this into a general error handling method
    /* Save old error handler */
    herr_t (*old_func)(hid_t, void*);
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
    /* Turn off error handling */
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    file_id = H5Fopen(file_name, H5F_ACC_RDWR, plist_id);
    /* Restore previous error handler */
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

    if(file_id < 0){
        file_id = H5Fcreate(file_name, H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
    }
    H5Pclose(plist_id);

    dimsf[0] = this->dims[0];
    dimsf[1] = this->dims[1];
    dimsf[2] = this->dims[2];
    dimsf[3] = _howmany;

    filespace = H5Screate_simple(rank, dimsf, NULL); 

    //dset_name = "test123";
    safe_create_group(file_id, dset_name);

    /* Save old error handler */
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
    /* Turn off error handling */
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    /* Restore previous error handler */
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

    if(dset_id < 0){
        dset_id = H5Dcreate(file_id, dset_name, data_type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // hid_t dcpl_id;
        // // Create the dataset with compression
        // dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        // H5Pset_chunk(dcpl_id, rank, dimsf);
        // H5Pset_deflate(dcpl_id, 9); // Compression level 9, 6 is best for size and speed
        // dset_id = H5Dcreate(file_id, dset_name, data_type, filespace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        // H5Pclose(dcpl_id);
    }
    
    count[0] = local_n0;
    count[1] = dimsf[1];
    count[2] = dimsf[2];
    count[3] = dimsf[3];
    offset[0] = local_0_start;
    offset[1] = 0;
    offset[2] = 0;
    offset[3] = 0;
    memspace = H5Screate_simple(rank, count, NULL);

    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
    status = H5Dwrite(dset_id, data_type, memspace, filespace, plist_id, data);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Fclose(file_id);
}


#endif