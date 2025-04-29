#ifndef GBDIFFUSION_H
#define GBDIFFUSION_H

#include "matmodel.h"
#include <Eigen/StdVector> // For Eigen's aligned_allocator

// Define protype for ReadAuxData4D here
void ReadAuxData4D(
    Reader            &r,           // Reader with existing 3D decomposition info
    const std::string &filename,    // HDF5 file name (not using r.ms_filename)
    const std::string &datasetName, // 4D dataset path inside the file (not using r.ms_datasetname)
    double           *&auxData4D    // Output pointer for the local 4D chunk (type double)
);

class GBDiffusion : public ThermalModel {
  public:
    GBDiffusion(Reader &reader)
        : ThermalModel(reader.l_e)
    {
        try {
            D_bulk             = reader.materialProperties["D_bulk"].get<double>();
            D_par              = reader.materialProperties["D_par"].get<double>();
            D_perp             = reader.materialProperties["D_perp"].get<double>();
            normal_filename    = reader.materialProperties["GBNormal_filename"].get<std::string>();
            normal_datasetname = reader.materialProperties["GBNormal_datasetname"].get<std::string>();

        } catch (const std::exception &e) {
            throw std::runtime_error("Missing material properties for the requested material model.");
        }

        // Read the GB normals from the HDF5 file
        ReadAuxData4D(reader, normal_filename, normal_datasetname, GBnormals);

        kapparef_mat = D_bulk * Matrix3d::Identity();
    }
    ~GBDiffusion()
    {
        if (GBnormals != nullptr) {
            delete[] GBnormals; // Clean up dynamically allocated memory
            GBnormals = nullptr;
        }
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        if (mat_index > -1) {
            // Bulk is isotropic
            sigma.block<3, 1>(i, 0) = D_bulk * eps.block<3, 1>(i, 0);
        } else if (mat_index == -1) {
            const ptrdiff_t base_idx = 3 * element_idx;
            double          nx       = GBnormals[base_idx];
            double          ny       = GBnormals[base_idx + 1];
            double          nz       = GBnormals[base_idx + 2];

            // Pre-compute normalization factor
            double inv_norm = 1.0 / std::sqrt(nx * nx + ny * ny + nz * nz);
            nx *= inv_norm;
            ny *= inv_norm;
            nz *= inv_norm;

            // Pre-compute products for the projector matrix (N⊗N)
            double nxnx = nx * nx;
            double nxny = nx * ny;
            double nxnz = nx * nz;
            double nyny = ny * ny;
            double nynz = ny * nz;
            double nznz = nz * nz;

            // Pre-compute coefficients
            double d_diff = D_par - D_perp;

            // Cache epsilon values to avoid repeated memory access
            double ex = eps(i, 0);
            double ey = eps(i + 1, 0);
            double ez = eps(i + 2, 0);

            // Calculate directly without constructing full matrices
            sigma(i, 0)     = D_par * ex - d_diff * (nxnx * ex + nxny * ey + nxnz * ez);
            sigma(i + 1, 0) = D_par * ey - d_diff * (nxny * ex + nyny * ey + nynz * ez);
            sigma(i + 2, 0) = D_par * ez - d_diff * (nxnz * ex + nynz * ey + nznz * ez);
        } else {
            throw std::runtime_error("GBDiffusion: Unknown material index");
        }
    }

    void postprocess(Solver<1> &solver, Reader &reader, const char *resultsFileName, int load_idx, int time_idx) override
    {

        if (find(reader.resultsToWrite.begin(), reader.resultsToWrite.end(), "GBnormals") != reader.resultsToWrite.end()) {
            for (int i = 0; i < solver.world_size; ++i) {
                if (i == solver.world_rank) {
                    char name[5096];
                    sprintf(name, "%s/load%i/time_step%i/GBnormals", reader.ms_datasetname, load_idx, time_idx);
                    reader.WriteSlab<double>(GBnormals, 3, resultsFileName, name);
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
    }

  private:
    double D_bulk;
    double D_par;
    double D_perp;
    string normal_filename;
    string normal_datasetname;

    double  *GBnormals = nullptr;
    Vector3d N;    // GB normal
    Matrix3d D_GB; // GB diffusion tensor
};

void ReadAuxData4D(
    Reader            &r,           // Reader with existing 3D decomposition info
    const std::string &filename,    // HDF5 file name (not using r.ms_filename)
    const std::string &datasetName, // 4D dataset path inside the file (not using r.ms_datasetname)
    double           *&auxData4D    // Output pointer for the local 4D chunk (type double)
)
{
    //--------------------------------------------------------------------------
    // 1. Open the HDF5 file (parallel or serial)
    //--------------------------------------------------------------------------
    MPI_Info info     = MPI_INFO_NULL;
    hid_t    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    // If you want parallel I/O:
    // H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    //--------------------------------------------------------------------------
    // 2. Open the 4D dataset
    //--------------------------------------------------------------------------
    // If you want collective reads:
    // plist_id = H5Pcreate(H5P_DATASET_XFER);
    // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    // Minimal approach:
    plist_id = H5P_DEFAULT;

    hid_t dset_id   = H5Dopen2(file_id, datasetName.c_str(), plist_id);
    hid_t filespace = H5Dget_space(dset_id);

    //--------------------------------------------------------------------------
    // 3. Confirm shape is Nx x Ny x Nz x k and Nx,Ny,Nz match r.dims[]
    //--------------------------------------------------------------------------
    hsize_t _dims[4];
    int     rank = H5Sget_simple_extent_dims(filespace, _dims, nullptr);
    if (rank != 4) {
        throw std::runtime_error("ReadAuxData4D: dataset is not 4D (rank=" + std::to_string(rank) + ")");
    }

    /* Verify the three spatial axes match r.dims in the appropriate order */
    const bool is_zyx = r.is_zyx; // ← flag set in Reader::ReadMS
    if (!is_zyx) {                /* on-disk order  X Y Z k */
        if (_dims[0] != (hsize_t) r.dims[0] ||
            _dims[1] != (hsize_t) r.dims[1] ||
            _dims[2] != (hsize_t) r.dims[2])
            throw std::runtime_error("ReadAuxData4D: XYZ dims do not match Reader grid.");
    } else {                                   /* on-disk order  Z Y X k */
        if (_dims[2] != (hsize_t) r.dims[0] || /* X */
            _dims[1] != (hsize_t) r.dims[1] || /* Y */
            _dims[0] != (hsize_t) r.dims[2])   /* Z */
            throw std::runtime_error("ReadAuxData4D: ZYX dims do not match Reader grid.");
    }
    const hsize_t kDim = _dims[3]; /* size of 4th axis */

    if (r.world_rank == 0) {
        printf("Reading 4-D dataset [%llu × %llu × %llu × %llu] (%s layout)…\n",
               (unsigned long long) _dims[0],
               (unsigned long long) _dims[1],
               (unsigned long long) _dims[2],
               (unsigned long long) _dims[3],
               is_zyx ? "ZYXk" : "XYZk");
    }

    /* ------------------------------------------------------------------ */
    /* 4. build FILE hyperslab (always slice along physical X)            */
    /* ------------------------------------------------------------------ */
    hsize_t fcount[4];
    hsize_t foffset[4];
    if (!is_zyx) { /* file dims:  X Y Z k   (slice file-dim 0) */
        fcount[0]  = static_cast<hsize_t>(r.local_n0);
        fcount[1]  = _dims[1]; /* Ny */
        fcount[2]  = _dims[2]; /* Nz */
        fcount[3]  = kDim;
        foffset[0] = static_cast<hsize_t>(r.local_0_start);
        foffset[1] = 0;
        foffset[2] = 0;
        foffset[3] = 0;
    } else {                                           /* file dims:  Z Y X k   (slice file-dim 2) */
        fcount[0]  = _dims[0];                         /* Nz */
        fcount[1]  = _dims[1];                         /* Ny */
        fcount[2]  = static_cast<hsize_t>(r.local_n0); /* Nx-slab */
        fcount[3]  = kDim;
        foffset[0] = 0;
        foffset[1] = 0;
        foffset[2] = static_cast<hsize_t>(r.local_0_start);
        foffset[3] = 0;
    }
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, foffset, nullptr, fcount, nullptr);

    /* ------------------------------------------------------------------ */
    /* 5. build MEMORY dataspace that matches the FILE slab exactly       */
    /*    (so HDF5 performs a pure copy, no hidden re-ordering)           */
    /* ------------------------------------------------------------------ */
    hid_t memspace = H5Screate_simple(4, fcount, nullptr);

    /* ------------------------------------------------------------------ */
    /* 6. read – with a transpose step only if ZYX was stored             */
    /* ------------------------------------------------------------------ */
    const size_t NxLoc = static_cast<size_t>(r.local_n0);
    const size_t Ny    = static_cast<size_t>(r.dims[1]);
    const size_t Nz    = static_cast<size_t>(r.dims[2]);

    size_t slabElems = static_cast<size_t>(fcount[0]) *
                       static_cast<size_t>(fcount[1]) *
                       static_cast<size_t>(fcount[2]) *
                       static_cast<size_t>(fcount[3]);

    if (!is_zyx) {
        /* XYZk : we can read straight into the final buffer */
        auxData4D = new double[slabElems];

        herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE,
                                memspace, filespace, plist_id, auxData4D);
        if (status < 0)
            throw std::runtime_error("ReadAuxData4D: H5Dread failed (XYZ)");
    } else {
        /* ZYXk : read into temp buffer shaped [Z][Y][X][k], then transpose
         *        to final logical order [X][Y][Z][k].                    */
        double *tmp = new double[slabElems];

        herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE,
                                memspace, filespace, plist_id, tmp);
        if (status < 0)
            throw std::runtime_error("ReadAuxData4D: H5Dread failed (ZYX)");

        auxData4D = new double[slabElems];

        /* temp indices:  idx_t = (((z*Ny)+y)*NxLoc + x)*k + d
         * final indices: idx_f = (((x*Ny)+y)*Nz    + z)*k + d                */
        for (size_t z = 0; z < Nz; ++z)
            for (size_t y = 0; y < Ny; ++y)
                for (size_t x = 0; x < NxLoc; ++x) {
                    size_t base_t = (((z * Ny) + y) * NxLoc + x) * kDim;
                    size_t base_f = (((x * Ny) + y) * Nz + z) * kDim;
                    std::memcpy(&auxData4D[base_f],
                                &tmp[base_t],
                                kDim * sizeof(double));
                }
        delete[] tmp;
    }

    /* ------------------------------------------------------------------ */
    /* 7. tidy up                                                         */
    /* ------------------------------------------------------------------ */
    H5Sclose(memspace);
    H5Dclose(dset_id);
    H5Fclose(file_id);

    MPI_Barrier(MPI_COMM_WORLD);

    if (r.world_rank == 0)
        printf("ReadAuxData4D: local slab read complete.\n");
}

#endif // GBDIFFUSION_H
