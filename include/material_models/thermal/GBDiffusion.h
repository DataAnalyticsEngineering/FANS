#ifndef GBDIFFUSION_H
#define GBDIFFUSION_H

#include "matmodel.h"
#include <Eigen/StdVector> // For Eigen's aligned_allocator

/**
 * @class GBDiffusion
 * @brief Material model for grain boundary diffusion in polycrystals
 *
 * This model implements diffusion in a polycrystalline material, differentiating between
 * bulk crystal diffusion and grain boundary diffusion (transversely isotropic).
 * The grains are characterized by a unit quaternion defining the local orientation and can
 * have diffeerent diffusion properties along the three main axes.
 * The grain boundaries are characterized by their normal vectors and can have different
 * diffusion properties parallel and perpendicular to the boundary plane.
 *
 * The model extends both ThermalModel and LinearModel<1, 3> to provide a linear diffusion
 * formulation that can be used in a thermal-like solver in FANS.
 *
 * @details The model:
 *   - Reads microstructure data containing grain boundaries from HDF5 files
 *   - Supports uniform or material-specific diffusivity values
 *   - Handles bulk regions with arbitrary rotations of reference diffusion tensor (D_bulk_11, D_bulk_12, ..., D_bulk_33)
 *   - Handles grain boundaries with transversely isotropic diffusion (D_par, D_perp)
 *   - Provides visualization of grain orientations and grain boundary normals in post-processing
 *
 * Required material parameters in JSON format:
 *   - material_uniformity: Boolean flag for uniform grain and GB properties
 *
 *   When material_uniformity is true (uniform properties):
 *   {
 *     "material_unformity": true,
 *     "D_bulk_11": 10.0,    // Diffusion coefficient (1,1) for all crystals in crystal system
 *     "D_bulk_12": 0.0,     // Diffusion coefficient (1,2) for all crystals in crystal system
 *     "D_bulk_13": 0.0,     // Diffusion coefficient (1,3) for all crystals in crystal system
 *     "D_bulk_22": 5.0,     // Diffusion coefficient (2,2) for all crystals in crystal system
 *     "D_bulk_23": 0.0,     // Diffusion coefficient (2,3) for all crystals in crystal system
 *     "D_bulk_33": 0.0,     // Diffusion coefficient (3,3) for all crystals in crystal system
 *     "D_par": 2.0,         // Diffusion coefficient parallel to the grain boundary for all GBs
 *     "D_perp": 0.5         // Diffusion coefficient perpendicular to the grain boundary for all GBs
 *   }
 *
 *   When material_unformity is false (tag-specific properties):
 *   {
 *     "material_unformity": false,
 *     "D_bulk_11": [...],    // Array of length (num_crystals + num_GB elements, but D_bulk_11 is only used for crystals (0 to num_crystals)
 *     "D_bulk_12": [...],    // Array of length (num_crystals + num_GB elements, but D_bulk_12 is only used for crystals (0 to num_crystals)
 *     "D_bulk_13": [...],    // Array of length (num_crystals + num_GB elements, but D_bulk_13 is only used for crystals (0 to num_crystals)
 *     "D_bulk_22": [...],    // Array of length (num_crystals + num_GB elements, but D_bulk_22 is only used for crystals (0 to num_crystals)
 *     "D_bulk_23": [...],    // Array of length (num_crystals + num_GB elements, but D_bulk_23 is only used for crystals (0 to num_crystals)
 *     "D_bulk_33": [...],    // Array of length (num_crystals + num_GB elements, but D_bulk_33 is only used for crystals (0 to num_crystals)
 *     "D_par":     [...],    // Array of length (num_crystals + num_GB elements), but D_par  is only used for GBs (num_crystals to num_crystals + num_GB)
 *     "D_perp":    [...]     // Array of length (num_crystals + num_GB elements), but D_perp is only used for GBs (num_crystals to num_crystals + num_GB)
 *   }
 */
class GBDiffusion : public ThermalModel, public LinearModel<1, 3> {
  public:
    GBDiffusion(const Reader &reader)
        : ThermalModel(reader)
    {
        vector<string> D_keys = {
            "D_bulk_11", "D_bulk_12", "D_bulk_13",
            "D_bulk_22", "D_bulk_23",
            "D_bulk_33"};

        try {
            // Read num_crystals, num_GB from the microstructure dataset attributes
            H5::H5File  file(reader.ms_filename, H5F_ACC_RDONLY);
            H5::DataSet ds = file.openDataSet(reader.ms_datasetname);
            ds.openAttribute("num_crystals").read(H5::PredType::NATIVE_INT64, &num_crystals);
            ds.openAttribute("num_GB").read(H5::PredType::NATIVE_INT64, &num_GB);

            n_mat       = num_crystals + num_GB;

            // Read rotation matrices for all grains (assumes dataset name 'rotation_matrices')
            grain_rot_matrices = FANS_malloc<double>(n_mat * 9);
            std::string rm_datasetname(reader.ms_datasetname);
            rm_datasetname.replace(rm_datasetname.find_last_of('/') + 1, std::string::npos, "rotation_matrices");
            file.openDataSet(rm_datasetname).read(grain_rot_matrices, H5::PredType::NATIVE_DOUBLE);

            // Read normal information for all GB interfaces (assumes dataset name 'normals)
            // NOTE: When indexing, shift index by num_crystals (starting from zero)
            GBnormals   = FANS_malloc<double>(num_GB * 3);
            std::string n_datasetname(reader.ms_datasetname);
            n_datasetname.replace(rm_datasetname.find_last_of('/') + 1, std::string::npos, "normals");
            file.openDataSet(n_datasetname).read(GBnormals, H5::PredType::NATIVE_DOUBLE);

            material_uniformity = reader.materialProperties["material_uniformity"].get<bool>();

            size_t num_constants = D_keys.size();
            D_bulk_constants.resize(num_constants, n_mat);

            D_par.resize(n_mat, 0.0);
            D_perp.resize(n_mat, 0.0);

            if (material_uniformity) {
                // Apply value to all bulk/GB phases
                D_bulk_constants.setZero();
                for (size_t k = 0; k < num_constants; ++k){
                    const double value = reader.materialProperties.at(D_keys[k]).get<double>();
                    fill_n(D_bulk_constants.row(k).begin(), num_crystals, value);
                }

                const double par_val  = reader.materialProperties["D_par"].get<double>();
                const double perp_val = reader.materialProperties["D_perp"].get<double>();

                fill_n(D_par.begin() + num_crystals, num_GB, par_val);
                fill_n(D_perp.begin() + num_crystals, num_GB, perp_val);
            } else {
                // Use phase-specific values
                for (int i = 0; i < n_mat; ++i) {
                    for (size_t k = 0; k < num_constants; ++k){
                        const auto &values = reader.materialProperties.at(D_keys[k]).get<vector<double>>();
                        if (values.size() != n_mat) {
                            throw std::runtime_error("Inconsistent size for material property: " + D_keys[k]);
                        }
                        D_bulk_constants.row(k) = Eigen::Map<const RowVectorXd>(values.data(), values.size());
                    }
                    D_par[i]  = reader.materialProperties["D_par"][i].get<double>();
                    D_perp[i] = reader.materialProperties["D_perp"][i].get<double>();
                }
            }

        } catch (const std::exception &e) {
            throw std::runtime_error("Error in GBDiffusion initialization: " + std::string(e.what()));
        }

        D_bulk_global.resize(num_crystals);
        kappa_average = Matrix3d::Zero();
        Matrix3d phase_kappa;
        phase_stiffness = new Matrix<double, 8, 8>[n_mat];

        for (size_t i = 0; i < n_mat; ++i) {
            phase_stiffness[i] = Matrix<double, 8, 8>::Zero();
            if (i < num_crystals) {
                // Reference grain diffusivity (crystal system)
                D_grain_ref << D_bulk_constants(0, i), D_bulk_constants(1, i), D_bulk_constants(2, i),
                    D_bulk_constants(1, i), D_bulk_constants(3, i), D_bulk_constants(4, i),
                    D_bulk_constants(2, i), D_bulk_constants(4, i), D_bulk_constants(5, i);

                // Rotation matrix
                rot_mat << grain_rot_matrices[9 * i + 0], grain_rot_matrices[9 * i + 1], grain_rot_matrices[9 * i + 2],
                        grain_rot_matrices[9 * i + 3], grain_rot_matrices[9 * i + 4], grain_rot_matrices[9 * i + 5],
                        grain_rot_matrices[9 * i + 6], grain_rot_matrices[9 * i + 7], grain_rot_matrices[9 * i + 8];

                // Store resulting bulk diffusivity tensor in global reference system
                D_bulk_global[i] = rot_mat * D_grain_ref * rot_mat.transpose();
                phase_kappa = D_bulk_global[i];
            } else if (i < n_mat) {
                // Grain boundary is transversely isotropic
                N           = Vector3d(GBnormals[3 * (i-num_crystals)+ 0], GBnormals[3 * (i-num_crystals) + 1], GBnormals[3 * (i-num_crystals) + 2]);
                N           = N.normalized();
                phase_kappa = D_par[i] * (Matrix3d::Identity() - N * N.transpose()) + D_perp[i] * N * N.transpose();
            } else {
                throw std::runtime_error("GBDiffusion: Unknown material index");
            }
            kappa_average += phase_kappa;
            for (int p = 0; p < n_gp; ++p) {
                phase_stiffness[i] += B_int[p].transpose() * phase_kappa * B_int[p] * v_e / n_gp;
            }
        }
        kappa_average = kappa_average / n_mat;
    }
    ~GBDiffusion() override
    {
        FANS_free(grain_rot_matrices);
        FANS_free(GBnormals);
        delete[] phase_stiffness;
        phase_stiffness = nullptr;
    }

    Matrix3d get_reference_stiffness() override
    {
        return kappa_average;
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        if (mat_index < num_crystals) {
            // Apply stored bulk diffusivities in global reference system
            sigma.block<3, 1>(i, 0) = D_bulk_global[mat_index] * eps.block<3, 1>(i, 0);

        } else if (mat_index < n_mat) {
            const ptrdiff_t base_idx = 3 * (mat_index - num_crystals);
            double          nx       = GBnormals[base_idx];
            double          ny       = GBnormals[base_idx + 1];
            double          nz       = GBnormals[base_idx + 2];

            // Pre-compute products for the projector matrix (N⊗N)
            double nxnx = nx * nx;
            double nxny = nx * ny;
            double nxnz = nx * nz;
            double nyny = ny * ny;
            double nynz = ny * nz;
            double nznz = nz * nz;

            // Pre-compute coefficients
            double d_diff = D_par[mat_index] - D_perp[mat_index];

            // Cache epsilon values to avoid repeated memory access
            double ex = eps(i, 0);
            double ey = eps(i + 1, 0);
            double ez = eps(i + 2, 0);

            // Calculate directly without constructing full matrices
            sigma(i, 0)     = D_par[mat_index] * ex - d_diff * (nxnx * ex + nxny * ey + nxnz * ez);
            sigma(i + 1, 0) = D_par[mat_index] * ey - d_diff * (nxny * ex + nyny * ey + nynz * ez);
            sigma(i + 2, 0) = D_par[mat_index] * ez - d_diff * (nxnz * ex + nynz * ey + nznz * ez);
        } else {
            throw std::runtime_error("GBDiffusion: Unknown material index");
        }
    }

    void postprocess(Solver<1, 3> &solver, Reader &reader, int load_idx, int time_idx) override
    {
        // Write first grain orientation axis to HDF5 file if requested
        if (find(reader.resultsToWrite.begin(), reader.resultsToWrite.end(), "grain_orientations_1") != reader.resultsToWrite.end()) {
            double *orientation_field = FANS_malloc<double>(solver.local_n0 * solver.n_y * solver.n_z * 3);

            for (ptrdiff_t element_idx = 0; element_idx < solver.local_n0 * solver.n_y * solver.n_z; ++element_idx) {
                int mat_index = solver.ms[element_idx];
                if (mat_index < num_crystals) {
                    orientation_field[element_idx * 3]         = grain_rot_matrices[9 * mat_index];
                    orientation_field[element_idx * 3 + 1]     = grain_rot_matrices[9 * mat_index + 1];
                    orientation_field[element_idx * 3 + 2]     = grain_rot_matrices[9 * mat_index + 2];
                }
            }
            reader.writeSlab("grain_orientations_1", load_idx, time_idx, orientation_field, {3});
            FANS_free(orientation_field);
        }

        // Write second grain orientation axis to HDF5 file if requested
        if (find(reader.resultsToWrite.begin(), reader.resultsToWrite.end(), "grain_orientations_2") != reader.resultsToWrite.end()) {
            double *orientation_field = FANS_malloc<double>(solver.local_n0 * solver.n_y * solver.n_z * 3);

            for (ptrdiff_t element_idx = 0; element_idx < solver.local_n0 * solver.n_y * solver.n_z; ++element_idx) {
                int mat_index = solver.ms[element_idx];
                if (mat_index < num_crystals) {
                    orientation_field[element_idx * 3]         = grain_rot_matrices[9 * mat_index + 3];
                    orientation_field[element_idx * 3 + 1]     = grain_rot_matrices[9 * mat_index + 4];
                    orientation_field[element_idx * 3 + 2]     = grain_rot_matrices[9 * mat_index + 5];
                }
            }
            reader.writeSlab("grain_orientations_2", load_idx, time_idx, orientation_field, {3});
            FANS_free(orientation_field);
        }

        // Write third grain orientation axis to HDF5 file if requested
        if (find(reader.resultsToWrite.begin(), reader.resultsToWrite.end(), "grain_orientations_3") != reader.resultsToWrite.end()) {
            double *orientation_field = FANS_malloc<double>(solver.local_n0 * solver.n_y * solver.n_z * 3);

            for (ptrdiff_t element_idx = 0; element_idx < solver.local_n0 * solver.n_y * solver.n_z; ++element_idx) {
                int mat_index = solver.ms[element_idx];
                if (mat_index < num_crystals) {
                    orientation_field[element_idx * 3]         = grain_rot_matrices[9 * mat_index + 6];
                    orientation_field[element_idx * 3 + 1]     = grain_rot_matrices[9 * mat_index + 7];
                    orientation_field[element_idx * 3 + 2]     = grain_rot_matrices[9 * mat_index + 8];
                }
            }
            reader.writeSlab("grain_orientations_3", load_idx, time_idx, orientation_field, {3});
            FANS_free(orientation_field);
        }

        // Write GBnormals to HDF5 file if requested
        if (find(reader.resultsToWrite.begin(), reader.resultsToWrite.end(), "GBnormals") != reader.resultsToWrite.end()) {
            double *GBnormals_field = FANS_malloc<double>(solver.local_n0 * solver.n_y * solver.n_z * 3);
            for (ptrdiff_t element_idx = 0; element_idx < solver.local_n0 * solver.n_y * solver.n_z; ++element_idx) {
                int mat_index_GB = solver.ms[element_idx] - num_crystals;
                if (mat_index_GB >= 0) {
                    GBnormals_field[element_idx * 3]     = GBnormals[3 * mat_index_GB];
                    GBnormals_field[element_idx * 3 + 1] = GBnormals[3 * mat_index_GB + 1];
                    GBnormals_field[element_idx * 3 + 2] = GBnormals[3 * mat_index_GB + 2];
                }
            }
            reader.writeSlab("GBnormals", load_idx, time_idx, GBnormals_field, {3});
            FANS_free(GBnormals_field);
        }
    }

  private:
    int  num_crystals = 0;
    int  num_GB       = 0;
    bool material_uniformity;

    MatrixXd D_bulk_constants;
    std::vector<Matrix3d, Eigen::aligned_allocator<Matrix3d>> D_bulk_global;
    vector<double> D_par;
    vector<double> D_perp;

    double  *grain_rot_matrices = nullptr;
    Matrix3d rot_mat;
    Matrix3d D_grain_ref;
    double  *GBnormals = nullptr;
    Vector3d N;
    Matrix3d kappa_average;
};

#endif // GBDIFFUSION_H
