#ifndef GBDIFFUSION_H
#define GBDIFFUSION_H

#include "matmodel.h"
#include <Eigen/StdVector> // For Eigen's aligned_allocator

/**
 * @class GBDiffusion
 * @brief Material model for grain boundary diffusion in polycrystals
 *
 * This model implements diffusion in a polycrystalline material, differentiating between
 * bulk crystal diffusion (isotropic) and grain boundary diffusion (transversely isotropic).
 * The grain boundaries are characterized by their normal vectors and can have different
 * diffusion properties parallel and perpendicular to the boundary plane.
 *
 * The model extends both ThermalModel and LinearModel<1, 3> to provide a linear diffusion
 * formulation that can be used in a thermal-like solver in FANS.
 *
 * @details The model:
 *   - Reads microstructure data containing grain boundaries from HDF5 files
 *   - Supports uniform or material-specific diffusivity values
 *   - Handles bulk regions with isotropic diffusion (D_bulk)
 *   - Handles grain boundaries with transversely isotropic diffusion (D_par, D_perp)
 *   - Provides visualization of grain boundary normals in post-processing
 *
 * Required material parameters in JSON format:
 *   - GB_unformity: Boolean flag for uniform GB properties
 *
 *   When GB_unformity is true (uniform properties):
 *   {
 *     "GB_unformity": true,
 *     "D_bulk": 1.0,        // Isotropic diffusion coefficient for all crystals
 *     "D_par": 2.0,         // Diffusion coefficient parallel to the grain boundary for all GBs
 *     "D_perp": 0.5         // Diffusion coefficient perpendicular to the grain boundary for all GBs
 *   }
 *
 *   When GB_unformity is false (tag-specific properties):
 *   {
 *     "GB_unformity": false,
 *     "D_bulk": [...],  // Array of length (num_crystals + num_GB elements), but D_bulk is only used for crystals (0 to num_crystals)
 *     "D_par":  [...],  // Array of length (num_crystals + num_GB elements), but D_par  is only used for GBs (num_crystals to num_crystals + num_GB)
 *     "D_perp": [...]   // Array of length (num_crystals + num_GB elements), but D_perp is only used for GBs (num_crystals to num_crystals + num_GB)
 *   }
 */
class GBDiffusion : public ThermalModel, public LinearModel<1, 3> {
  public:
    GBDiffusion(const Reader &reader)
        : ThermalModel(reader)
    {
        try {
            // Read num_crystals, num_GB and GBVoxelInfo from the microstructure dataset attributes
            H5::H5File  file(reader.ms_filename, H5F_ACC_RDONLY);
            H5::DataSet ds = file.openDataSet(reader.ms_datasetname);
            ds.openAttribute("num_crystals").read(H5::PredType::NATIVE_INT64, &num_crystals);
            ds.openAttribute("num_GB").read(H5::PredType::NATIVE_INT64, &num_GB);
            std::string   json_text;
            H5::Attribute attr    = ds.openAttribute("GBVoxelInfo");
            H5::StrType   strType = attr.getStrType();
            attr.read(strType, json_text);

            n_mat       = num_crystals + num_GB;
            GBnormals   = FANS_malloc<double>(n_mat * 3);
            auto gbInfo = json::parse(json_text);
            for (auto &kv : gbInfo.items()) {
                int   tag    = kv.value().at("GB_tag").get<int>();
                auto &normal = kv.value()["GB_normal"];

                GBnormals[(tag) * 3]     = normal[0].get<double>();
                GBnormals[(tag) * 3 + 1] = normal[1].get<double>();
                GBnormals[(tag) * 3 + 2] = normal[2].get<double>();
            }
            GB_unformity = reader.materialProperties["GB_unformity"].get<bool>();

            D_bulk.resize(n_mat, 0.0);
            D_par.resize(n_mat, 0.0);
            D_perp.resize(n_mat, 0.0);

            if (GB_unformity) {
                double bulk_val = reader.materialProperties["D_bulk"].get<double>();
                double par_val  = reader.materialProperties["D_par"].get<double>();
                double perp_val = reader.materialProperties["D_perp"].get<double>();

                fill_n(D_bulk.begin(), num_crystals, bulk_val);
                fill_n(D_par.begin() + num_crystals, num_GB, par_val);
                fill_n(D_perp.begin() + num_crystals, num_GB, perp_val);
            } else {
                for (int i = 0; i < n_mat; ++i) {
                    D_bulk[i] = reader.materialProperties["D_bulk"][i].get<double>();
                    D_par[i]  = reader.materialProperties["D_par"][i].get<double>();
                    D_perp[i] = reader.materialProperties["D_perp"][i].get<double>();
                }
            }

        } catch (const std::exception &e) {
            throw std::runtime_error("Error in GBDiffusion initialization: " + std::string(e.what()));
        }

        kappa_average = Matrix3d::Zero();
        Matrix3d phase_kappa;
        phase_stiffness = new Matrix<double, 8, 8>[n_mat];

        for (size_t i = 0; i < n_mat; ++i) {
            phase_stiffness[i] = Matrix<double, 8, 8>::Zero();
            if (i < num_crystals) {
                // Bulk is isotropic
                phase_kappa = D_bulk[i] * Matrix3d::Identity();
            } else if (i < n_mat) {
                // Grain boundary is transversely isotropic
                N           = Vector3d(GBnormals[3 * i + 0], GBnormals[3 * i + 1], GBnormals[3 * i + 2]);
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
    ~GBDiffusion()
    {
        FANS_free(GBnormals);
        delete[] phase_stiffness;
        phase_stiffness = nullptr;
    }

    Matrix3d get_reference_stiffness() const override
    {
        return kappa_average;
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        if (mat_index < num_crystals) {
            sigma.block<3, 1>(i, 0) = D_bulk[mat_index] * eps.block<3, 1>(i, 0);
        } else if (mat_index < n_mat) {
            const ptrdiff_t base_idx = 3 * mat_index;
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
        // Write GBnormals to HDF5 file if requested
        if (find(reader.resultsToWrite.begin(), reader.resultsToWrite.end(), "GBnormals") != reader.resultsToWrite.end()) {
            double *GBnormals_field = FANS_malloc<double>(solver.local_n0 * solver.n_y * solver.n_z * 3);
            for (ptrdiff_t element_idx = 0; element_idx < solver.local_n0 * solver.n_y * solver.n_z; ++element_idx) {
                int mat_index = solver.ms[element_idx];
                if (mat_index >= num_crystals) {
                    GBnormals_field[element_idx * 3]     = GBnormals[3 * mat_index];
                    GBnormals_field[element_idx * 3 + 1] = GBnormals[3 * mat_index + 1];
                    GBnormals_field[element_idx * 3 + 2] = GBnormals[3 * mat_index + 2];
                }
            }
            reader.writeSlab("GBnormals", load_idx, time_idx, GBnormals_field, 3);
            FANS_free(GBnormals_field);
        }
    }

  private:
    int  num_crystals = 0;
    int  num_GB       = 0;
    bool GB_unformity;

    vector<double> D_bulk;
    vector<double> D_par;
    vector<double> D_perp;

    double  *GBnormals = nullptr;
    Vector3d N;
    Matrix3d kappa_average;
};

#endif // GBDIFFUSION_H
