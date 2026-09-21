#ifndef GBDIFFUSION_H
#define GBDIFFUSION_H

#include "matmodel.h"
#include <array>
#include <cstdint>
#include <Eigen/StdVector> // For Eigen's aligned_allocator

/**
 * @class GBDiffusion
 * @brief Material model for grain boundary diffusion in polycrystals
 *
 * This model implements diffusion in a polycrystalline material, differentiating between
 * bulk crystal diffusion and grain boundary diffusion (transversely isotropic).
 * The grains are characterized by an active crystal-to-sample rotation matrix and can
 * have different diffusion properties along the three crystal axes.
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
 *
 * Required material parameters in JSON format:
 *   - material_uniformity: Boolean flag for uniform grain and GB properties
 *
 *   When material_uniformity is true (uniform properties):
 *   {
 *     "material_uniformity": true,
 *     "D_bulk_11": 10.0,    // Diffusion coefficient (1,1) for all crystals in crystal system
 *     "D_bulk_12": 0.0,     // Diffusion coefficient (1,2) for all crystals in crystal system
 *     "D_bulk_13": 0.0,     // Diffusion coefficient (1,3) for all crystals in crystal system
 *     "D_bulk_22": 5.0,     // Diffusion coefficient (2,2) for all crystals in crystal system
 *     "D_bulk_23": 0.0,     // Diffusion coefficient (2,3) for all crystals in crystal system
 *     "D_bulk_33": 1.0,     // Diffusion coefficient (3,3) for all crystals in crystal system
 *     "D_par": 2.0,         // Diffusion coefficient parallel to the grain boundary for all GBs
 *     "D_perp": 0.5         // Diffusion coefficient perpendicular to the grain boundary for all GBs
 *   }
 *
 *   When material_uniformity is false (tag-specific properties):
 *   {
 *     "material_uniformity": false,
 *     "D_bulk_11": [...],    // One value per phase; only crystal entries are used
 *     "D_bulk_12": [...],
 *     "D_bulk_13": [...],
 *     "D_bulk_22": [...],
 *     "D_bulk_23": [...],
 *     "D_bulk_33": [...],
 *     "D_par":     [...],    // One value per phase; only GB entries are used
 *     "D_perp":    [...]
 *   }
 */
class GBDiffusion : public ThermalModel, public LinearModel<1, 3> {
  public:
    GBDiffusion(const Reader &reader)
        : ThermalModel(reader)
    {
        static constexpr std::array<const char *, 6> D_keys = {
            "D_bulk_11", "D_bulk_12", "D_bulk_13",
            "D_bulk_22", "D_bulk_23",
            "D_bulk_33"};

        try {
            H5::H5File   file(reader.ms_filename, H5F_ACC_RDONLY);
            H5::DataSet  ds = file.openDataSet(reader.ms_datasetname);
            std::int64_t crystal_count, boundary_count;
            ds.openAttribute("num_crystals").read(H5::PredType::NATIVE_INT64, &crystal_count);
            ds.openAttribute("num_GB").read(H5::PredType::NATIVE_INT64, &boundary_count);
            const int num_crystals = static_cast<int>(crystal_count);
            const int num_GB       = static_cast<int>(boundary_count);
            n_mat                  = num_crystals + num_GB;

            auto sibling = [&](const char *name) {
                std::string path(reader.ms_datasetname);
                path.replace(path.find_last_of('/') + 1, std::string::npos, name);
                return path;
            };
            vector<double> grain_rot_matrices(9 * num_crystals), GB_normals(3 * n_mat);
            file.openDataSet(sibling("rotation_matrices")).read(grain_rot_matrices.data(), H5::PredType::NATIVE_DOUBLE);
            file.openDataSet(sibling("GB_normals")).read(GB_normals.data(), H5::PredType::NATIVE_DOUBLE);

            Matrix<double, 6, Dynamic> D_bulk_constants = Matrix<double, 6, Dynamic>::Zero(6, n_mat);
            VectorXd                   D_par = VectorXd::Zero(n_mat), D_perp = VectorXd::Zero(n_mat);
            if (reader.materialProperties["material_uniformity"].get<bool>()) {
                for (size_t k = 0; k < D_keys.size(); ++k)
                    D_bulk_constants.row(k).head(num_crystals).setConstant(reader.materialProperties.at(D_keys[k]).get<double>());
                D_par.tail(num_GB).setConstant(reader.materialProperties["D_par"].get<double>());
                D_perp.tail(num_GB).setConstant(reader.materialProperties["D_perp"].get<double>());
            } else {
                auto values = [&](const char *name) {
                    auto result = reader.materialProperties.at(name).get<vector<double>>();
                    if (result.size() != static_cast<size_t>(n_mat))
                        throw std::runtime_error("Inconsistent size for material property: " + string(name));
                    return result;
                };
                for (size_t k = 0; k < D_keys.size(); ++k) {
                    const auto data         = values(D_keys[k]);
                    D_bulk_constants.row(k) = Map<const RowVectorXd>(data.data(), n_mat);
                }
                const auto par = values("D_par"), perp = values("D_perp");
                D_par  = Map<const VectorXd>(par.data(), n_mat);
                D_perp = Map<const VectorXd>(perp.data(), n_mat);
            }

            phase_diffusivities.resize(n_mat);
            phase_stiffness_storage.resize(n_mat);
            phase_stiffness = phase_stiffness_storage.data();
            kappa_average.setZero();

            for (int phase = 0; phase < n_mat; ++phase) {
                auto &phase_kappa = phase_diffusivities[phase];
                if (phase < num_crystals) {
                    Matrix3d D_grain_ref;
                    D_grain_ref << D_bulk_constants(0, phase), D_bulk_constants(1, phase), D_bulk_constants(2, phase),
                        D_bulk_constants(1, phase), D_bulk_constants(3, phase), D_bulk_constants(4, phase),
                        D_bulk_constants(2, phase), D_bulk_constants(4, phase), D_bulk_constants(5, phase);
                    using RowMatrix3d = Matrix<double, 3, 3, RowMajor>;
                    const Map<const RowMatrix3d> rot_mat(grain_rot_matrices.data() + 9 * phase);
                    phase_kappa.noalias() = rot_mat * D_grain_ref * rot_mat.transpose();
                } else {
                    const Vector3d normal = Map<const Vector3d>(GB_normals.data() + 3 * phase).normalized();
                    phase_kappa           = D_par[phase] * (Matrix3d::Identity() - normal * normal.transpose()) + D_perp[phase] * normal * normal.transpose();
                }

                kappa_average += phase_kappa;
                phase_stiffness[phase].setZero();
                for (const auto &B : B_int)
                    phase_stiffness[phase].noalias() += B.transpose() * phase_kappa * B * v_e / n_gp;
            }
            kappa_average /= n_mat;
        } catch (const std::exception &e) {
            throw std::runtime_error("Error in GBDiffusion initialization: " + std::string(e.what()));
        }
    }

    Matrix3d get_reference_stiffness() override
    {
        return kappa_average;
    }

    void get_sigma(int i, int mat_index, ptrdiff_t) override
    {
        sigma.segment<3>(i).noalias() = phase_diffusivities[mat_index] * eps.segment<3>(i);
    }

  private:
    std::vector<Matrix3d, Eigen::aligned_allocator<Matrix3d>>                         phase_diffusivities;
    std::vector<Matrix<double, 8, 8>, Eigen::aligned_allocator<Matrix<double, 8, 8>>> phase_stiffness_storage;
    Matrix3d                                                                          kappa_average;
};

#endif // GBDIFFUSION_H
