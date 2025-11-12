#ifndef MATERIALMANAGER_H
#define MATERIALMANAGER_H

#include "general.h"
#include "matmodel.h"

template <int howmany, int n_str>
Matmodel<howmany, n_str> *createMatmodel(const Reader &reader);

/**
 * @brief Lightweight struct for fast material model lookup in hot loops
 */
template <int howmany, int n_str>
struct MaterialInfo {
    Matmodel<howmany, n_str>    *model;
    LinearModel<howmany, n_str> *linear_model; // Pre-cast pointer (nullptr if nonlinear)
    int                          local_mat_id;
    bool                         is_linear;
};

/**
 * @brief MaterialManager: Manages multiple material models for multi-material simulations
 *
 * Enables different phases in the microstructure to use different constitutive models.
 *
 */
template <int howmany, int n_str>
class MaterialManager {
  private:
    MaterialInfo<howmany, n_str> *phase_to_info{nullptr}; // [n_phases] - HOT DATA
    int                           n_phases;

  public:
    vector<Matmodel<howmany, n_str> *> models;           // vector of unique material models
    Matrix<double, n_str, n_str>       kapparef_mat;     // Reference stiffness for fundamental solution
    bool                               all_linear{true}; // True if ALL phases use linear models

    MaterialManager(const Reader &reader)
    {
        const auto &mats = reader.inputJson.at("materials");
        if (!mats.is_array() || mats.empty())
            throw std::runtime_error("MaterialManager: 'materials' must be non-empty array");

        int max_phase = -1;
        for (const auto &mg : mats) {
            if (!mg.contains("phases") || !mg.contains("matmodel") || !mg.contains("material_properties"))
                throw std::runtime_error("MaterialManager: material group missing required fields");

            for (int p : mg["phases"].get<vector<int>>())
                max_phase = std::max(max_phase, p);
        }

        n_phases = max_phase + 1;
        if (n_phases == 0)
            throw std::runtime_error("MaterialManager: No phases defined");

        phase_to_info = new MaterialInfo<howmany, n_str>[n_phases]();

        // Create models and map phases
        for (const auto &mg : mats) {
            auto *model = create_material_model_from_json(mg, reader);
            models.push_back(model);

            auto *linear_model = dynamic_cast<LinearModel<howmany, n_str> *>(model);
            bool  is_linear    = (linear_model != nullptr);
            if (!is_linear)
                all_linear = false;

            auto phases = mg["phases"].get<vector<int>>();
            for (size_t i = 0; i < phases.size(); ++i) {
                int p = phases[i];
                if (p < 0 || p >= n_phases || phase_to_info[p].model)
                    throw std::runtime_error("MaterialManager: Invalid or duplicate phase " + std::to_string(p));
                phase_to_info[p] = {model, linear_model, static_cast<int>(i), is_linear};
            }
        }

        // Verify all phases assigned
        for (int p = 0; p < n_phases; ++p)
            if (!phase_to_info[p].model)
                throw std::runtime_error("MaterialManager: Phase " + std::to_string(p) + " not assigned");

        compute_reference_stiffness(reader);

        // Print detailed information about material configuration for logging
        if (reader.world_rank == 0) {
            printf("\n# MaterialManager initialized:\n");
            printf("#   Number of material models: %zu\n", models.size());
            printf("#   Number of phases: %d\n#\n", n_phases);

            for (size_t i = 0; i < mats.size(); ++i) {
                const auto &mg = mats[i];
                printf("# Material Model %zu: %s\n", i + 1, mg["matmodel"].get<string>().c_str());

                // Print phases
                auto phases = mg["phases"].get<vector<int>>();
                printf("#   Phases: [");
                for (size_t j = 0; j < phases.size(); ++j) {
                    printf("%d", phases[j]);
                    if (j < phases.size() - 1)
                        printf(", ");
                }
                printf("]\n");

                // Print material properties
                printf("#   Material properties:\n");
                const auto &props = mg["material_properties"];
                for (auto it = props.begin(); it != props.end(); ++it) {
                    printf("#     %s: ", it.key().c_str());
                    if (it.value().is_array()) {
                        printf("[");
                        for (size_t k = 0; k < it.value().size(); ++k) {
                            if (it.value()[k].is_number()) {
                                printf("%.5g", it.value()[k].get<double>());
                            } else if (it.value()[k].is_string()) {
                                printf("\"%s\"", it.value()[k].get<string>().c_str());
                            }
                            if (k < it.value().size() - 1)
                                printf(", ");
                        }
                        printf("]");
                    } else if (it.value().is_number()) {
                        printf("%.5g", it.value().get<double>());
                    } else if (it.value().is_string()) {
                        printf("\"%s\"", it.value().get<string>().c_str());
                    }
                    printf("\n");
                }
                printf("#\n");
            }
        }
    }

    ~MaterialManager()
    {
        for (auto *m : models)
            delete m;
        delete[] phase_to_info;
    }

    // Initialize internal variables for all material models
    void initialize_internal_variables(ptrdiff_t num_elements, int num_gauss_points)
    {
        for (auto *model : models) {
            model->initializeInternalVariables(num_elements, num_gauss_points);
        }
    }

    // Compute reference stiffness for the fundamental solution
    void compute_reference_stiffness(const Reader &reader)
    {
        if (reader.inputJson.contains("reference_material")) {
            // Check for optional user-defined reference material
            auto ref_mat = reader.inputJson["reference_material"].get<vector<vector<double>>>();
            if (ref_mat.size() != n_str || ref_mat[0].size() != n_str) {
                throw std::invalid_argument("reference_material must be " + std::to_string(n_str) + "x" + std::to_string(n_str));
            }

            for (int i = 0; i < n_str; ++i) {
                kapparef_mat.row(i) = Eigen::Map<const Eigen::RowVectorXd>(ref_mat[i].data(), n_str);
            }
            Eigen::LLT<Matrix<double, n_str, n_str>> llt(kapparef_mat);
            if (llt.info() != Eigen::Success) {
                throw std::invalid_argument("reference_material must be symmetric positive definite");
            }

            if (reader.world_rank == 0) {
                cout << "# Using user-defined reference material for fundamental solution." << endl;
            }
        } else {
            // Default: simple average of reference stiffness across all material models
            kapparef_mat = Matrix<double, n_str, n_str>::Zero();
            for (auto *model : models) {
                kapparef_mat += model->get_reference_stiffness();
            }
            kapparef_mat /= static_cast<double>(models.size());
        }
    }

    // Update internal variables after converged time step
    void update_internal_variables()
    {
        for (auto *model : models) {
            model->updateInternalVariables();
        }
    }

    // Set macroscale loading gradient for all models
    void set_gradient(vector<double> g0)
    {
        for (auto *model : models) {
            model->setGradient(g0);
        }
    }

    // Post-processing for all material models
    void postprocess(Solver<howmany, n_str> &solver, Reader &reader, int load_idx, int time_idx)
    {
        for (auto *model : models) {
            model->postprocess(solver, reader, load_idx, time_idx);
        }
    }

    // Get material info for a phase
    inline const MaterialInfo<howmany, n_str> &get_info(int phase_id) const
    {
        return phase_to_info[phase_id];
    }
    inline size_t get_num_models() const
    {
        return models.size();
    } // Get number of unique material models
    inline int get_num_phases() const
    {
        return n_phases;
    } // Get number of phases in microstructure

  private:
    // Helper function to create a material model from a material group in JSON
    Matmodel<howmany, n_str> *create_material_model_from_json(
        const json   &mat_group,
        const Reader &base_reader)
    {
        // Create a minimal temporary reader
        Reader temp_reader;

        // Copy safe members
        temp_reader.world_rank  = base_reader.world_rank;
        temp_reader.world_size  = base_reader.world_size;
        temp_reader.FE_type     = base_reader.FE_type;
        temp_reader.strain_type = base_reader.strain_type;
        temp_reader.problemType = base_reader.problemType;
        temp_reader.method      = base_reader.method;
        temp_reader.l_e         = base_reader.l_e;
        temp_reader.dims        = base_reader.dims;

        // Override material properties and n_mat for this specific material group
        temp_reader.materialProperties = mat_group["material_properties"];

        // n_mat for this model is the number of phases using it
        auto phases       = mat_group["phases"].get<vector<int>>();
        temp_reader.n_mat = phases.size();

        // Override matmodel name for factory function
        temp_reader.matmodel = mat_group["matmodel"].get<string>();

        return createMatmodel<howmany, n_str>(temp_reader);
    }
};

#endif // MATERIALMANAGER_H
