/**
 * @file PseudoPlastic.h
 * @brief This file contains the declaration of the PseudoPlastic class and its derived classes.
 *
 * The PseudoPlastic class is a base class that represents a pseudo-plastic material model.
 * The models implemented are described in https://doi.org/10.1016/j.euromechsol.2017.11.007 -> Appendix A.1 and A.2.
 * It contains common properties and methods for all pseudo-plastic material models. The derived classes,
 * PseudoPlasticLinearHardening and PseudoPlasticNonLinearHardening, implement specific variations
 * of the pseudo-plastic material model with linear and nonlinear hardening, respectively.
 *
 * The PseudoPlastic class and its derived classes are used in the FANS for simulating
 * mechanical behavior of materials.
 */

#ifndef PSEUDOPLASTIC_H
#define PSEUDOPLASTIC_H

#include "matmodel.h"
#include "solver.h"

class PseudoPlastic : public SmallStrainMechModel {
  public:
    PseudoPlastic(const Reader &reader)
        : SmallStrainMechModel(reader)
    {
        try {
            bulk_modulus  = reader.materialProperties["bulk_modulus"].get<vector<double>>();
            shear_modulus = reader.materialProperties["shear_modulus"].get<vector<double>>();
            yield_stress  = reader.materialProperties["yield_stress"].get<vector<double>>();
        } catch (const std::exception &e) {
            throw std::runtime_error("Missing material properties for the requested material model.");
        }
        n_mat = bulk_modulus.size();

        if (kapparef_mat.isZero()) {
            const double Kbar      = std::accumulate(bulk_modulus.begin(), bulk_modulus.end(), 0.0) / static_cast<double>(n_mat);
            const double Gbar      = std::accumulate(shear_modulus.begin(), shear_modulus.end(), 0.0) / static_cast<double>(n_mat);
            const double lambdabar = Kbar - 2.0 * Gbar / 3.0;
            kapparef_mat.setZero();
            kapparef_mat.topLeftCorner(3, 3).setConstant(lambdabar);
            kapparef_mat.diagonal().array() += 2.0 * Gbar;
        }
    }

    void initializeInternalVariables(ptrdiff_t num_elements, int num_gauss_points) override
    {
        plastic_flag.resize(num_elements, VectorXi::Zero(num_gauss_points));
    }

    virtual void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override = 0; // Pure virtual method

    void postprocess(Solver<3, 6> &solver, Reader &reader, int load_idx, int time_idx) override
    {
        VectorXf element_plastic_flag = VectorXf::Zero(solver.local_n0 * solver.n_y * solver.n_z);
        for (ptrdiff_t elem_idx = 0; elem_idx < solver.local_n0 * solver.n_y * solver.n_z; ++elem_idx) {
            element_plastic_flag(elem_idx) = plastic_flag[elem_idx].cast<float>().mean();
        }

        reader.writeSlab("plastic_flag", load_idx, time_idx, element_plastic_flag.data(), 1);
    }

  protected:
    vector<double>       bulk_modulus;
    vector<double>       shear_modulus;
    vector<double>       yield_stress;
    vector<double>       eps_crit;
    vector<VectorXi>     plastic_flag;
    Matrix<double, 6, 1> dev_eps;
    double               treps, norm_dev_eps, buf1, buf2;
};

class PseudoPlasticLinearHardening : public PseudoPlastic {
  public:
    PseudoPlasticLinearHardening(const Reader &reader)
        : PseudoPlastic(reader)
    {
        try {
            hardening_parameter = reader.materialProperties["hardening_parameter"].get<vector<double>>();
        } catch (const std::exception &e) {
            throw std::runtime_error("Missing material properties for the requested material model.");
        }

        E_s.resize(n_mat);
        eps_crit.resize(n_mat);

        for (int i = 0; i < n_mat; ++i) {
            eps_crit[i] = sqrt(2. / 3.) * yield_stress[i] / (2. * shear_modulus[i]);
            E_s[i]      = (3. * shear_modulus[i]) / (3. * shear_modulus[i] + hardening_parameter[i]);
        }
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        treps             = eps.block<3, 1>(i, 0).sum();
        dev_eps.head<3>() = eps.block<3, 1>(i, 0) - (1.0 / 3.0) * treps * Matrix<double, 3, 1>::Ones();
        dev_eps.tail<3>() = eps.block<3, 1>(i + 3, 0);

        norm_dev_eps = dev_eps.norm();
        buf1         = bulk_modulus[mat_index] * treps;

        if (norm_dev_eps <= eps_crit[mat_index]) {
            buf2                             = 2.0 * shear_modulus[mat_index];
            plastic_flag[element_idx](i / 6) = mat_index;
        } else {
            buf2 = (b * yield_stress[mat_index] + a * E_s[mat_index] * hardening_parameter[mat_index] *
                                                      (norm_dev_eps - eps_crit[mat_index])) /
                   norm_dev_eps;
            plastic_flag[element_idx](i / 6) = this->n_mat + mat_index;
        }
        sigma.block<3, 1>(i, 0).setConstant(buf1);
        sigma.block<3, 1>(i, 0) += buf2 * dev_eps.head<3>();
        sigma.block<3, 1>(i + 3, 0) = buf2 * dev_eps.tail<3>();
    }

  private:
    vector<double> hardening_parameter;
    vector<double> E_s;
    double         a = 2. / 3;
    double         b = sqrt(a);
};

class PseudoPlasticNonLinearHardening : public PseudoPlastic {
  public:
    PseudoPlasticNonLinearHardening(const Reader &reader)
        : PseudoPlastic(reader)
    {
        try {
            hardening_exponent = reader.materialProperties["hardening_exponent"].get<vector<double>>();
            eps_0              = reader.materialProperties["eps_0"].get<vector<double>>(); // ε0 parameter
        } catch (const std::exception &e) {
            throw std::runtime_error("Missing material properties for the requested material model.");
        }

        eps_crit.resize(n_mat);
        for (int i = 0; i < n_mat; ++i) {
            eps_crit[i] = eps_0[i] * pow(yield_stress[i] / (3.0 * shear_modulus[i] * eps_0[i]), 1.0 / (1.0 - hardening_exponent[i]));
        }
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        treps             = eps.block<3, 1>(i, 0).sum();
        dev_eps.head<3>() = eps.block<3, 1>(i, 0) - (1.0 / 3.0) * treps * Matrix<double, 3, 1>::Ones();
        dev_eps.tail<3>() = eps.block<3, 1>(i + 3, 0);

        norm_dev_eps = sqrt(2.0 / 3.0) * dev_eps.norm(); // ε_eq

        buf1 = bulk_modulus[mat_index] * treps;
        sigma.block<6, 1>(i, 0).setConstant(0);
        if (norm_dev_eps <= eps_crit[mat_index]) {
            buf2 = 2.0 * shear_modulus[mat_index];
            sigma.block<3, 1>(i, 0).setConstant(buf1);
            sigma.block<3, 1>(i, 0) += buf2 * dev_eps.head<3>();
            sigma.block<3, 1>(i + 3, 0) = buf2 * dev_eps.tail<3>();

            plastic_flag[element_idx](i / 6) = mat_index;
        } else {
            buf2 = sqrt(2.0 / 3.0) * yield_stress[mat_index] *
                   pow(norm_dev_eps / eps_0[mat_index], hardening_exponent[mat_index]);
            sigma.block<3, 1>(i, 0).setConstant(buf1);
            sigma.block<6, 1>(i, 0) += buf2 * dev_eps / dev_eps.norm();

            plastic_flag[element_idx](i / 6) = this->n_mat + mat_index;
        }
    }

  private:
    vector<double> hardening_exponent;
    vector<double> eps_0;
};

#endif // PSEUDOPLASTIC_H
