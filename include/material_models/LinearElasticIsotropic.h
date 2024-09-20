#ifndef LINEARELASTICISOTROPIC_H
#define LINEARELASTICISOTROPIC_H

#include "matmodel.h"

class LinearElasticIsotropic : public MechModel, public LinearModel<3> {
  public:
    LinearElasticIsotropic(vector<double> l_e, map<string, vector<double>> materialProperties)
        : MechModel(l_e)
    {
        try {
            bulk_modulus = materialProperties["bulk_modulus"];
            mu           = materialProperties["shear_modulus"];
        } catch (const std::exception &e) {
            throw std::runtime_error("Missing material properties for the requested material model.");
        }

        n_mat = bulk_modulus.size();
        lambda.resize(n_mat);
        mu.resize(n_mat);

        for (int i = 0; i < n_mat; ++i) {
            lambda[i] = bulk_modulus[i] - (2.0 / 3.0) * mu[i];
        }

        double lambda_ref = (*max_element(lambda.begin(), lambda.end()) +
                             *min_element(lambda.begin(), lambda.end())) /
                            2;
        double mu_ref = (*max_element(mu.begin(), mu.end()) +
                         *min_element(mu.begin(), mu.end())) /
                        2;

        kapparef_mat = Matrix<double, 6, 6>::Zero();
        kapparef_mat.topLeftCorner(3, 3).setConstant(lambda_ref);
        kapparef_mat += 2 * mu_ref * Matrix<double, 6, 6>::Identity();

        phase_stiffness = new Matrix<double, 24, 24>[n_mat];
        Matrix<double, 6, 6> phase_kappa;

        for (int i = 0; i < n_mat; i++) {
            phase_kappa.setZero();
            phase_stiffness[i] = Matrix<double, 24, 24>::Zero();

            phase_kappa.topLeftCorner(3, 3).setConstant(lambda[i]);
            phase_kappa += 2 * mu[i] * Matrix<double, 6, 6>::Identity();

            for (int p = 0; p < 8; ++p) {
                phase_stiffness[i] += B_int[p].transpose() * phase_kappa * B_int[p] * v_e * 0.1250;
            }
        }
    }

    void get_sigma(int i, int mat_index, ptrdiff_t element_idx) override
    {
        double buf1     = lambda[mat_index] * (eps(i, 0) + eps(i + 1, 0) + eps(i + 2, 0));
        double buf2     = 2 * mu[mat_index];
        sigma(i + 0, 0) = buf1 + buf2 * eps(i + 0, 0);
        sigma(i + 1, 0) = buf1 + buf2 * eps(i + 1, 0);
        sigma(i + 2, 0) = buf1 + buf2 * eps(i + 2, 0);
        sigma(i + 3, 0) = buf2 * eps(i + 3, 0);
        sigma(i + 4, 0) = buf2 * eps(i + 4, 0);
        sigma(i + 5, 0) = buf2 * eps(i + 5, 0);
    }

  private:
    vector<double> bulk_modulus;
    vector<double> lambda;
    vector<double> mu;
};

#endif // LINEARELASTICISOTROPIC_H
