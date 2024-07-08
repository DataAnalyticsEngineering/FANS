#include "matmodel.h"

#include <Eigen/Dense>
using Eigen::MatrixXd;
using namespace Eigen;


ThermalLinear :: ThermalLinear(vector<double> l_e, map<string, vector<double>> materialProperties) : ThermalModel(l_e){

    conductivity = materialProperties["conductivity"];
    n_mat = conductivity.size();

    double kappa_ref = (*max_element(conductivity.begin(), conductivity.end()) + *min_element(conductivity.begin(), conductivity.end())) / 2;
    kapparef_mat = kappa_ref * Matrix3d::Identity();

    Matrix3d phase_kappa;
    phase_stiffness = new Matrix<double, 8, 8>[n_mat];

    for (size_t i = 0; i < n_mat; ++i) {
        phase_stiffness[i] = Matrix<double, 8, 8>::Zero();
        phase_kappa = conductivity[i] * Matrix3d::Identity();
    
        for (int p = 0; p < 8; ++p) {
            phase_stiffness[i] += B_int[p].transpose() * phase_kappa * B_int[p] * v_e * 0.1250;
        }
    }
}

MechLinear :: MechLinear(vector<double> l_e, map<string, vector<double>> materialProperties) : MechModel(l_e) {
    
    vector<double> young_modulus = materialProperties["young_modulus"];
    vector<double> poisson_ratio = materialProperties["poisson_ratio"];
    n_mat = young_modulus.size();
    lambda.resize(n_mat);
    mu.resize(n_mat);

    for (int i = 0; i < n_mat; ++i) {
        lambda[i] = (young_modulus[i]*poisson_ratio[i])/((1+poisson_ratio[i])*(1-(2*poisson_ratio[i])));
        mu[i]     = young_modulus[i]/(2*(1+poisson_ratio[i]));
        // if( verbosity > 8 ) printf("lambda[%i] = %16.8e   /   mu[%i] = %16.8e\n", i, lambda[i], i, mu[i]);
    }

    
    double lambda_ref = (*max_element(lambda.begin(), lambda.end()) + *min_element(lambda.begin(), lambda.end())) / 2;
    double mu_ref     = (*max_element(mu.begin(), mu.end()) + *min_element(mu.begin(), mu.end())) / 2;

    // if( verbosity > 8 ) printf("# lambda_ref: %16.8e, mu_ref: %16.8e\n", lambda_ref, mu_ref);


    kapparef_mat = Matrix<double, 6, 6>::Zero();
    kapparef_mat.topLeftCorner(3, 3).setConstant(lambda_ref);
    kapparef_mat += 2*mu_ref * Matrix<double, 6, 6>::Identity();

    phase_stiffness = new Matrix<double, 24, 24>[n_mat];
    Matrix<double, 6, 6> phase_kappa;

    for(int i = 0; i < n_mat; i++){
        phase_kappa.setZero();
        phase_stiffness[i] = Matrix<double, 24, 24>::Zero();

        phase_kappa.topLeftCorner(3, 3).setConstant(lambda[i]);
        phase_kappa += 2*mu[i] * Matrix<double, 6, 6>::Identity();
        
        for (int p = 0; p < 8; ++p) {
            phase_stiffness[i] += B_int[p].transpose() * phase_kappa * B_int[p] * v_e * 0.1250;
        }
    }    
}


HyperElastic :: HyperElastic(vector<double> l_e, map<string, vector<double>> materialProperties) : MechModel(l_e) {
    
    bulk_modulus = materialProperties["bulk_modulus"];
    shear_modulus = materialProperties["shear_modulus"];
    critical_stress = materialProperties["critical_stress"];
    hardening_parameter = materialProperties["hardening_parameter"];

    n_mat = bulk_modulus.size();
    eps_crit.resize(n_mat);
    E_s.resize(n_mat);
    Matrix<double, 6, 6>* Ce = new Matrix<double, 6, 6>[n_mat];

    Matrix<double, 6, 6> topLeft = Matrix<double, 6, 6>::Zero();
    topLeft.topLeftCorner(3, 3).setConstant(1);

    for (int i = 0; i < n_mat; ++i) {
        eps_crit[i] = pow(2./3.,.5)*critical_stress[i]/(2.*shear_modulus[i]);
        E_s[i]      = (3.*shear_modulus[i])/(3.*shear_modulus[i] + hardening_parameter[i]);

        Ce[i] = 3*bulk_modulus[i] * topLeft + 2*shear_modulus[i] * (-1./3 * topLeft + Matrix<double, 6, 6>::Identity());
    }

    kapparef_mat = 0.5*(Ce[0] + Ce[1]); //TODO: only works for two materials
}