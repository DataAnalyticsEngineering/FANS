#include "solverCG.h"
#include "solverFP.h"

// Thermal models
#include "material_models/LinearThermalIsotropic.h"

// Mechanical models
#include "material_models/LinearElasticIsotropic.h"
#include "material_models/PseudoPlastic.h"
#include "material_models/VonMisesPlasticLinearIsotropicHardening.h"
#include "material_models/RateIndependentJ2Plasticity.h"

template <int howmany>
Matmodel<howmany> *createMatmodel(const Reader &reader);

template <>
Matmodel<1> *createMatmodel(const Reader &reader)
{
    if (reader.matmodel == "LinearThermalIsotropic") {
        return new LinearThermalIsotropic(reader.l_e, reader.materialProperties);
    } else {
        throw std::invalid_argument(reader.matmodel + " is not a valid matmodel for thermal problem");
    }
}

template <>
Matmodel<3> *createMatmodel(const Reader &reader)
{
    if (reader.matmodel == "LinearElasticIsotropic") {
        return new LinearElasticIsotropic(reader.l_e, reader.materialProperties);

    // PseudoPlastic models
    } else if (reader.matmodel == "PseudoPlasticLinearHardening") {
        return new PseudoPlasticLinearHardening(reader.l_e, reader.materialProperties);
    } else if (reader.matmodel == "PseudoPlasticNonLinearHardening") {
        return new PseudoPlasticNonLinearHardening(reader.l_e, reader.materialProperties);


    // J2 Plastic models
    } else if (reader.matmodel == "VonMisesPlasticLinearIsotropicHardening") {
        return new VonMisesPlasticLinearIsotropicHardening(reader.l_e, reader.materialProperties);


    } else if (reader.matmodel == "RateIndependentJ2PlasticityLinearIsotropicHardening") {
        return new RateIndependentJ2PlasticityLinearIsotropicHardening(reader.l_e, reader.materialProperties);
    } else if (reader.matmodel == "RateIndependentJ2PlasticityNonLinearIsotropicHardening") {
        return new RateIndependentJ2PlasticityNonLinearIsotropicHardening(reader.l_e, reader.materialProperties);

    } else {
        throw std::invalid_argument(reader.matmodel + " is not a valid matmodel for mechanical problem");
    }
}

template <int howmany>
Solver<howmany> *createSolver(const Reader &reader, Matmodel<howmany> *matmodel)
{
    if (reader.method == "fp") {
        return new SolverFP<howmany>(reader, matmodel);
    } else if (reader.method == "cg") {
        return new SolverCG<howmany>(reader, matmodel);
    } else {
        throw std::invalid_argument(reader.method + " is not a valid method");
    }
}
