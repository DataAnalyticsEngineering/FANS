#include "solverCG.h"
#include "solverFP.h"

// Thermal models
#include "material_models/LinearThermal.h"
#include "material_models/GBDiffusion.h"

// Small strain mechanical models
#include "material_models/LinearElastic.h"
#include "material_models/PseudoPlastic.h"
#include "material_models/J2Plasticity.h"
#include "material_models/J2PlasticityNew.h"

// Large strain mechanical models
#include "material_models/SaintVenantKirchhoff.h"
#include "material_models/CompressibleNeoHookean.h"

template <int howmany, int n_str>
Matmodel<howmany, n_str> *createMatmodel(const Reader &reader);

template <>
Matmodel<1, 3> *createMatmodel<1, 3>(const Reader &reader)
{
    if (reader.matmodel == "LinearThermalIsotropic") {
        return new LinearThermalIsotropic(reader);
    } else if (reader.matmodel == "LinearThermalTriclinic") {
        return new LinearThermalTriclinic(reader);
    } else if (reader.matmodel == "GBDiffusion") {
        return new GBDiffusion(reader);
    } else {
        throw std::invalid_argument(reader.matmodel + " is not a valid matmodel for thermal problem");
    }
}

template <>
Matmodel<3, 6> *createMatmodel<3, 6>(const Reader &reader)
{
    // Linear Elastic models
    if (reader.matmodel == "LinearElasticIsotropic") {
        return new LinearElasticIsotropic(reader);
    } else if (reader.matmodel == "LinearElasticTriclinic") {
        return new LinearElasticTriclinic(reader);

        // Pseudo Plastic models
    } else if (reader.matmodel == "PseudoPlasticLinearHardening") {
        return new PseudoPlasticLinearHardening(reader);
    } else if (reader.matmodel == "PseudoPlasticNonLinearHardening") {
        return new PseudoPlasticNonLinearHardening(reader);

        // J2 Plastic models
    } else if (reader.matmodel == "J2ViscoPlastic_LinearIsotropicHardening") {
        return new J2ViscoPlastic_LinearIsotropicHardening(reader);
    } else if (reader.matmodel == "J2ViscoPlastic_NonLinearIsotropicHardening") {
        return new J2ViscoPlastic_NonLinearIsotropicHardening(reader);
    } else if (reader.matmodel == "J2PlasticityNew_LinearIsotropicHardening") {
        return new J2PlasticityNew_LinearIsotropicHardening(reader);

    } else {
        throw std::invalid_argument(reader.matmodel + " is not a valid small strain material model");
    }
}

template <>
Matmodel<3, 9> *createMatmodel<3, 9>(const Reader &reader)
{
    if (reader.matmodel == "SaintVenantKirchhoff") {
        return new SaintVenantKirchhoff(reader);
    } else if (reader.matmodel == "CompressibleNeoHookean") {
        return new CompressibleNeoHookean(reader);
    } else {
        throw std::invalid_argument(reader.matmodel + " is not a valid large strain material model");
    }
}

template <int howmany, int n_str>
Solver<howmany, n_str> *createSolver(Reader &reader, MaterialManager<howmany, n_str> *matmanager)
{
    if (reader.method == "fp") {
        return new SolverFP<howmany, n_str>(reader, matmanager);
    } else if (reader.method == "cg") {
        return new SolverCG<howmany, n_str>(reader, matmanager);
    } else {
        throw std::invalid_argument(reader.method + " is not a valid method");
    }
}

template <int howmany, int n_str>
MaterialManager<howmany, n_str> *createMaterialManager(const Reader &reader)
{
    return new MaterialManager<howmany, n_str>(reader);
}
