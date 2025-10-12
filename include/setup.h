#include "solverCG.h"
#include "solverFP.h"

// Thermal models
#include "material_models/LinearThermal.h"
#include "material_models/GBDiffusion.h"

// Small strain mechanical models
#include "material_models/LinearElastic.h"
#include "material_models/PseudoPlastic.h"
#include "material_models/J2Plasticity.h"

// Large strain mechanical models
#include "material_models/SaintVenantKirchhoff.h"
#include "material_models/CompressibleNeoHookean.h"

template <int howmany, int n_str>
Matmodel<howmany, n_str> *createMatmodel(const Reader &reader);

template <>
Matmodel<1, 3> *createMatmodel<1, 3>(const Reader &reader)
{
    if (reader.matmodel == "LinearThermalIsotropic") {
        return new LinearThermalIsotropic(reader.l_e, reader.materialProperties);
    } else if (reader.matmodel == "LinearThermalTriclinic") {
        return new LinearThermalTriclinic(reader.l_e, reader.materialProperties);
    } else if (reader.matmodel == "GBDiffusion") {
        return new GBDiffusion(const_cast<Reader &>(reader));
    } else {
        throw std::invalid_argument(reader.matmodel + " is not a valid matmodel for thermal problem");
    }
}

template <>
Matmodel<3, 6> *createMatmodel<3, 6>(const Reader &reader)
{
    // Linear Elastic models
    if (reader.matmodel == "LinearElasticIsotropic") {
        return new LinearElasticIsotropic(reader.l_e, reader.materialProperties);
    } else if (reader.matmodel == "LinearElasticTriclinic") {
        return new LinearElasticTriclinic(reader.l_e, reader.materialProperties);

        // Pseudo Plastic models
    } else if (reader.matmodel == "PseudoPlasticLinearHardening") {
        return new PseudoPlasticLinearHardening(reader.l_e, reader.materialProperties);
    } else if (reader.matmodel == "PseudoPlasticNonLinearHardening") {
        return new PseudoPlasticNonLinearHardening(reader.l_e, reader.materialProperties);

        // J2 Plastic models
    } else if (reader.matmodel == "J2ViscoPlastic_LinearIsotropicHardening") {
        return new J2ViscoPlastic_LinearIsotropicHardening(reader.l_e, reader.materialProperties);
    } else if (reader.matmodel == "J2ViscoPlastic_NonLinearIsotropicHardening") {
        return new J2ViscoPlastic_NonLinearIsotropicHardening(reader.l_e, reader.materialProperties);

    } else {
        throw std::invalid_argument(reader.matmodel + " is not a valid small strain material model");
    }
}

template <>
Matmodel<3, 9> *createMatmodel<3, 9>(const Reader &reader)
{
    if (reader.matmodel == "SaintVenantKirchhoff") {
        return new SaintVenantKirchhoff(reader.l_e, reader.materialProperties);
    } else if (reader.matmodel == "CompressibleNeoHookean") {
        return new CompressibleNeoHookean(reader.l_e, reader.materialProperties);
    } else {
        throw std::invalid_argument(reader.matmodel + " is not a valid large strain material model");
    }
}

template <int howmany, int n_str>
Solver<howmany, n_str> *createSolver(Reader &reader, Matmodel<howmany, n_str> *matmodel)
{
    if (reader.method == "fp") {
        return new SolverFP<howmany, n_str>(reader, matmodel);
    } else if (reader.method == "cg") {
        return new SolverCG<howmany, n_str>(reader, matmodel);
    } else {
        throw std::invalid_argument(reader.method + " is not a valid method");
    }
}
