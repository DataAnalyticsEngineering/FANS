[Tensors]
  [a]
    type = Python
    expr = 'Scalar(torch.tensor(1.0, dtype=torch.float64))'
  []
  [sdirs]
    type = Python
    expr = 'MillerIndex.fill(1, 1, 0)'
  []
  [splanes]
    type = Python
    expr = 'MillerIndex.fill(1, 1, 1)'
  []
  [no_rotation]
    type = Python
    expr = 'MRP(torch.zeros(3, dtype=torch.float64))'
  []
[]

[Solvers]
  [newton]
    type = Newton
    abs_tol = 1e-12
    rel_tol = 0
    max_its = 100
  []
[]

[Data]
  [crystal_geometry]
    type = CubicCrystal
    lattice_parameter = 'a'
    slip_directions = 'sdirs'
    slip_planes = 'splanes'
  []
[]

[Models]
  [elastic_strain]
    type = SR2LinearCombination
    from = 'forces/E state/internal/plastic_strain'
    to = 'state/elastic_strain'
    weights = '1 -1'
  []
  [elastic_strain_full]
    type = SR2ToR2
    input = 'state/elastic_strain'
    output = 'state/internal/elastic_strain_full'
  []
  [QT_elastic_strain]
    type = R2Multiplication
    A = 'orientation'
    B = 'state/internal/elastic_strain_full'
    to = 'state/internal/QT_elastic_strain'
    transpose_A = true
  []
  [QT_elastic_strain_Q]
    type = R2Multiplication
    A = 'state/internal/QT_elastic_strain'
    B = 'orientation'
    to = 'state/internal/QT_elastic_strain_Q'
  []
  [crystal_elastic_strain]
    type = R2ToSR2
    input = 'state/internal/QT_elastic_strain_Q'
    output = 'state/internal/crystal_elastic_strain'
  []
  [identity]
    type = MRPConstantParameter
    value = 'no_rotation'
    parameter = 'state/internal/no_rotation'
  []
  [elastic_tensor]
    type = CubicElasticityTensor
    coefficients = '100 0.25 50'
    coefficient_types = 'YOUNGS_MODULUS POISSONS_RATIO SHEAR_MODULUS'
  []
  [elasticity]
    type = GeneralElasticity
    elastic_stiffness_tensor = 'elastic_tensor'
    strain = 'state/internal/crystal_elastic_strain'
    orientation = 'state/internal/no_rotation'
    stress = 'state/internal/crystal_cauchy_stress'
  []
  [crystal_cauchy_stress_full]
    type = SR2ToR2
    input = 'state/internal/crystal_cauchy_stress'
    output = 'state/internal/crystal_cauchy_stress_full'
  []
  [Q_cauchy_stress]
    type = R2Multiplication
    A = 'orientation'
    B = 'state/internal/crystal_cauchy_stress_full'
    to = 'state/internal/Q_cauchy_stress'
  []
  [Q_cauchy_stress_QT]
    type = R2Multiplication
    A = 'state/internal/Q_cauchy_stress'
    B = 'orientation'
    to = 'state/internal/Q_cauchy_stress_QT'
    transpose_B = true
  []
  [cauchy_stress]
    type = R2ToSR2
    input = 'state/internal/Q_cauchy_stress_QT'
    output = 'state/S'
  []
  [resolved_shear]
    type = ResolvedShear
    stress = 'state/S'
    orientation_matrix = 'orientation'
    resolved_shears = 'state/internal/resolved_shears'
  []
  [plastic_deformation_rate]
    type = PlasticDeformationRate
    orientation_matrix = 'orientation'
    slip_rates = 'state/internal/slip_rates'
    plastic_deformation_rate = 'state/internal/plastic_strain_rate'
  []
  [sum_slip_rates]
    type = SumSlipRates
    slip_rates = 'state/internal/slip_rates'
    sum_slip_rates = 'state/internal/sum_slip_rates'
  []
  [slip_rule]
    type = PowerLawSlipRule
    resolved_shears = 'state/internal/resolved_shears'
    slip_strengths = 'state/internal/slip_strengths'
    slip_rates = 'state/internal/slip_rates'
    n = 8.0
    gamma0 = 2.0e-1
  []
  [slip_strength]
    type = SingleSlipStrengthMap
    constant_strength = 0.05
    slip_hardening = 'state/internal/slip_hardening'
    slip_strengths = 'state/internal/slip_strengths'
  []
  [voce_hardening]
    type = VoceSingleSlipHardeningRule
    initial_slope = 0.5
    saturated_hardening = 0.05
    slip_hardening = 'state/internal/slip_hardening'
    sum_slip_rates = 'state/internal/sum_slip_rates'
    slip_hardening_rate = 'state/internal/slip_hardening_rate'
  []
  [integrate_slip_hardening]
    type = ScalarBackwardEulerTimeIntegration
    variable = 'state/internal/slip_hardening'
  []
  [integrate_plastic_strain]
    type = SR2BackwardEulerTimeIntegration
    variable = 'state/internal/plastic_strain'
  []
  [implicit_rate]
    type = ComposedModel
    models = 'elastic_strain elastic_strain_full QT_elastic_strain QT_elastic_strain_Q crystal_elastic_strain identity
              elasticity crystal_cauchy_stress_full Q_cauchy_stress Q_cauchy_stress_QT cauchy_stress resolved_shear
              plastic_deformation_rate sum_slip_rates slip_rule slip_strength voce_hardening
              integrate_slip_hardening integrate_plastic_strain'
  []
  [model]
    type = ImplicitUpdate
    equation_system = 'system'
    solver = 'newton'
  []
  [cp]
    type = ComposedModel
    models = 'model elastic_strain elastic_strain_full QT_elastic_strain QT_elastic_strain_Q crystal_elastic_strain
              identity elasticity crystal_cauchy_stress_full Q_cauchy_stress Q_cauchy_stress_QT cauchy_stress'
    additional_outputs = 'state/internal/plastic_strain state/internal/slip_hardening'
  []
[]

[EquationSystems]
  [system]
    type = NonlinearSystem
    model = 'implicit_rate'
    unknowns = 'state/internal/plastic_strain state/internal/slip_hardening'
    residuals = 'state/internal/plastic_strain_residual state/internal/slip_hardening_residual'
  []
[]
