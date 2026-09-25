[Models]
  [Ee_tr]
    type = SR2LinearCombination
    from = 'forces/E state/internal/Ep~1'
    to = 'state/internal/Ee_tr'
    weights = '1 -1'
  []
  [S_tr]
    type = LinearIsotropicElasticity
    coefficients = '75.00010399997781 0.2999997226667258'
    coefficient_types = 'YOUNGS_MODULUS POISSONS_RATIO'
    strain = 'state/internal/Ee_tr'
    stress = 'state/internal/S_tr'
  []
  [X_old]
    type = LinearKinematicHardening
    kinematic_plastic_strain = 'state/internal/Ep~1'
    back_stress = 'state/internal/X'
    hardening_modulus = 2.0
  []
  [O_tr]
    type = SR2LinearCombination
    from = 'state/internal/S_tr state/internal/X'
    to = 'state/internal/O_tr'
    weights = '1 -1'
  []
  [vm_tr]
    type = SR2Invariant
    tensor = 'state/internal/O_tr'
    invariant = 'state/internal/sm_tr'
    invariant_type = 'VONMISES'
  []
  [N]
    type = AssociativeJ2FlowDirection
    mandel_stress = 'state/internal/O_tr'
    flow_direction = 'state/internal/N'
  []

  [dep]
    type = ScalarLinearCombination
    from = 'state/internal/ep state/internal/ep~1'
    to = 'state/internal/dep'
    weights = '1 -1'
  []
  [sm]
    type = ScalarLinearCombination
    from = 'state/internal/sm_tr state/internal/ep state/internal/ep~1'
    to = 'state/internal/sm'
    weights = '1 -89.5386 89.5386'
  []
  [isoharden]
    type = LinearIsotropicHardening
    equivalent_plastic_strain = 'state/internal/ep'
    isotropic_hardening = 'state/internal/k'
    hardening_modulus = 5.0
  []
  [yieldfn]
    type = YieldFunction
    effective_stress = 'state/internal/sm'
    isotropic_hardening = 'state/internal/k'
    yield_function = 'state/internal/fp'
    yield_stress = 0.1
  []
  [consistency]
    type = MinMapComplementarity
    a = 'state/internal/fp'
    b = 'state/internal/dep'
    complementarity = 'residual/internal/ep'
  []

  [residuals]
    type = ComposedModel
    models = 'Ee_tr S_tr X_old O_tr vm_tr dep sm isoharden yieldfn consistency'
  []
  [guess]
    type = ConstantExtrapolationPredictor
    unknowns_Scalar = 'state/internal/ep'
  []
  [solve]
    type = ImplicitUpdate
    equation_system = 'system'
    solver = 'newton'
    predictor = 'guess'
  []

  [dEp]
    type = AssociativePlasticFlow
    flow_rate = 'state/internal/dep'
    flow_direction = 'state/internal/N'
    plastic_strain_rate = 'state/internal/dEp'
  []
  [Ep]
    type = SR2LinearCombination
    from = 'state/internal/Ep~1 state/internal/dEp'
    to = 'state/internal/Ep'
    weights = '1 1'
  []
  [Ee]
    type = SR2LinearCombination
    from = 'forces/E state/internal/Ep'
    to = 'state/internal/Ee'
    weights = '1 -1'
  []
  [S]
    type = LinearIsotropicElasticity
    coefficients = '75.00010399997781 0.2999997226667258'
    coefficient_types = 'YOUNGS_MODULUS POISSONS_RATIO'
    strain = 'state/internal/Ee'
    stress = 'state/S'
  []
  [j2]
    type = ComposedModel
    models = 'solve Ee_tr S_tr X_old O_tr N dep dEp Ep Ee S'
    additional_outputs = 'state/internal/Ep state/internal/ep'
  []
[]

[EquationSystems]
  [system]
    type = NonlinearSystem
    model = 'residuals'
    unknowns = 'state/internal/ep'
    residuals = 'residual/internal/ep'
  []
[]

[Solvers]
  [newton]
    type = Newton
    abs_tol = 1e-10
    rel_tol = 1e-8
    max_its = 100
  []
[]
