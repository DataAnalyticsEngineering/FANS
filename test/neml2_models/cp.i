[Tensors]
  [sdirs]
    type = Python
    expr = 'MillerIndex.fill(1, 1, 0)'
  []
  [splanes]
    type = Python
    expr = 'MillerIndex.fill(1, 1, 1)'
  []
  [a]
    type = Python
    expr = 'Scalar(torch.tensor(1.0, dtype=torch.float64))'
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
  [Ee]
    type = SR2LinearCombination
    from = 'forces/E state/internal/Ep'
    to = 'state/internal/Ee'
    weights = '1 -1'
  []
  [elasticity]
    type = LinearIsotropicElasticity
    coefficients = '75.00010399997781 0.2999997226667258'
    coefficient_types = 'YOUNGS_MODULUS POISSONS_RATIO'
    strain = 'state/internal/Ee'
    stress = 'state/S'
  []
  [rss]
    type = ResolvedShear
    stress = 'state/S'
    orientation_matrix = 'orientation'
    resolved_shears = 'state/internal/rss'
  []
  [strength]
    type = ScalarConstantParameter
    value = 0.05
    parameter = 'state/internal/tau_hat'
  []
  [slip]
    type = PowerLawSlipRule
    resolved_shears = 'state/internal/rss'
    slip_strengths = 'state/internal/tau_hat'
    slip_rates = 'state/internal/dgamma'
    gamma0 = 1e-2
    n = 10
  []
  [dEp]
    type = PlasticDeformationRate
    orientation_matrix = 'orientation'
    slip_rates = 'state/internal/dgamma'
    plastic_deformation_rate = 'state/internal/dEp'
  []
  [residual]
    type = SR2LinearCombination
    from = 'state/internal/Ep state/internal/Ep~1 state/internal/dEp'
    to = 'residual/internal/Ep'
    weights = '1 -1 -1'
  []
  [residuals]
    type = ComposedModel
    models = 'Ee elasticity rss strength slip dEp residual'
  []
  [guess]
    type = ConstantExtrapolationPredictor
    unknowns_SR2 = 'state/internal/Ep'
  []
  [solve]
    type = ImplicitUpdate
    equation_system = 'system'
    solver = 'newton'
    predictor = 'guess'
  []
  [cp]
    type = ComposedModel
    models = 'solve Ee elasticity'
    additional_outputs = 'state/internal/Ep'
  []
[]

[EquationSystems]
  [system]
    type = NonlinearSystem
    model = 'residuals'
    unknowns = 'state/internal/Ep'
    residuals = 'residual/internal/Ep'
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
