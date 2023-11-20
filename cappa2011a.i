# Working
# Goals:
# 1. implement Hydro model only
# 2. implement mechanical part
# 3. implement plastic model with mohr coulomb


gravity = -9.81
solid_density = 2260
fluid_density = 1059
porosity0 = 0.1

[Mesh]
  [file_mesh]
    type = FileMeshGenerator
    file = 'cappa2011_simple.msh'
  []
[]

[GlobalParams]
  PorousFlowDictator = 'dictator'
  gravity = '0 -9.81 0'
  displacements = 'disp_x disp_y'
[]

[AuxVariables]
  [xnacl]
    initial_condition = 0.01
  []
  [temp0]
    order = CONSTANT
    family = MONOMIAL
  []
  [massfrac_ph0_sp0]
    initial_condition = 1 # all H20 in phase=0
  []
  [massfrac_ph1_sp0]
    initial_condition = 0 # no H2O in phase=1
  []
[]

[AuxKernels]
  [temp]
    type = FunctionAux
    variable = temp0
    function = 'temp_ic'
    execute_on = 'INITIAL'
  []
[]

[Variables]
  [pwater]
    #scaling = 1E6
    [InitialCondition]
      type = FunctionIC
      function = 'p0'
    []
  []
  [sgas]
    initial_condition = 0.0
  []
[]
[ICs]
  [temperature_ic]
    type = FunctionIC
    function = 'temp_ic'
    variable = temp0
  []
[]

[Modules/TensorMechanics/Master]
  [./plane_strain]
    strain = SMALL
    out_of_plane_direction = z
    add_variables = true
    planar_formulation = PLANE_STRAIN
    #eigenstrain_names = 'ini_strain'
    scaling = 1e-10
    generate_output = 'stress_xx stress_xy stress_yy stress_zz strain_xx strain_xy strain_yy strain_zz'
  [../]
[]


[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pwater
  []
  [flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pwater
  []
  [mass1]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = sgas
  []
  [flux1]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = sgas
  []
  [poro_x]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 1.0
    variable = disp_x
    component = 0
  []
  [poro_y]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 1.0
    variable = disp_y
    component = 1
  []
  [vol_strain_rate_co2]
    type = PorousFlowMassVolumetricExpansion
    fluid_component = 1
    variable = sgas
  []
  [vol_strain_rate_water]
     type = PorousFlowMassVolumetricExpansion
     fluid_component = 0
     variable = pwater
   []

[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pwater sgas disp_x disp_y'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    alpha = 5.099e-5 #1e1
    m = 0.457
    sat_lr = 0.0
    pc_max = 1e7
  []

  # [./mc_coh] #Cohesion
  #   type = TensorMechanicsHardeningConstant
  #   value = 0
  #   convert_to_radians = true
  # [../]
  #
  # [./mc_phi] # Friction angle
  #   type = TensorMechanicsHardeningExponential
  #   value_0 = 0.436332 # 25 deg
  #   value_residual = 0.174533 # 10deg
  #   rate = 1E7
  # [../]
  # [./mc_psi] # Dilation angle
  #   type = TensorMechanicsHardeningConstant
  #   value = 20 # degree
  #   convert_to_radians = true
  # [../]
  # [./mc] # Mohr-Coulomb
  #   type = TensorMechanicsPlasticMohrCoulomb
  #   cohesion = mc_coh
  #   friction_angle = mc_phi
  #   dilation_angle = mc_psi
  #   mc_tip_smoother = 4
  #   mc_edge_smoother = 20
  #   yield_function_tolerance = 1E-3
  #   internal_constraint_tolerance = 1E-9
  # [../]
[]
[FluidProperties]
  [true_water]
    type = Water97FluidProperties
  []
  [tabulated_water]
    type = TabulatedFluidProperties
    fp = true_water
    temperature_min = 273.15
    temperature_max = 426
    pressure_max = 1E8
    interpolated_properties = 'density viscosity enthalpy internal_energy'
    fluid_property_file = water97_tabulated_11.csv
  []
  [true_co2]
    type = CO2FluidProperties
  []
  [tabulated_co2]
    type = TabulatedFluidProperties
    fp = true_co2
    temperature_min = 315
    temperature_max = 326
    pressure_max = 1E8
    interpolated_properties = 'density viscosity enthalpy internal_energy'
    fluid_property_file = co2_tabulated_11.csv
  []
[]

[Materials]
  [eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
  []
  [temperature]
    type = PorousFlowTemperature
    temperature = 'temp0'
  []
  [ppss]
    type = PorousFlow2PhasePS
    phase0_porepressure = pwater
    phase1_saturation = sgas
    capillary_pressure = pc
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'massfrac_ph0_sp0 massfrac_ph1_sp0'
  []
  [brine]
    type = PorousFlowBrine
    water_fp = tabulated_water
    temperature_unit = Kelvin
    xnacl = 0.1047
    phase = 0
  []

  [gas]
    type = PorousFlowSingleComponentFluid
    fp = tabulated_co2
    phase = 1
  []

  [porosity_aquifer]
    type = PorousFlowPorosity
    block = 'aquifer'
    porosity_zero = 0.1
    biot_coefficient = 1.0
    solid_bulk = 1.0 # Required but irrelevant when biot_coefficient is unity
    mechanical = true
  []
  [porosity_cap]
    type = PorousFlowPorosityConst
    block = 'caps'
    porosity = '0.01'
  []
  [permeability_aquifer]
    block = 'aquifer'
    type = PorousFlowPermeabilityKozenyCarman
    poroperm_function = kozeny_carman_phi0
    k0 = 1e-11
    m = 5
    n = 5
    phi0 = 0.1
  []
  [permeability_cap]
    type = PorousFlowPermeabilityConst
    block = 'caps'
    permeability = '1e-19 0 0 0 1e-19 0 0 0 0'
  []
  [relperm_water]
     type = PorousFlowRelativePermeabilityVG
     m = 0.457
     phase = 0
     s_res = 0.3
     sum_s_res = 0.35
     wetting = true
   []
  [relperm_gas]
    type = PorousFlowRelativePermeabilityBC
    phase = 1
    s_res = 0.05
    sum_s_res = 0.35
    nw_phase = true
    lambda = 2
  []
  [elasticity_tensor_other]
     type = ComputeIsotropicElasticityTensor
     youngs_modulus = 10E9 # 10 GPa, 5GPa
     poissons_ratio = 0.25
   []
   # [ini_strain]
   #   type = ComputeEigenstrainFromInitialStress
   #   initial_stress = 'sxx0 0 0  0 syy0 0  0 0 sxx0'
   #   eigenstrain_name = 'ini_strain'
   # []
   [volstrain]
     type = PorousFlowVolumetricStrain
   []
   [stress]
    type = ComputeLinearElasticStress
  []
  #  [./mc] # Compute stress as per Mohr-Coulomb
  #   type = ComputeMultiPlasticityStress
  #   block = 'aquifer'
  #   ep_plastic_tolerance = 1E-10
  #   plastic_models = mc
  #   max_NR_iterations = 1000
  #   debug_fspb = crash
  # [../]
[]

[BCs]
  # Left pgas = No flow
  [injection_area]
    type = PorousFlowSink
    boundary = 'injection_area'
    variable = sgas
    fluid_phase = 1 #gas phase
    flux_function = 'min(t/100.0,1)*(-0.001)' # kg/s/m  #kg/m3 * m/s -> kg/m2/s
    use_relperm = false
  []
  [top_pwater]
    type = DirichletBC
    variable = pwater
    value = 12.853e6  # # as per gradient
    boundary = 'top'
  []
  [bottom_pwater]
    type = DirichletBC
    variable = pwater
    value = 16.777e6  # as per gradient
    boundary = 'bottom'
  []
  [right_pwater]
    type = FunctionDirichletBC
    variable = pwater
    function = p0 # same as initial pressure gradient
    boundary = 'right'
  []

  [right_sgas]
    type = PorousFlowOutflowBC
    variable = sgas
    mass_fraction_component = 1
    boundary = 'right'
  []

  [bottom_roller]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [left_roller]
    type = DirichletBC
    variable = disp_x
    boundary = 'left injection_area'
    value = 0
  []

[]
[Preconditioning]
  active = 'smp'
  [smp]
    type = SMP
    full = true
#    petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
#    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
#    petsc_options_value = 'gmres      asm      lu           NONZERO                   2               1Eq1       1E-3        500'
  []
  [mumps]
    type = SMP
    full = true
    #petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    #petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -snes_rtol -snes_atol -snes_max_it'
    #petsc_options_value = 'gmres      lu       mumps                         NONZERO               1E-4       1E-1       50'
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 3456000
  nl_max_its = 50
  l_max_its = 20
  dtmax = 3600
  l_abs_tol = 1e-4
  nl_abs_tol = 1e-2
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'
  #automatic_scaling = true

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 10
  []
  dtmin = 1
[]

[Functions]
  [p0]
    type = ParsedFunction
    expression = '0.1e6 - 9.81e3 * y' # -ve y 9.81 MPa/km = 9.81e6/1000m = 9.81e3Pa/m
    execute_on = INITIAL
  []
  [temp_ic]
    type = ParsedFunction
    expression =  '273.15+(10 - 25e-3 * y)' # -ve y 25C/km = 25/1000m = 25e-3C/m
    execute_on = INITIAL
  []
  [sxx0]
    type = ParsedFunction
    expression = '-0.7*${gravity} * y * (${solid_density} - ${fluid_density}) * (1.0 - ${porosity0})'  # initial effective stress that should result from weight force
  []
  [syy0]
    type = ParsedFunction
    expression = '-${gravity} * y * (${solid_density} - ${fluid_density}) * (1.0 - ${porosity0})'  # initial effective stress that should result from weight force
  []
[]

[Outputs]
  exodus = true
  interval = 2
[]
