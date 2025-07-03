################################################################################
## Griffin Main Application input file                                        ##
## Steady state neutronics model                                              ##
## Neutron diffusion with delayed precursor source, no equivalence            ##
################################################################################

[GlobalParams]
  library_name = 'serpent_MCRE_xs'
  library_file = 'serpent_MCRE_xs_new.xml'
  is_meter = true
  plus = true
  isotopes = 'pseudo'
  densities = '1.0'
  # No displacement modeled
  # fixed_meshes = true
[]

# [TransportSystems]
#   particle = neutron
#   equation_type = transient
#   restart_transport_system = true
#   scaling_eigenkernels = 1.16922 # To correct initialization drift 1.1675

#   G = 9

#   VacuumBoundary = 'wall-reflector reactor_out'

#   [diff]
#     scheme = CFEM-Diffusion
#     n_delay_groups = 6
#     external_dnp_variable = 'dnp'
#     family = LAGRANGE
#     order = FIRST
#     fission_source_aux = true
#   []
# []

[PowerDensity]
  #power = 25e3
  power = 1e6
  power_density_variable = power_density
  integrated_power_postprocessor = integrated_power
  power_scaling_postprocessor = power_scaling
  family = MONOMIAL
  order = CONSTANT
  #execute_on = 'transfer'
[]

################################################################################
# AUXILIARY SYSTEM
################################################################################

[AuxVariables]
  [tfuel]
    order = CONSTANT
    family = MONOMIAL
    block = 'reactor pump pipe'
  []
  [trefl]
    order = CONSTANT
    family = MONOMIAL
    block = 'reflector'
  []
  [c1]
    order = CONSTANT
    family = MONOMIAL
    block = 'reactor pump pipe'
  []
  [c2]
    order = CONSTANT
    family = MONOMIAL
    block = 'reactor pump pipe'
  []
  [c3]
    order = CONSTANT
    family = MONOMIAL
    block = 'reactor pump pipe'
  []
  [c4]
    order = CONSTANT
    family = MONOMIAL
    block = 'reactor pump pipe'
  []
  [c5]
    order = CONSTANT
    family = MONOMIAL
    block = 'reactor pump pipe'
  []
  [c6]
    order = CONSTANT
    family = MONOMIAL
    block = 'reactor pump pipe'
  []

[]

[AuxKernels]
  [build_dnp]
    type = BuildArrayVariableAux
    variable = dnp
    component_variables = 'c1 c2 c3 c4 c5 c6'
    execute_on = 'timestep_begin final'
  []
[]

################################################################################
# CROSS SECTIONS
################################################################################

[Materials]
  [reactor_material]
    type = CoupledFeedbackNeutronicsMaterial
    grid_names = 'tfuel'
    grid_variables = 'tfuel'
    material_id = 1
    block = 'reactor pump pipe'
  []
  [reflector_material]
    type = CoupledFeedbackNeutronicsMaterial
    grid_names = 'trefl'
    grid_variables = 'trefl'
    material_id = 2
    block = 'reflector'
  []
[]

[Debug]
  show_var_residual_norms = true
[]

################################################################################
# POST-PROCESSORS
################################################################################

[Postprocessors]
  [power]
    type = ElementIntegralVariablePostprocessor
    variable = power_density
    execute_on = 'initial timestep_begin timestep_end transfer'
    outputs = all
    block = 'reactor pump pipe'
  []

  [flux_0_average_reac]
    type = ElementAverageValue
    variable = sflux_g0
    block = 'reactor'
    outputs = out
  []
  [flux_1_average_reac]
    type = ElementAverageValue
    variable = sflux_g1
    block = 'reactor'
    outputs = out
  []
  [flux_2_average_reac]
    type = ElementAverageValue
    variable = sflux_g2
    block = 'reactor'
    outputs = out
  []
  [flux_3_average_reac]
    type = ElementAverageValue
    variable = sflux_g3
    block = 'reactor'
    outputs = out
  []
  [flux_4_average_reac]
    type = ElementAverageValue
    variable = sflux_g4
    block = 'reactor'
    outputs = out
  []
  [flux_5_average_reac]
    type = ElementAverageValue
    variable = sflux_g5
    block = 'reactor'
    outputs = out
  []
  [flux_6_average_reac]
    type = ElementAverageValue
    variable = sflux_g6
    block = 'reactor'
    outputs = out
  []
  [flux_7_average_reac]
    type = ElementAverageValue
    variable = sflux_g7
    block = 'reactor'
    outputs = out
  []
  [flux_0_average_ref]
    type = ElementAverageValue
    variable = sflux_g0
    block = 'reactor'
    outputs = out
  []
  [flux_1_average_ref]
    type = ElementAverageValue
    variable = sflux_g1
    block = 'reflector'
    outputs = out
  []
  [flux_2_average_ref]
    type = ElementAverageValue
    variable = sflux_g2
    block = 'reflector'
    outputs = out
  []
  [flux_3_average_ref]
    type = ElementAverageValue
    variable = sflux_g3
    block = 'reflector'
    outputs = out
  []
  [flux_4_average_ref]
    type = ElementAverageValue
    variable = sflux_g4
    block = 'reflector'
    outputs = out
  []
  [flux_5_average_ref]
    type = ElementAverageValue
    variable = sflux_g5
    block = 'reflector'
    outputs = out
  []
  [flux_6_average_ref]
    type = ElementAverageValue
    variable = sflux_g6
    block = 'reflector'
    outputs = out
  []
  [flux_7_average_ref]
    type = ElementAverageValue
    variable = sflux_g7
    block = 'reflector'
    outputs = out
  []
  [flux_0_average_reac_max]
    type = ElementExtremeValue
    variable = sflux_g0
    block = 'reactor'
    outputs = out
  []
  [flux_1_average_reac_max]
    type = ElementExtremeValue
    variable = sflux_g1
    block = 'reactor'
    outputs = out
  []
  [flux_2_average_reac_max]
    type = ElementExtremeValue
    variable = sflux_g2
    block = 'reactor'
    outputs = out
  []
  [flux_3_average_reac_max]
    type = ElementExtremeValue
    variable = sflux_g3
    block = 'reactor'
    outputs = out
  []
  [flux_4_average_reac_max]
    type = ElementExtremeValue
    variable = sflux_g4
    block = 'reactor'
    outputs = out
  []
  [flux_5_average_reac_max]
    type = ElementExtremeValue
    variable = sflux_g5
    block = 'reactor'
    outputs = out
  []
  [flux_6_average_reac_max]
    type = ElementExtremeValue
    variable = sflux_g6
    block = 'reactor'
    outputs = out
  []
  [flux_7_average_reac_max]
    type = ElementExtremeValue
    variable = sflux_g7
    block = 'reactor'
    outputs = out
  []
  [flux_0_average_ref_max]
    type = ElementExtremeValue
    variable = sflux_g0
    block = 'reactor'
    outputs = out
  []
  [flux_1_average_ref_max]
    type = ElementExtremeValue
    variable = sflux_g1
    block = 'reflector'
    outputs = out
  []
  [flux_2_average_ref_max]
    type = ElementExtremeValue
    variable = sflux_g2
    block = 'reflector'
    outputs = out
  []
  [flux_3_average_ref_max]
    type = ElementExtremeValue
    variable = sflux_g3
    block = 'reflector'
    outputs = out
  []
  [flux_4_average_ref_max]
    type = ElementExtremeValue
    variable = sflux_g4
    block = 'reflector'
    outputs = out
  []
  [flux_5_average_ref_max]
    type = ElementExtremeValue
    variable = sflux_g5
    block = 'reflector'
    outputs = out
  []
  [flux_6_average_ref_max]
    type = ElementExtremeValue
    variable = sflux_g6
    block = 'reflector'
    outputs = out
  []
  [flux_7_average_ref_max]
    type = ElementExtremeValue
    variable = sflux_g7
    block = 'reflector'
    outputs = out
  []
  [c1_average_reac]
    type = ElementAverageValue
    variable = c1
    block = 'reactor'
    outputs = out
  []
  [c2_average_reac]
    type = ElementAverageValue
    variable = c2
    block = 'reactor'
    outputs = out
  []
  [c3_average_reac]
    type = ElementAverageValue
    variable = c3
    block = 'reactor'
    outputs = out
  []
  [c4_average_reac]
    type = ElementAverageValue
    variable = c4
    block = 'reactor'
    outputs = out
  []
  [c5_average_reac]
    type = ElementAverageValue
    variable = c5
    block = 'reactor'
    outputs = out
  []
  [c6_average_reac]
    type = ElementAverageValue
    variable = c6
    block = 'reactor'
    outputs = out
  []
  [c1_average_pipe]
    type = ElementAverageValue
    variable = c1
    block = 'pipe pump'
    outputs = out
  []
  [c2_average_pipe]
    type = ElementAverageValue
    variable = c2
    block = 'pipe pump'
    outputs = out
  []
  [c3_average_pipe]
    type = ElementAverageValue
    variable = c3
    block = 'pipe pump'
    outputs = out
  []
  [c4_average_pipe]
    type = ElementAverageValue
    variable = c4
    block = 'pipe pump'
    outputs = out
  []
  [c5_average_pipe]
    type = ElementAverageValue
    variable = c5
    block = 'pipe pump'
    outputs = out
  []
  [c6_average_pipe]
    type = ElementAverageValue
    variable = c6
    block = 'pipe pump'
    outputs = out
  []
[]

# [VectorPostprocessors]
#   [vpp]
#     type = VectorOfPostprocessors
#     postprocessors = 'reactor-power reactor-power-func
#                       flow max_v pdrop-mixing-plate
#                       heat-loss-pipe heat-loss-ref heat-input-pipe-pp heat-input-ref-pp
#                       T_reac_bot T_reac_top T_reac_ave T_reac_max T_pipe_ave
#                       T_ref_wall T_ref_ave T_ref_max
#                       c1_average_reac c2_average_reac c3_average_reac c4_average_reac c5_average_reac c6_average_reac
#                       c1_average_pipe c2_average_pipe c3_average_pipe c4_average_pipe c5_average_pipe c6_average_pipe'
#   []
# []

################################################################################
# SIMULATION OUTPUTS
################################################################################

[Outputs]
  exodus = true
  [out]
    type = CSV
    execute_on = 'initial timestep_end final'
  []
  checkpoint = true
  # [restart]                             #Originally Implemented
  #   type = Exodus
  #   execute_on = 'final'
  #   file_base = 'run_neutronics_restart'    #FileNotAvailable
  # []
  print_linear_converged_reason = false
  print_linear_residuals = false
  print_nonlinear_converged_reason = false
  # hide = 'dt_limit'
[]

################################################################################
# MULTIPHYSICS
################################################################################

[Transfers]
  [fission_source]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = ns
    source_variable = fission_source
    variable = fission_source
    error_on_miss = false
  []

  [c1]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = ns
    source_variable = 'c1'
    variable = 'c1'
  []
  [c2]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = ns
    source_variable = 'c2'
    variable = 'c2'
  []
  [c3]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = ns
    source_variable = 'c3'
    variable = 'c3'
  []
  [c4]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = ns
    source_variable = 'c4'
    variable = 'c4'
  []
  [c5]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = ns
    source_variable = 'c5'
    variable = 'c5'
  []
  [c6]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = ns
    source_variable = 'c6'
    variable = 'c6'
  []
  [T]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = ns
    source_variable = 'T'
    variable = 'tfuel'
  []
  [T_ref]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = ns
    source_variable = 'T_ref'
    variable = 'trefl'
  []
[]
