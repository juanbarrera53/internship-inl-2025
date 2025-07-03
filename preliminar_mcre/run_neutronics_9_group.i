################################################################################
## Griffin Main Application input file                                        ##
## Steady state neutronics model                                              ##
## Neutron diffusion with delayed precursor source, no equivalence            ##
################################################################################

!include 'template_common/common_settings_neutronics.i'

[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = 'mcre_mesh.e'
    #coord_type = 'XYZ'
  []
  [reactor_out_bdry]
    type = ParsedGenerateSideset
    input = 'fmg'
    new_sideset_name = 'reactor_out'
    combinatorial_geometry = '(y > 775) | ((y <= 775) & (x > 1150 | x < -250))'
    include_only_external_sides = true
  []
  [scaling]
    type = TransformGenerator
    input = reactor_out_bdry
    transform = 'SCALE'
    vector_value = '0.001 0.001 0.001'
  []
[]

[TransportSystems]
  particle = neutron
  equation_type = eigenvalue

  G = 9

  VacuumBoundary = 'wall-reflector reactor_out'

  [diff]
    scheme = CFEM-Diffusion
    n_delay_groups = 6
    external_dnp_variable = 'dnp'
    family = LAGRANGE
    order = FIRST
    fission_source_aux = true

    # For PJFNKMO
    assemble_scattering_jacobian = true
    assemble_fission_jacobian = true

    initialize_scalar_flux = true
  []
[]

[AuxVariables]
  # The dnp variable needs to be checked in a null-transient. In principle, it should not be necessary to restart it,
  # and a restart of each group and the AuxKernel should reconstruct the dnp value
  [dnp]
    order = CONSTANT
    family = MONOMIAL
    components = 6
    block = 'reactor pump pipe'
  []
[]

################################################################################
# AUXILIARY SYSTEM
################################################################################

[ICs]
  [tfuel_ic]
    type = ConstantIC
    variable = tfuel
    value = 900.0
  []
  [trefl_ic]
    type = ConstantIC
    variable = trefl
    value = 760.0
  []
[]

################################################################################
# EXECUTION / SOLVE
################################################################################

[Executioner]
  type = Eigenvalue
  solve_type = PJFNKMO

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart '
  petsc_options_value = 'hypre boomeramg 50'
  #l_max_its = 100

  free_power_iterations = 5 # important to obtain fundamental mode eigenvalue

  nl_abs_tol = 1e-9

  # # MultiApp Parameters
  # fixed_point_abs_tol = 1e-3
  # fixed_point_algorithm = picard
  # fixed_point_min_its = 4
  # fixed_point_max_its = 30
  # relaxation_factor = 0.7

  normalization = fission_source_integral

[]

################################################################################
# POST-PROCESSORS
################################################################################

# [Postprocessors]
#   [dt_limit]
#     type = Receiver
#     value = 0.1
#   []
# []

################################################################################
# MULTIPHYSICS
################################################################################

[MultiApps]
  [ns]
    type = FullSolveMultiApp
    input_files = 'run_ns_initial_res.i'
    execute_on = 'timestep_end'
    no_backup_and_restore = true
    keep_solution_during_restore = true
  []
[]

[Transfers]                                                       #Originally On
  [power_density]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = ns
    source_variable = power_density
    variable = power_density
    from_postprocessors_to_be_preserved = power                 #Originally Commented Out
    to_postprocessors_to_be_preserved = reactor-power           #Originally Commented Out
    error_on_miss = false
  []
[]

# [Transfers]
#   [power_density]
#     type = MultiAppGeneralFieldNearestLocationTransfer
#     to_multi_app = ns
#     source_variable = power_density
#     variable = power_density
#   []
# []
