################################################################################
## Molten Salt Fast Reactor - Euratom EVOL + Rosatom MARS Design              ##
## Pronghorn input file to initialize velocity fields                         ##
## This runs a slow relaxation to steady state while ramping down the fluid   ##
## viscosity.                                                                 ##
################################################################################

!include 'template_common/common_settings_ns.i'

################################################################################
# GEOMETRY
################################################################################

[Mesh]
  [fmg]
    type = FileMeshGenerator
    use_for_exodus_restart = true
    file = 'run_ns_initial_res_restart.e'
  []
[]

[Problem]
  allow_initial_conditions_with_restart = true
[]

################################################################################
# EQUATIONS: VARIABLES, KERNELS & BCS
################################################################################

[Variables]
  [pressure]
    initial_from_file_var = pressure
  []
  [superficial_vel_x]
    initial_from_file_var = superficial_vel_x
  []
  [superficial_vel_y]
    initial_from_file_var = superficial_vel_y
  []
  [superficial_vel_z]
    initial_from_file_var = superficial_vel_z
  []
  [T]
    #initial_from_file_var = T
    initial_condition = 900
  []
  [T_ref]
    initial_from_file_var = T_ref
  []
[]

[FVKernels]
  # inactive = 'u_time v_time w_time heat_time heat_time_ref'

  [heat_src]
    type = FVBodyForce
    variable = T
    function = cosine_guess
    value = '${fparse power/0.21757}' #Volume integral of cosine shape is 0.21757
    block = 'reactor'
  []
[]

################################################################################
# EXECUTION / SOLVE
################################################################################

[Executioner]
  type = Transient

  # Time-stepping parameters
  start_time = 1e7
  end_time = 2e7
  # dt = 1e6

  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 10
    dt = 10
    #dt = 0.01
    timestep_limiting_postprocessor = 'dt_limit'
  []

  # [TimeStepper]
  #   type = SolutionTimeAdaptiveDT
  #   dt = 0.1
  # []

  # [TimeStepper]
  #   type = FunctionDT
  #   function = dts
  # []

  # Solver parameters
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -ksp_gmres_restart'
  petsc_options_value = 'lu NONZERO 20'
  line_search = 'none'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-4
  nl_max_its = 10 # fail early and try again with a shorter time step
  l_max_its = 80
  automatic_scaling = true

  # MultiApp
  relaxation_factor = 1.0

[]

################################################################################
# MULTIAPPS and TRANSFERS for precursors transport
################################################################################

[MultiApps]
  [prec_transport]
    type = TransientMultiApp
    input_files = 'run_prec_transport.i'
    execute_on = 'timestep_end'
    #no_backup_and_restore = true
    sub_cycling = true
  []
[]
