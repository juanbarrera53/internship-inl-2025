################################################################################
## Molten Salt Fast Reactor - Euratom EVOL + Rosatom MARS Design              ##
## Pronghorn input file to initialize velocity fields                         ##
## This runs a slow relaxation to steady state while ramping down the fluid   ##
## viscosity.                                                                 ##
################################################################################

!include 'template_common/common_settings_prec.i'

################################################################################
# GEOMETRY
################################################################################

[Mesh]
  [restart]
    type = FileMeshGenerator
    use_for_exodus_restart = true
    file = 'run_ns_initial_res_restart.e'
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
  #dt = 1

  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 10
    dt = 1
  []

  # [TimeStepper]
  #   type = FunctionDT
  #   function = dts
  # []

  # Solver parameters
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -ksp_gmres_restart'
  petsc_options_value = 'lu NONZERO 20'
  line_search = 'none'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  nl_max_its = 10 # fail early and try again with a shorter time step
  l_max_its = 80
  automatic_scaling = true
[]
