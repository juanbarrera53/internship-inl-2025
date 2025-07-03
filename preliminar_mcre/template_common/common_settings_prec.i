################################################################################
## Molten Salt Fast Reactor - Euratom EVOL + Rosatom MARS Design              ##
## Pronghorn input file to initialize velocity fields                         ##
## This runs a slow relaxation to steady state while ramping down the fluid   ##
## viscosity.                                                                 ##
################################################################################

# Mass flow rate tuning
porosity = 1.0
mu = 0.0166

# Numerical scheme parameters
advected_interp_method = 'upwind'
velocity_interp_method = 'rc'

# Dynamic scaling paramters
mixing_length_pipe_callibrated = '${fparse 0.07 * 0.1 * 0.06}'
mixing_length_reactor_callibrated = '${fparse 0.07 * 0.1 * 2}'
Sc_t = 0.9

# Delayed neutron precursor parameters. Lambda values are decay constants in
# [1 / s]. Beta values are production fractions.
lambda1 = 0.0124667
lambda2 = 0.0282917
lambda3 = 0.0425244
lambda4 = 0.133042
lambda5 = 0.292467
lambda6 = 0.666488
beta1 = 0.000250489
beta2 = 0.00104893
beta3 = 0.000726812
beta4 = 0.00143736
beta5 = 0.0022503
beta6 = 0.000680667

[GlobalParams]
  rhie_chow_user_object = 'pins_rhie_chow_interpolator'

  two_term_boundary_expansion = true
  advected_interp_method = ${advected_interp_method}
  velocity_interp_method = ${velocity_interp_method}
  u = superficial_vel_x
  v = superficial_vel_y
  w = superficial_vel_z
  porosity = porosity
  pressure = pressure

  mixing_length = 'mixing_length'
  schmidt_number = ${Sc_t}

  block = 'reactor pipe pump mixing-plate'
[]

[Problem]
  kernel_coverage_check = false
[]

################################################################################
# EQUATIONS: VARIABLES, KERNELS & BCS
################################################################################

[UserObjects]
  [pins_rhie_chow_interpolator]
    type = PINSFVRhieChowInterpolator
    a_u = a_u_var
    a_v = a_v_var
    a_w = a_w_var
  []
[]

[Variables]
  [c1]
    type = MooseVariableFVReal
  []
  [c2]
    type = MooseVariableFVReal
  []
  [c3]
    type = MooseVariableFVReal
  []
  [c4]
    type = MooseVariableFVReal
  []
  [c5]
    type = MooseVariableFVReal
  []
  [c6]
    type = MooseVariableFVReal
  []
[]

[AuxVariables]
  [power_density]
    type = MooseVariableFVReal
  []
  [fission_source]
    type = MooseVariableFVReal
  []
  [superficial_vel_x]
    type = PINSFVSuperficialVelocityVariable
  []
  [superficial_vel_y]
    type = PINSFVSuperficialVelocityVariable
  []
  [superficial_vel_z]
    type = PINSFVSuperficialVelocityVariable
  []
  [pressure]
    type = INSFVPressureVariable
  []
  [a_u_var]
    type = MooseVariableFVReal
  []
  [a_v_var]
    type = MooseVariableFVReal
  []
  [a_w_var]
    type = MooseVariableFVReal
  []
[]

[FVKernels]
  [c1_advection]
    type = INSFVScalarFieldAdvection
    variable = c1
  []
  [c2_advection]
    type = INSFVScalarFieldAdvection
    variable = c2
  []
  [c3_advection]
    type = INSFVScalarFieldAdvection
    variable = c3
  []
  [c4_advection]
    type = INSFVScalarFieldAdvection
    variable = c4
  []
  [c5_advection]
    type = INSFVScalarFieldAdvection
    variable = c5
  []
  [c6_advection]
    type = INSFVScalarFieldAdvection
    variable = c6
  []
  [c1_diffusion]
    type = FVDiffusion
    coeff = '${mu}'
    variable = c1
  []
  [c2_diffusion]
    type = FVDiffusion
    coeff = '${mu}'
    variable = c2
  []
  [c3_diffusion]
    type = FVDiffusion
    coeff = '${mu}'
    variable = c3
  []
  [c4_diffusion]
    type = FVDiffusion
    coeff = '${mu}'
    variable = c4
  []
  [c5_diffusion]
    type = FVDiffusion
    coeff = '${mu}'
    variable = c5
  []
  [c6_diffusion]
    type = FVDiffusion
    coeff = '${mu}'
    variable = c6
  []
  [c1_diffusion_turb]
    type = INSFVMixingLengthScalarDiffusion
    variable = c1
  []
  [c2_diffusion_turb]
    type = INSFVMixingLengthScalarDiffusion
    variable = c2
  []
  [c3_diffusion_turb]
    type = INSFVMixingLengthScalarDiffusion
    variable = c3
  []
  [c4_diffusion_turb]
    type = INSFVMixingLengthScalarDiffusion
    variable = c4
  []
  [c5_diffusion_turb]
    type = INSFVMixingLengthScalarDiffusion
    variable = c5
  []
  [c6_diffusion_turb]
    type = INSFVMixingLengthScalarDiffusion
    variable = c6
  []
  [c1_src]
    type = FVCoupledForce
    variable = c1
    v = fission_source
    coef = ${beta1}
  []
  [c2_src]
    type = FVCoupledForce
    variable = c2
    v = fission_source
    coef = ${beta2}
  []
  [c3_src]
    type = FVCoupledForce
    variable = c3
    v = fission_source
    coef = ${beta3}
  []
  [c4_src]
    type = FVCoupledForce
    variable = c4
    v = fission_source
    coef = ${beta4}
  []
  [c5_src]
    type = FVCoupledForce
    variable = c5
    v = fission_source
    coef = ${beta5}
  []
  [c6_src]
    type = FVCoupledForce
    variable = c6
    v = fission_source
    coef = ${beta6}
  []
  [c1_decay]
    type = FVReaction
    variable = c1
    rate = ${lambda1}
  []
  [c2_decay]
    type = FVReaction
    variable = c2
    rate = ${lambda2}
  []
  [c3_decay]
    type = FVReaction
    variable = c3
    rate = ${lambda3}
  []
  [c4_decay]
    type = FVReaction
    variable = c4
    rate = ${lambda4}
  []
  [c5_decay]
    type = FVReaction
    variable = c5
    rate = ${lambda5}
  []
  [c6_decay]
    type = FVReaction
    variable = c6
    rate = ${lambda6}
  []
[]

[FVBCs]
  [no-slip-u]
    type = INSFVNoSlipWallBC
    boundary = 'wall-reactor wall-pipe wall-pump wall-reactor-full'
    variable = superficial_vel_x
    function = 0
  []
  [no-slip-v]
    type = INSFVNoSlipWallBC
    boundary = 'wall-reactor wall-pipe wall-pump wall-reactor-full'
    variable = superficial_vel_y
    function = 0
  []
  [no-slip-w]
    type = INSFVNoSlipWallBC
    boundary = 'wall-reactor wall-pipe wall-pump wall-reactor-full'
    variable = superficial_vel_z
    function = 0
  []
[]

################################################################################
# MATERIALS
################################################################################

[Materials]
  [porous_mat]
    type = ADGenericFunctorMaterial #defines mu artificially for numerical convergence
    prop_names = 'porosity' #it converges to the real mu eventually.
    prop_values = '${porosity}'
    block = 'reactor pipe pump mixing-plate reflector'
  []
  [mixing_length_mat]
    type = ADPiecewiseByBlockFunctorMaterial
    prop_name = 'mixing_length'
    subdomain_to_prop_value = 'reactor      ${mixing_length_reactor_callibrated}
                               mixing-plate ${mixing_length_reactor_callibrated}
                               pipe         ${mixing_length_pipe_callibrated}
                               pump         ${mixing_length_pipe_callibrated}'
  []
[]

################################################################################
# UTILITY FUNCTIONS
################################################################################

[Functions]
  [dts]
    type = PiecewiseLinear
    x = '0   1e7 ${fparse 1e7+0.1} ${fparse 1e7+0.2}'
    y = '0.1 0.1 1e6               1e6'
  []
[]

[Debug]
  show_var_residual_norms = true
[]

################################################################################
# SIMULATION OUTPUTS
################################################################################

[Outputs]
  csv = true
  hide = 'dt_limit'
  [restart]
    type = Exodus
    execute_on = 'timestep_end final'
  []
  # Reduce base output
  print_linear_converged_reason = false
  print_linear_residuals = false
  print_nonlinear_converged_reason = false
[]

[Postprocessors]
  [dt_limit]
    type = Receiver
    default = 1
  []
[]
