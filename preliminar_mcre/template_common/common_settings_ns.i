 ################################################################################
## Molten Salt Fast Reactor - Euratom EVOL + Rosatom MARS Design              ##
## Pronghorn input file to initialize velocity fields                         ##
## This runs a slow relaxation to steady state while ramping down the fluid   ##
## viscosity.                                                                 ##
################################################################################

# Material properties fuel

# rho = 3279. # density [kg / m^3]  (@1000K)
# mu = 0.005926  # viscosity [Pa s]
# k = 0.38
#cp = 640.

rho = 'rho'
mu = 'mu'
k = 'k'
cp = 'cp'
mu_true = 0.005926  # viscosity [Pa s]


# rho = ${FluidProperties/fluid_properties/rho} # density [kg / m^3]  (@1000K)
# mu = ${FluidProperties/fluid_properties/mu} # viscosity [Pa s]
# k = ${FluidProperties/fluid_properties/k}
# cp = ${FluidProperties/fluid_properties/cp}

# Material properties reflector
k_ref = 30.
cp_ref = 880.
rho_ref = 3580.


#power = 25e3
power = 1e6

#alpha_b = '${fparse 1.0/rho}'

# Mass flow rate tuning
friction = 11.0 # [kg / m^4]
pump_force = '${fparse 0.06*0.25*4.0e6}' # [N / m^3]
porosity = 1.0
#T_hx = 592
Reactor_Area = '${fparse 3.14159*0.2575*0.2575}'


#Pr = '${fparse mu*cp/k}'

# Numerical scheme parameters
advected_interp_method = 'upwind'
velocity_interp_method = 'rc'

mixing_length_pipe_callibrated = '${fparse 0.07 * 0.1 * 0.06}'
mixing_length_reactor_callibrated = '${fparse 0.07 * 0.1 * 2}'

[GlobalParams]
  rhie_chow_user_object = 'pins_rhie_chow_interpolator'

  two_term_boundary_expansion = true

  advected_interp_method = ${advected_interp_method}
  velocity_interp_method = ${velocity_interp_method}

  u = superficial_vel_x
  v = superficial_vel_y
  w = superficial_vel_z

  pressure = pressure
  porosity = porosity

  rho = rho
  mu = mu
  speed = speed

  mixing_length = 'mixing_length'

  #gravity = '-9.81 0 0'
[]

################################################################################
# EQUATIONS: VARIABLES, KERNELS & BCS
################################################################################

[UserObjects]
  [pins_rhie_chow_interpolator]
    type = PINSFVRhieChowInterpolator
    block = 'reactor pipe pump mixing-plate'
  []
  [pin_pressure]
    type = NSFVPressurePin
    pin_type = 'point-value'
    point = '0 0 0'
    phi0 = 'phi0-pp'
    block = 'reactor pipe pump mixing-plate'
  []
[]

[Variables]
  [pressure]
    type = INSFVPressureVariable
    block = 'reactor pipe pump mixing-plate'
  []
  [superficial_vel_x]
    type = PINSFVSuperficialVelocityVariable
    block = 'reactor pipe pump mixing-plate'
  []
  [superficial_vel_y]
    type = PINSFVSuperficialVelocityVariable
    block = 'reactor pipe pump mixing-plate'
  []
  [superficial_vel_z]
    type = PINSFVSuperficialVelocityVariable
    block = 'reactor pipe pump mixing-plate'
  []
  [T]
    type = INSFVEnergyVariable
    block = 'reactor pipe pump mixing-plate'
  []
  [T_ref]
    type = INSFVEnergyVariable
    block = 'reflector'
    scaling = 0.001
  []
[]

[FVKernels]
  [mass]
    type = PINSFVMassAdvection
    variable = pressure
    advected_interp_method = 'skewness-corrected'
    velocity_interp_method = 'rc'
    rho = ${rho}
  []

  [u_time]
    type = PINSFVMomentumTimeDerivative
    variable = superficial_vel_x
    momentum_component = 'x'
    rho = ${rho}
  []
  [u_advection]
    type = PINSFVMomentumAdvection
    variable = superficial_vel_x
    momentum_component = 'x'
    rho = ${rho}
  []
  [u_viscosity]
    type = PINSFVMomentumDiffusion
    variable = superficial_vel_x
    momentum_component = 'x'
  []
  [u_viscosity_rans_pipe]
    type = INSFVMixingLengthReynoldsStress
    variable = superficial_vel_x
    momentum_component = 'x'
    rho = ${rho}
  []
  [u_pressure]
    type = PINSFVMomentumPressure
    variable = superficial_vel_x
    momentum_component = 'x'
    pressure = pressure
  []
  [u_friction_pump]
    type = PINSFVMomentumFriction
    variable = superficial_vel_x
    momentum_component = 'x'
    Darcy_name = 'DFC'
    rho = ${rho}
    Forchheimer_name = 'FFC'
    block = 'pump'
  []
  [u_friction_pump_correction]
    type = PINSFVMomentumFrictionCorrection
    variable = superficial_vel_x
    momentum_component = 'x'
    rho = ${rho}
    Darcy_name = 'DFC'
    Forchheimer_name = 'FFC'
    consistent_scaling = 100.
    block = 'pump'
  []
  [u_friction_mixing_plate]
    type = PINSFVMomentumFriction
    variable = superficial_vel_x
    momentum_component = 'x'
    Darcy_name = 'DFC_plate'
    Forchheimer_name = 'FFC_plate'
    block = 'mixing-plate'
  []
  [u_friction_mixing_plate_correction]
    type = PINSFVMomentumFrictionCorrection
    variable = superficial_vel_x
    momentum_component = 'x'
    rho = ${rho}
    Darcy_name = 'DFC_plate'
    Forchheimer_name = 'FFC_plate'
    consistent_scaling = 100.
    block = 'mixing-plate'
  []

  [v_time]
    type = PINSFVMomentumTimeDerivative
    variable = superficial_vel_y
    momentum_component = 'y'
    rho = ${rho}
  []
  [v_advection]
    type = PINSFVMomentumAdvection
    variable = superficial_vel_y
    momentum_component = 'y'
    rho = ${rho}
  []
  [v_viscosity]
    type = PINSFVMomentumDiffusion
    variable = superficial_vel_y
    momentum_component = 'y'
  []
  [v_viscosity_rans_pipe]
    type = INSFVMixingLengthReynoldsStress
    variable = superficial_vel_y
    momentum_component = 'y'
    rho = ${rho}
    mixing_length = ${mixing_length_pipe_callibrated}
    block = 'pipe pump'
  []
  [v_viscosity_rans]
    type = INSFVMixingLengthReynoldsStress
    variable = superficial_vel_y
    momentum_component = 'y'
    rho = ${rho}
  []
  [v_pressure]
    type = PINSFVMomentumPressure
    variable = superficial_vel_y
    momentum_component = 'y'
    pressure = pressure
  []
  [pump]
    type = INSFVBodyForce
    variable = superficial_vel_y
    momentum_component = 'y'
    functor = ${pump_force}
    block = 'pump'
  []
  [v_friction_pump]
    type = PINSFVMomentumFriction
    variable = superficial_vel_y
    momentum_component = 'y'
    Darcy_name = 'DFC'
    Forchheimer_name = 'FFC'
    block = 'pump'
  []
  [v_friction_pump_correction]
    type = PINSFVMomentumFrictionCorrection
    variable = superficial_vel_y
    momentum_component = 'y'
    rho = ${rho}
    Darcy_name = 'DFC'
    Forchheimer_name = 'FFC'
    consistent_scaling = 100.
    block = 'pump'
  []
  [v_friction_mixing_plate]
    type = PINSFVMomentumFriction
    variable = superficial_vel_y
    momentum_component = 'y'
    Darcy_name = 'DFC_plate'
    Forchheimer_name = 'FFC_plate'
    block = 'mixing-plate'
  []
  [v_friction_mixing_plate_correction]
    type = PINSFVMomentumFrictionCorrection
    variable = superficial_vel_y
    momentum_component = 'y'
    rho = ${rho}
    Darcy_name = 'DFC_plate'
    Forchheimer_name = 'FFC_plate'
    consistent_scaling = 100.
    block = 'mixing-plate'
  []

  [w_time]
    type = PINSFVMomentumTimeDerivative
    variable = superficial_vel_z
    momentum_component = 'z'
    rho = ${rho}

  []
  [w_advection]
    type = PINSFVMomentumAdvection
    variable = superficial_vel_z
    momentum_component = 'z'
    rho = ${rho}
  []
  [w_viscosity]
    type = PINSFVMomentumDiffusion
    variable = superficial_vel_z
    momentum_component = 'z'
  []
  [w_viscosity_rans]
    type = INSFVMixingLengthReynoldsStress
    variable = superficial_vel_z
    momentum_component = 'z'
    rho = ${rho}
  []
  [w_pressure]
    type = PINSFVMomentumPressure
    variable = superficial_vel_z
    momentum_component = 'z'
    pressure = pressure
  []
  [w_friction_pump]
    type = PINSFVMomentumFriction
    variable = superficial_vel_z
    momentum_component = 'z'
    Darcy_name = 'DFC'
    Forchheimer_name = 'FFC'
    block = 'pump'
  []
  [w_friction_pump_correction]
    type = PINSFVMomentumFrictionCorrection
    variable = superficial_vel_z
    momentum_component = 'z'
    rho = ${rho}
    Darcy_name = 'DFC'
    Forchheimer_name = 'FFC'
    consistent_scaling = 100.
    block = 'pump'
  []
  [w_friction_mixing_plate]
    type = PINSFVMomentumFriction
    variable = superficial_vel_z
    momentum_component = 'z'
    Darcy_name = 'DFC_plate'
    Forchheimer_name = 'FFC_plate'
    block = 'mixing-plate'
  []
  [w_friction_mixing_plate_correction]
    type = PINSFVMomentumFrictionCorrection
    variable = superficial_vel_z
    momentum_component = 'z'
    rho = ${rho}
    Darcy_name = 'DFC_plate'
    Forchheimer_name = 'FFC_plate'
    consistent_scaling = 100.
    block = 'mixing-plate'
  []

  ####### FUEL ENERGY EQUATION #######

  [heat_time]
    type = PINSFVEnergyTimeDerivative
    variable = T
    is_solid = false
    cp = cp
  []
  [heat_advection]
    type = PINSFVEnergyAdvection
    variable = T
  []
  [heat_diffusion]
    type = PINSFVEnergyDiffusion
    variable = T
    k = ${k}
  []
  [heat_turbulent_diffusion]
    type = WCNSFVMixingLengthEnergyDiffusion
    variable = T
    schmidt_number = 0.9
    cp = cp
  []

  ####### REFLECTOR ENERGY EQUATION #######

  [heat_time_ref]
    type = INSFVEnergyTimeDerivative
    variable = T_ref
    rho = ${rho_ref}
    dh_dt = dh_solid_dt
  []
  [heat_diffusion_ref]
    type = FVDiffusion
    variable = T_ref
    coeff = ${k_ref}
  []
[]

[FVInterfaceKernels]
  [convection]
    type = FVConvectionCorrelationInterface
    variable1 = T
    variable2 = T_ref
    boundary = 'wall-reactor-reflector'
    h = htc
    T_solid = T_ref
    T_fluid = T
    subdomain1 = reactor
    subdomain2 = reflector
    wall_cell_is_bulk = true
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
  [heat-losses-outshield]
    type = FVFunctorConvectiveHeatFluxBC
    variable = T
    T_bulk = T
    T_solid = 300.
    is_solid = false
    heat_transfer_coefficient = 3. #50.
    boundary = 'heat-loss-section-outshield'
  []
  [heated-outshield-pipe]
    type = FVFunctorNeumannBC
    variable = T
    boundary = 'heat-loss-section-outshield'
    functor = heat-input-pipe
  []
  [heat-losses-reflector]
    type = FVFunctorConvectiveHeatFluxBC
    variable = T_ref
    T_bulk = 350.
    T_solid = T_ref
    is_solid = true
    heat_transfer_coefficient = htc_rad_ref
    boundary = 'wall-reflector'
  []
  [heated-reflector-walls]
    type = FVFunctorNeumannBC
    variable = T_ref
    boundary = 'heated-reflector-walls'
    functor = heat-input-ref
  []
  [heated-inshield-pipe]
    type = FVNeumannBC
    variable = T
    boundary = 'heat-loss-section-inshield'
    value = 0.0 #500.
  []
[]

[AuxVariables]
  [h_DeltaT]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
  []
  [h_DeltaT_rad_aux]
    type = MooseVariableFVReal
    block = 'reflector'
  []
  [h_DeltaT_rad]
    type = MooseVariableFVReal
    block = 'reflector'
  []
  [a_u_var]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
  []
  [a_v_var]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
  []
  [a_w_var]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
  []
  [power_density]
    type = MooseVariableFVReal
    # [InitialCondition]
    #   type = FunctionIC
    #   function = 'power-density-func'
    # []
  []
  [fission_source]
    type = MooseVariableFVReal
    # Fission source is re-initalized by a transfer from neutronics
    # [InitialCondition]
    #   type = FunctionIC
    #   function = 'power-density-func'
    # []
  []
  [c1]
    type = MooseVariableFVReal
    initial_from_file_var = c1
    block = 'reactor pipe pump mixing-plate'
  []
  [c2]
    type = MooseVariableFVReal
    initial_from_file_var = c2
    block = 'reactor pipe pump mixing-plate'
  []
  [c3]
    type = MooseVariableFVReal
    initial_from_file_var = c3
    block = 'reactor pipe pump mixing-plate'
  []
  [c4]
    type = MooseVariableFVReal
    initial_from_file_var = c4
    block = 'reactor pipe pump mixing-plate'
  []
  [c5]
    type = MooseVariableFVReal
    initial_from_file_var = c5
    block = 'reactor pipe pump mixing-plate'
  []
  [c6]
    type = MooseVariableFVReal
    initial_from_file_var = c6
    block = 'reactor pipe pump mixing-plate'
  []
  [pressure_initial]
    type = MooseVariableFVReal
    initial_condition = 101325
  []
[]

[AuxKernels]
  [h_DeltaT_out]
    type = ParsedAux
    variable = h_DeltaT
    coupled_variables = 'T'
    expression = '3.*(T-300.)'
  []
  [h_DeltaT_rad_out_pre]
    type = ParsedAux
    variable = h_DeltaT_rad_aux
    coupled_variables = 'T_ref'
    expression = 'T_ref-350.'
  []
  [h_DeltaT_rad_out]
    type = FunctorAux
    functor = 'h_DeltaT_rad_aux'
    variable = h_DeltaT_rad
    factor = 'htc_rad_ref'
  []
  [comp_a_u]
    type = FunctorAux
    functor = 'ax'
    variable = 'a_u_var'
    block = 'reactor pipe pump mixing-plate'
    execute_on = 'timestep_end'
  []
  [comp_a_v]
    type = FunctorAux
    functor = 'ay'
    variable = 'a_v_var'
    block = 'reactor pipe pump mixing-plate'
    execute_on = 'timestep_end'
  []
  [comp_a_w]
    type = FunctorAux
    functor = 'az'
    variable = 'a_w_var'
    block = 'reactor pipe pump mixing-plate'
    execute_on = 'timestep_end'
  []
[]

################################################################################
# MATERIALS
################################################################################

[Functions]
  [heat-input-ref]
    type = ParsedFunction
    expression = '0.'
  []
  [heat-input-pipe]
    type = ParsedFunction
    expression = '0.'
  []
  [Re_reactor]
    type = ParsedFunction
    expression = 'flow_hx_bot/Reactor_Area * (2*0.2575) / mu '
    symbol_names = 'flow_hx_bot mu Reactor_Area'
    #symbol_values = 'flow_hx_bot ${mu} ${Reactor_Area}'
    symbol_values = '25.2 viscosity_avg ${Reactor_Area}'
  []
  [Pr_reactor]
    type = ParsedFunction
    expression = 'mu * cp / k'
    symbol_names = 'mu cp k'
    symbol_values = 'viscosity_avg cp_avg k_avg'
  []
  [htc]
    type = ParsedFunction
    expression = '600 *k/(2*0.2575)* 0.023 * Re_reactor^0.8 * Pr^0.3'
    symbol_names = 'k Re_reactor Pr'
    symbol_values = 'k_avg Re_reactor Pr_reactor'
  []
  [htc_rad_ref]
    type = ParsedFunction
    expression = '(T_ref_wall^2+350.^2)*(T_ref_wall+350.)*5.67e-8 / (1/0.18+1/0.35-1.)' #https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node136.html #e_mgo=0.18 #e_steel=0.35
    symbol_names = 'T_ref_wall'
    symbol_values = 'T_ref_wall'
  []
  [ad_rampdown_mu_func]
    type = ADParsedFunction
    expression = mu*(0.1*exp(-3*t)+1)
    symbol_names = 'mu'
    symbol_values = ${mu_true}
  []
  [cosine_guess]
    type = ParsedFunction
    expression = 'max(0, cos(x*pi/2/1.5))*max(0, cos(y*pi/2/1.5))*max(0, cos(z*pi/2/1.5))'
  []
  [power-density-func]
    type = ParsedFunction
    expression = '${power}/0.21757 * cg'
    symbol_names = 'cg'
    symbol_values = 'cosine_guess'
  []
  [dts]
    type = PiecewiseLinear
    x = '0   1e7 ${fparse 1e7+0.05} ${fparse 1e7+0.2}'
    y = '0.1 0.1 1e7                1e7'
  []
[]

[FunctorMaterials]
  [test]
    type = GeneralFunctorFluidProps
    T_fluid = T 
    characteristic_length = 1
    fp = fluid_properties
    outputs = 'all'
    output_properties = 'rho mu cp k'
    force_define_density = true
    block = 'reactor pipe pump mixing-plate'
  []
[]

# [Materials]
#   [to_vars]
#     type = FluidPropertiesMaterialPT
#     fp = fluid_properties
#     pressure = pressure_initial
#     temperature = T
#     outputs = 'all'
#     output_properties = 'density cp k viscosity'
#     block ='reactor pipe pump mixing-plate'
#     compute_entropy = false
#     compute_sound_speed = false
#   []
# []

[FluidProperties]
  [fluid_properties]
    type = TemperaturePressureFunctionFluidProperties
    rho = rho_fluid
    cp = cp_fluid
    k = k_fluid
    mu = mu_fluid
  []
[]


[Functions]
  [rho_fluid]
    type = ParsedFunction
    expression = '(4212.6-1.0686 * x)'
  []
  [cp_fluid]
    type = ParsedFunction
    symbol_names = 'perturb'
    symbol_values = 'cp_perturb_sub'
    expression = '(8900.439 + (-1.377936e1 * x) + (6.400369e-3 * x^2) + (-8.443758e8 / (x^2))) * perturb'
  []
  [k_fluid]
    type = ParsedFunction
    symbol_names = 'perturb'
    symbol_values = 'k_perturb_sub'
    expression = '(5.6820) + (-8.7832e-3 * x) + (4.0967e-6 * x^2 ) + (-5.7642e5/(x^2)) * perturb'
  []
  [mu_fluid]
    type = ParsedFunction
    symbol_names = 'perturb'
    symbol_values = 'mu_perturb_sub'
    expression = '(1.505e-4 * exp (2.666e4/(8.314*x)))*perturb'
  []
[]


[Postprocessors]
  [average_T]
    type = ElementAverageValue
    variable = T
    block = 'pipe reactor pump mixing-plate'
    execute_on = 'initial timestep_end'
  []
  [rho_avg]
    type = ElementAverageValue
    variable = 'rho_out'
    block ='pipe pump reactor mixing-plate'
    execute_on = 'initial timestep_end'
  []
  [cp_avg]
    type = ElementAverageValue
    variable = 'cp_out'
    block ='pipe pump reactor mixing-plate'
    execute_on = 'initial timestep_end'
  []
  [k_avg]
    type = ElementAverageValue
    variable = 'k_out'
    block ='pipe pump reactor mixing-plate'
    execute_on = 'initial timestep_end'
  []
  [viscosity_avg]
    type = ElementAverageValue
    variable = 'mu_out'
    block ='pipe pump reactor mixing-plate'
    execute_on = 'initial timestep_end'
  []
  [rho_perturb_sub]
    type = Receiver
    default = 1
    outputs = none
  []
  [cp_perturb_sub]
    type = Receiver
    default = 1
    outputs = none
  []
  [k_perturb_sub]
    type = Receiver
    default = 1
    outputs = none
  []
  [mu_perturb_sub]
    type = Receiver
    default = 1
    outputs = none
  []
[]



[FunctorMaterials]

  [constant_porosity]
    type = ADGenericFunctorMaterial #defines mu artificially for numerical convergence
    prop_names = 'porosity' #it converges to the real mu eventually.
    prop_values = '${porosity}'
  []
  # [mu_spatial]
  #   type = ADPiecewiseByBlockFunctorMaterial
  #   prop_name = 'mu'
  #   subdomain_to_prop_value = 'pipe         ad_rampdown_mu_func
  #                              pump         ad_rampdown_mu_func
  #                              mixing-plate ad_rampdown_mu_func
  #                              reactor      ad_rampdown_mu_func'
  # []
  [friction_material_pump]
    type = ADGenericVectorFunctorMaterial #defines mu artificially for numerical convergence
    prop_names = 'DFC FFC' #it converges to the real mu eventually.
    prop_values = '${fparse 1*friction} ${fparse 1*friction} ${fparse 1*friction}
                   ${fparse 1*friction} ${fparse 1*friction} ${fparse 1*friction}'
  []
  [friction_material_mixing_plate]
    type = ADGenericVectorFunctorMaterial #defines mu artificially for numerical convergence
    prop_names = 'DFC_plate FFC_plate' #it converges to the real mu eventually.
    prop_values = '${fparse 1*friction} ${fparse 1*friction} ${fparse 1*friction}
                   ${fparse 1*friction} ${fparse 1*friction} ${fparse 1*friction}'
  []
  [ins_fv]
    type = INSFVEnthalpyFunctorMaterial
    rho = ${rho}
    cp = cp
    temperature = 'T'
    block = 'reactor pipe pump mixing-plate'
  []
  [mixing_length_mat]
    type = ADPiecewiseByBlockFunctorMaterial
    prop_name = 'mixing_length'
    subdomain_to_prop_value = 'reactor      ${mixing_length_reactor_callibrated}
                               mixing-plate ${mixing_length_reactor_callibrated}
                               pipe         ${mixing_length_pipe_callibrated}
                               pump         ${mixing_length_pipe_callibrated}'
  []
  [ins_fv_solid]
    type = INSFVEnthalpyFunctorMaterial
    temperature = 'T_ref'
    rho = ${rho_ref}
    cp = ${cp_ref}
    h = h_solid
    rho_h = rho_h_solid
    block = 'reflector'
  []
  [speed]
    type = PINSFVSpeedFunctorMaterial
    superficial_vel_x = superficial_vel_x
    superficial_vel_y = superficial_vel_y
    superficial_vel_z = superficial_vel_z
  []
[]

[Debug]
  show_var_residual_norms = true
[]

################################################################################
# SIMULATION OUTPUTS
################################################################################

[Outputs]
  [out]
    type = CSV
    execute_on = 'initial timestep_end final'
  []
  [restart]
    type = Exodus
    execute_on = 'timestep_end final'
  []
  # Reduce base output
  print_linear_converged_reason = false
  print_linear_residuals = false
  print_nonlinear_converged_reason = false
  hide = 'dt_limit'
[]

[Postprocessors]
  [phi0-pp]
    type = Receiver
    default = 0.0
    outputs = out
  []

  [reactor-power]
    type = ElementIntegralVariablePostprocessor
    variable = power_density
    execute_on = 'initial timestep_begin timestep_end transfer'
    outputs = all
    block = 'reactor pump pipe'
  []
  [reactor-power-func]
    type = ElementIntegralFunctorPostprocessor
    functor = 'power-density-func'
    block = 'reactor'
    outputs = out
  []

  [flow]
    type = VolumetricFlowRate
    boundary = 'reactor_bot'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    vel_z = superficial_vel_z
    advected_quantity = ${rho}
    #outputs = out
  []
  [max_v]
    type = ElementExtremeValue
    variable = superficial_vel_x
    value_type = max
    block = 'reactor pump pipe'
    outputs = out
  []
  [pdrop-mixing-plate]
    type = PressureDrop
    pressure = pressure
    upstream_boundary = 'reactor_bot'
    downstream_boundary = 'mixing-plate-downstream'
    boundary = 'mixing-plate-downstream reactor_bot'
    outputs = out
  []

  [heat-loss-pipe]
    type = SideIntegralVariablePostprocessor
    boundary = 'heat-loss-section-outshield'
    variable = h_DeltaT
    outputs = out
  []
  [heat-loss-ref]
    type = SideIntegralVariablePostprocessor
    boundary = 'wall-reflector'
    variable = h_DeltaT_rad
    outputs = out
  []
  [heat-input-pipe-pp]
    type = FunctionSideIntegral
    boundary = 'heat-loss-section-outshield'
    function = heat-input-pipe
    outputs = out
  []
  [heat-input-ref-pp]
    type = FunctionSideIntegral
    boundary = 'heated-reflector-walls'
    function = heat-input-ref
    outputs = out
  []

  [T_reac_bot]
    type = SideAverageValue
    boundary = 'reactor_bot'
    variable = T
  []
  [T_reac_top]
    type = SideAverageValue
    boundary = 'reactor_top'
    variable = T
  []
  [T_reac_ave]
    type = ElementAverageValue
    variable = T
    block = 'reactor'
    outputs = out
  []
  [T_reac_max]
    type = ElementExtremeValue
    variable = T
    block = 'reactor'
    outputs = out
  []

  [T_pipe_ave]
    type = ElementAverageValue
    variable = T
    block = 'pipe'
    outputs = out
  []

  [T_ref_wall]
    type = SideAverageValue
    boundary = 'wall-reflector'
    variable = T_ref

  []
  [T_ref_ave]
    type = ElementAverageValue
    variable = T_ref
    block = 'reflector'
  []
  [T_ref_max]
    type = ElementExtremeValue
    variable = T_ref
    block = 'reflector'
    outputs = out
  []
  # [rho_T_average]
  #   type = FunctionElementAverage
  #   function = rho_T
  #   block = 'reactor'
  # []
  # [cp_T_average]
  #   type = FunctionElementAverage
  #   function = cp_T
  #   block = 'reactor'
  # []
  # [k_T_average]
  #   type = FunctionElementAverage
  #   function = k_T
  #   block = 'reactor'
  # []
  # [mu_T_average]
  #   type = FunctionElementAverage
  #   function = mu_T
  #   block = 'reactor'
  # [] 
  # [tfuel]
  #   type = Receiver
  #   default = 900
  # []
  [c1_outlet]
    type = VolumetricFlowRate
    boundary = 'reactor_top'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    vel_z = superficial_vel_z
    advected_quantity = c1
    outputs = out
  []
  [c2_outlet]
    type = VolumetricFlowRate
    boundary = 'reactor_top'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    vel_z = superficial_vel_z
    advected_quantity = c2
    outputs = out
  []
  [c3_outlet]
    type = VolumetricFlowRate
    boundary = 'reactor_top'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    vel_z = superficial_vel_z
    advected_quantity = c3
    outputs = out
  []
  [c4_outlet]
    type = VolumetricFlowRate
    boundary = 'reactor_top'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    vel_z = superficial_vel_z
    advected_quantity = c4
    outputs = out
  []
  [c5_outlet]
    type = VolumetricFlowRate
    boundary = 'reactor_top'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    vel_z = superficial_vel_z
    advected_quantity = c5
    outputs = out
  []
  [c6_outlet]
    type = VolumetricFlowRate
    boundary = 'reactor_top'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    vel_z = superficial_vel_z
    advected_quantity = c6
    outputs = out
  []

  [c1_inlet]
    type = VolumetricFlowRate
    boundary = 'reactor_bot'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    vel_z = superficial_vel_z
    advected_quantity = c1
    outputs = out
  []
  [c2_inlet]
    type = VolumetricFlowRate
    boundary = 'reactor_bot'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    vel_z = superficial_vel_z
    advected_quantity = c2
    outputs = out
  []
  [c3_inlet]
    type = VolumetricFlowRate
    boundary = 'reactor_bot'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    vel_z = superficial_vel_z
    advected_quantity = c3
    outputs = out
  []
  [c4_inlet]
    type = VolumetricFlowRate
    boundary = 'reactor_bot'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    vel_z = superficial_vel_z
    advected_quantity = c4
    outputs = out
  []
  [c5_inlet]
    type = VolumetricFlowRate
    boundary = 'reactor_bot'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    vel_z = superficial_vel_z
    advected_quantity = c5
    outputs = out
  []
  [c6_inlet]
    type = VolumetricFlowRate
    boundary = 'reactor_bot'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    vel_z = superficial_vel_z
    advected_quantity = c6
    outputs = out
  []

  [dt_limit]
    type = Receiver
    default = 1e6
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
#                       c1_inlet c2_inlet c3_inlet c4_inlet c5_inlet c6_inlet
#                       c1_outlet c2_outlet c3_outlet c4_outlet c5_outlet c6_outlet'
#   []
# []

################################################################################
# MULTIAPPS and TRANSFERS for precursors transport
################################################################################

[MultiApps]
  active = ''
[]

[Transfers]
  active = ''
  [power_density]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = prec_transport
    source_variable = power_density
    variable = power_density
  []
  [fission_source]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = prec_transport
    source_variable = fission_source
    variable = fission_source
  []
  [u_x]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = prec_transport
    source_variable = superficial_vel_x
    variable = superficial_vel_x
  []
  [u_y]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = prec_transport
    source_variable = superficial_vel_y
    variable = superficial_vel_y
  []
  [u_z]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = prec_transport
    source_variable = superficial_vel_z
    variable = superficial_vel_z
  []
  [a_u]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = prec_transport
    source_variable = a_u_var
    variable = a_u_var
  []
  [a_v]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = prec_transport
    source_variable = a_v_var
    variable = a_v_var
  []
  [a_w]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = prec_transport
    source_variable = a_w_var
    variable = a_w_var
  []

  [c1]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = prec_transport
    source_variable = 'c1'
    variable = 'c1'
  []
  [c2]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = prec_transport
    source_variable = 'c2'
    variable = 'c2'
  []
  [c3]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = prec_transport
    source_variable = 'c3'
    variable = 'c3'
  []
  [c4]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = prec_transport
    source_variable = 'c4'
    variable = 'c4'
  []
  [c5]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = prec_transport
    source_variable = 'c5'
    variable = 'c5'
  []
  [c6]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = prec_transport
    source_variable = 'c6'
    variable = 'c6'
  []
[]
