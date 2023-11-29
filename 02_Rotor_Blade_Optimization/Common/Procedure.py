# Procedure.py
# 
# Created: Jun 2022, E. Botero
# Modified: 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import MARC
from MARC.Core import Units
from MARC.Analyses.Process import Process
from MARC.Optimization.write_optimization_outputs import write_optimization_outputs
import MARC.Components.Energy.Converters.Propeller as prop_class

import numpy as np
from jax import jacfwd as jacobian
import jax.numpy as jnp
import jax.numpy.linalg as LA
from jax import jit
from MARC.Analyses.Mission.Segments.Segment   import Segment 
from MARC.Methods.Noise.Fidelity_Zero.Propeller.propeller_mid_fidelity import propeller_mid_fidelity

# ----------------------------------------------------------------------        
#   Setup
# ----------------------------------------------------------------------   

def setup(rs):
    
    # Setup variables in the problem formulation
    
    
    # ------------------------------------------------------------------
    #   Analysis Procedure
    # ------------------------------------------------------------------ 
    
    # size the base config
    procedure = Process()
    
    # If the discretization is the new style we'll need to update it
    procedure.set_rotor = set_rotor
    
    # Run the rotor in hover
    procedure.hover = run_rotor_hover
    
    # Run the rotor in hover
    procedure.OEI = run_rotor_OEI    
    
    # If we have alpha == 1 we don't need noise
    if rs.alpha !=1:
        procedure.hover_noise = run_hover_noise
        
    # If we have a prop rotor
    if rs.prop_rotor:
        procedure.run_rotor_forward = run_rotor_forward
        
        # If we have alpha == 1 we don't need noise
        if rs.alpha !=1:
            procedure.forward_noise = run_forward_noise
        

    # post process the results
    procedure.post_process = post_process
        
    return procedure


# ----------------------------------------------------------------------        
#   Set Rotor
# ----------------------------------------------------------------------   
@jit
def set_rotor(nexus):
    
    # unpacks
    rs             = nexus.rotor_settings
    vehicle_hover  = nexus.vehicle_configurations.hover 
    vehicle_cruise = nexus.vehicle_configurations.cruise 
    vehicle_oei    = nexus.vehicle_configurations.oei 
    
    rotor_h   = vehicle_hover.networks.Battery_Electric_Rotor.rotors.prop_rotor
    rotor_c   = vehicle_cruise.networks.Battery_Electric_Rotor.rotors.prop_rotor
    rotor_oei = vehicle_oei.networks.Battery_Electric_Rotor.rotors.prop_rotor
        
    if not rs.traditional_format:     
            
        airfoil_geometry_data = rotor_h.airfoil_geometry_data 
            
        # Update geometry of blade 
        c       = updated_blade_geometry(rotor_h.radius_distribution/rotor_h.tip_radius ,rotor_h.chord_r,rotor_h.chord_p,rotor_h.chord_q,rotor_h.chord_t)     
        beta    = updated_blade_geometry(rotor_h.radius_distribution/rotor_h.tip_radius ,rotor_h.twist_r,rotor_h.twist_p,rotor_h.twist_q,rotor_h.twist_t)   
        
    
        # compute max thickness distribution   
        t_max  = airfoil_geometry_data.max_thickness*c   
         
        rotor_h.chord_distribution          = c
        rotor_h.twist_distribution          = beta  + rs.root_incidence
        rotor_h.mid_chord_alignment         = c/4. - c[0]/4.
        rotor_h.max_thickness_distribution  = t_max 
        
        rotor_oei.chord_distribution         = rotor_h.chord_distribution
        rotor_oei.twist_distribution         = rotor_h.twist_distribution
        rotor_oei.mid_chord_alignment        = rotor_h.mid_chord_alignment  
        rotor_oei.max_thickness_distribution = rotor_h.max_thickness_distribution
        rotor_oei.thickness_to_chord         = rotor_h.thickness_to_chord    
          
        rotor_c.chord_distribution         = rotor_h.chord_distribution
        rotor_c.twist_distribution         = rotor_h.twist_distribution
        rotor_c.mid_chord_alignment        = rotor_h.mid_chord_alignment  
        rotor_c.max_thickness_distribution = rotor_h.max_thickness_distribution
        rotor_c.thickness_to_chord         = rotor_h.thickness_to_chord  
        
    # Pack the twists and chords
    rotor_h.twist_distribution   = rotor_h.twist_distribution + rs.root_incidence
    rotor_c.twist_distribution   = rotor_h.twist_distribution
    rotor_c.chord_distribution   = rotor_h.chord_distribution
    rotor_oei.twist_distribution = rotor_h.twist_distribution
    rotor_oei.chord_distribution = rotor_h.chord_distribution 

    vehicle_hover.store_diff()  
    vehicle_oei.store_diff()     
    vehicle_cruise.store_diff()  
    
    return nexus
# ----------------------------------------------------------------------
#   Update blade geometry 
# ---------------------------------------------------------------------- 
@jit
def updated_blade_geometry(chi,c_r,p,q,c_t):
    """ Computes planform function of twist and chord distributron using hyperparameters  
          
          Inputs:  
             chi - prop-rotor radius distribution [None]
             c_r - hyperparameter no. 1           [None]
             p   - hyperparameter no. 2           [None]
             q   - hyperparameter no. 3           [None]
             c_t - hyperparameter no. 4           [None]
                   
          Outputs:       
             x_lin  - function distribution       [None]
              
          Assumptions: 
             N/A 
        
          Source:
              Traub, Lance W., et al. "Effect of taper ratio at low reynolds number."
              Journal of Aircraft 52.3 (2015): 734-747.
              
    """           

    n       = jnp.linspace(len(chi)-1,0,len(chi))          
    theta_n = n*(jnp.pi/2)/len(chi)              
    y_n     = chi[-1]*jnp.cos(theta_n)          
    eta_n   = jnp.abs(y_n/chi[-1])            
    x_cos   = c_r*(1 - eta_n**p)**q + c_t*eta_n  
    x_lin   = jnp.interp(chi,eta_n, x_cos)  
    return x_lin 


# ----------------------------------------------------------------------
#   Run the Rotor Hover
# ---------------------------------------------------------------------- 
@jit
def run_rotor_OEI(nexus):
    
    # unpack
    rotor = nexus.vehicle_configurations.oei.networks.Battery_Electric_Rotor.rotors.prop_rotor
    rs    = nexus.rotor_settings

    # Setup Test conditions
    speed    = rs.OEI.design_freestream_velocity 
    altitude = np.array([rs.OEI.design_altitude])
    TD       = rs.temperature_deviation
    R        = rotor.tip_radius
    TM       = nexus.tip_mach_OEI
    

    # Calculate the atmospheric properties
    atmosphere            = MARC.Analyses.Atmospheric.US_Standard_1976()
    atmosphere_conditions = atmosphere.compute_values(altitude,temperature_deviation=TD)

    # Pack everything up
    conditions                                          = MARC.Analyses.Mission.Segments.Conditions.Aerodynamics()
    conditions.freestream.update(atmosphere_conditions)
    conditions.frames.inertial.velocity_vector          = jnp.array([[0.,0.,speed]])
    conditions.propulsion.throttle                      = jnp.array([[0.8]])
    conditions.frames.body.transform_to_inertial        = jnp.array([[[1., 0., 0.],[0., 1., 0.],[0., 0., -1.]]])
    
    # Calculate the RPM
    tip_speed = atmosphere_conditions.speed_of_sound*TM
    RPM       = tip_speed/R
    
    # Set the RPM
    rotor.inputs.omega = jnp.array(RPM)    
    
    # Run the rotor
    F, Q, P, Cp, outputs, etap  = prop_class.spin(rotor,conditions)
    
    # Pack the results
    nexus.results.OEI.thrust       = -F[0,-1]
    nexus.results.OEI.torque       = Q
    nexus.results.OEI.power        = P
    nexus.results.OEI.power_c      = Cp
    nexus.results.OEI.collective   = rotor.inputs.pitch_command  
    nexus.results.OEI.full_results = outputs
    nexus.results.OEI.efficiency   = etap 
    nexus.results.OEI.conditions   = conditions
    
    return nexus

# ----------------------------------------------------------------------
#   Run the Rotor Hover
# ---------------------------------------------------------------------- 
@jit
def run_rotor_hover(nexus):
    
    # unpack
    rotor = nexus.vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor
    rs    = nexus.rotor_settings

    # Setup Test conditions
    speed    = rs.hover.design_freestream_velocity 
    altitude = np.array([rs.hover.design_altitude])
    TD       = rs.temperature_deviation
    R        = rotor.tip_radius
    TM       = nexus.tip_mach_hover
    

    # Calculate the atmospheric properties
    atmosphere            = MARC.Analyses.Atmospheric.US_Standard_1976()
    atmosphere_conditions = atmosphere.compute_values(altitude,temperature_deviation=TD)

    # Pack everything up
    conditions                                          = MARC.Analyses.Mission.Segments.Conditions.Aerodynamics()
    conditions.freestream.update(atmosphere_conditions)
    conditions.frames.inertial.velocity_vector          = jnp.array([[0.,0.,speed]])
    conditions.propulsion.throttle                      = jnp.array([[0.8]])
    conditions.frames.body.transform_to_inertial        = jnp.array([[[1., 0., 0.],[0., 1., 0.],[0., 0., -1.]]])

    # Calculate the RPM
    tip_speed = atmosphere_conditions.speed_of_sound*TM
    RPM       = tip_speed/R
    
    # Set the RPM
    rotor.inputs.omega = jnp.array(RPM)    
    
    # Run the rotor
    F, Q, P, Cp, outputs, etap  = prop_class.spin(rotor,conditions)
    
    # Pack the results
    nexus.results.hover.thrust       = -F[0,-1]
    nexus.results.hover.torque       = Q
    nexus.results.hover.power        = P
    nexus.results.hover.power_c      = Cp
    nexus.results.hover.collective   = rotor.inputs.pitch_command 
    nexus.results.hover.full_results = outputs
    nexus.results.hover.efficiency   = etap 
    nexus.results.hover.conditions   = conditions
    
    return nexus


# ----------------------------------------------------------------------
#   Run the Rotor Forward
# ---------------------------------------------------------------------- 
@jit
def run_rotor_forward(nexus):
    
    # unpack
    rotor = nexus.vehicle_configurations.cruise.networks.Battery_Electric_Rotor.rotors.prop_rotor
    rs    = nexus.rotor_settings

    # Setup Test conditions
    speed    = rs.cruise.freeestream_velocity 
    altitude = np.array([rs.cruise.design_altitude])
    TD       = rs.temperature_deviation
    R        = rotor.tip_radius
    TM       = nexus.tip_mach_cruise
    
    # Calculate the atmospheric properties
    atmosphere            = MARC.Analyses.Atmospheric.US_Standard_1976()
    atmosphere_conditions = atmosphere.compute_values(altitude,temperature_deviation=TD)

    # Pack everything up
    conditions                                          = MARC.Analyses.Mission.Segments.Conditions.Aerodynamics()
    conditions.freestream.update(atmosphere_conditions)
    conditions.frames.inertial.velocity_vector          = jnp.array([[0.,0.,speed]])
    conditions.propulsion.throttle                      = jnp.array([[0.8]])
    conditions.frames.body.transform_to_inertial        = jnp.array([[[1., 0., 0.],[0., 1., 0.],[0., 0., -1.]]])
    
    # Calculate the RPM
    tip_speed = atmosphere_conditions.speed_of_sound*TM
    RPM       = tip_speed/R
    
    # Set the RPM
    rotor.inputs.omega = jnp.array(RPM)    
    
    # Run the rotor
    F, Q, P, Cp, outputs, etap  = prop_class.spin(rotor,conditions)
    
    # Pack the results
    nexus.results.cruise.thrust       = -F[0,-1]
    nexus.results.cruise.torque       = Q
    nexus.results.cruise.power        = P
    nexus.results.cruise.power_c      = Cp
    nexus.results.cruise.full_results = outputs
    nexus.results.cruise.collective   = rotor.inputs.pitch_command 
    nexus.results.cruise.efficiency   = etap 
    nexus.results.cruise.conditions   = conditions

    return nexus

# ----------------------------------------------------------------------
#   Run the hover noise
# ---------------------------------------------------------------------- 
@jit
def run_hover_noise(nexus):
    
    # unpack
    rotor        = nexus.vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor
    rotors       = nexus.vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors
    rs           = nexus.rotor_settings    
    conditions   = nexus.results.hover.conditions 
    full_results = nexus.results.hover.full_results

    # microphone locations
    altitude            = rs.hover.design_altitude
    ctrl_pts            = 1 
    theta               = rs.microphone_evaluation_angle
    S_hover             = jnp.maximum(altitude,20*Units.feet)  
    mic_positions_hover = jnp.array([[0.0 , S_hover*jnp.sin(theta)  ,S_hover*jnp.cos(theta)]])      
    
    # Run noise model  
    conditions.noise.total_microphone_locations      = jnp.repeat(mic_positions_hover[ jnp.newaxis,:,: ],1,axis=0)
    conditions.aerodynamics.angle_of_attack          = np.ones((ctrl_pts,1))* 0. * Units.degrees 
    segment                                          = Segment() 
    segment.state.conditions                         = conditions
    #segment.state.conditions.expand_rows(ctrl_pts)  
    noise                                            = MARC.Analyses.Noise.Fidelity_Zero() 
    settings                                         = noise.settings   
    num_mic                                          = len(conditions.noise.total_microphone_locations[0])  
    conditions.noise.number_of_microphones           = num_mic   
    
    sk                                               = list(settings.keys())
    settings.static_keys                             = []
    for k in sk:
        if not isinstance(settings[k],bool):
            settings.static_keys.append(k)
    
    propeller_noise_hover                            = propeller_mid_fidelity(rotors,full_results,segment,settings)   
    
    mean_SPL_hover                                   = jnp.mean(propeller_noise_hover.SPL_dBA)    
    
    # Pack
    nexus.results.hover.mean_SPL   = mean_SPL_hover 
    nexus.results.hover.noise_data = propeller_noise_hover 

    return nexus


# ----------------------------------------------------------------------
#   Run the forward noise
# ---------------------------------------------------------------------- 
@jit
def run_forward_noise(nexus):
    
    # unpack
    rotor        = nexus.vehicle_configurations.cruise.networks.Battery_Electric_Rotor.rotors.prop_rotor
    rs           = nexus.rotor_settings    
    conditions   = nexus.results.cruise.conditions 
    full_results = nexus.results.cruise.full_results

    # microphone locations
    altitude            = rs.cruise.design_altitude
    ctrl_pts            = 1 
    theta               = rs.microphone_evaluation_angle  
    S_hover             = jnp.maximum(altitude,20*Units.feet)  
    mic_positions_hover = jnp.array([[0.0 , S_hover*np.sin(theta)  ,S_hover*np.cos(theta)]])      
    
    # Run noise model  
    conditions.noise.total_microphone_locations      = jnp.repeat(mic_positions_hover[ jnp.newaxis,:,: ],1,axis=0)
    conditions.aerodynamics.angle_of_attack          = np.ones((ctrl_pts,1))* 0. * Units.degrees 
    segment                                          = Segment() 
    segment.state.conditions                         = conditions
    segment.state.conditions.expand_rows(ctrl_pts)  
    noise                                            = MARC.Analyses.Noise.Fidelity_Zero() 
    settings                                         = noise.settings   
    num_mic                                          = len(conditions.noise.total_microphone_locations[0])  
    conditions.noise.number_of_microphones           = num_mic   
    propeller_noise_hover                            = propeller_mid_fidelity(rotor,full_results,segment,settings)   
    mean_SPL_hover                                   = jnp.mean(propeller_noise_hover.SPL_dBA)    
    
    # Pack
    nexus.results.cruise.mean_SPL = mean_SPL_hover 
    nexus.results.cruise.noise_data = propeller_noise_hover 
    
    return nexus

# ----------------------------------------------------------------------
#   Post Process Results to give back to the optimizer
# ----------------------------------------------------------------------   
@jit
def post_process(nexus):
    
    #
    summary = nexus.summary
    rotor   = nexus.vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor
    full_results_hover  = nexus.results.hover.full_results
    full_results_cruise = nexus.results.cruise.full_results
    etap_cruise         = nexus.results.cruise.efficiency
    rs                  = nexus.rotor_settings
    
    # Use a P-norm on the Cl
    summary.max_sectional_cl_hover  = jnp.sum(full_results_hover.lift_coefficient**8)/8
    summary.max_sectional_cl_cruise = jnp.sum(full_results_cruise.lift_coefficient**8)/8
    mean_CL_hover                   = jnp.mean(full_results_hover.lift_coefficient[0])
    mean_CL_cruise                  = jnp.mean(full_results_cruise.lift_coefficient[0])
    
    # pack cl
    rotor.design_Cl_cruise = mean_CL_cruise
    rotor.design_Cl_hover  = mean_CL_hover    
    
    # blade taper consraint 
    c                               = rotor.chord_distribution
    blade_taper                     = c[-1]/c[0]

    # blade twist constraint  
    beta_blade                     = rotor.twist_distribution 

    # -------------------------------------------------------
    # OBJECTIVE FUNCTION
    # ------------------------------------------------------- 
    # Unpack the coefficients
    alpha = rs.alpha
    beta  = rs.beta
    gamma = rs.gamma
    
    # Efficiency and FOM
    ideal_eta_cruise = 1 
    ideal_FM_hover   = 1
    FM_hover         = full_results_hover.figure_of_merit
    eta_cruise       = etap_cruise[0][0] 
    
    # Acoustics
    Acoustic_Metric_cruise = nexus.results.cruise.mean_SPL
    Acoustic_Metric_hover  = nexus.results.hover.mean_SPL
    ideal_SPL              = rs.ideal_SPL_dBA
    
    # Calculate the two objectives
    aero_objective  = LA.norm((FM_hover - ideal_FM_hover)*100/(ideal_FM_hover*100))*beta \
        +  LA.norm((eta_cruise - ideal_eta_cruise)*100/(ideal_eta_cruise*100))*(1-beta) 
    
    acous_objective = LA.norm((Acoustic_Metric_hover - ideal_SPL)/ideal_SPL)*gamma + \
         LA.norm((Acoustic_Metric_cruise - ideal_SPL)/ideal_SPL)*(1-gamma)
         

        
    #Final Packing  
    summary.blade_taper_constraint = blade_taper 
    summary.blade_twist_constraint = beta_blade[0] - beta_blade[-1] 
    summary.objective              = aero_objective*10*alpha + acous_objective*10*(1-alpha)    
    # q to p ratios 
    summary.chord_p_to_q_ratio = rotor.chord_p/rotor.chord_q
    summary.twist_p_to_q_ratio = rotor.twist_p/rotor.twist_q    

   
    return nexus    
