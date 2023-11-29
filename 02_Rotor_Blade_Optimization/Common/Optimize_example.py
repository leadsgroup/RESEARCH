# Optimize.py
# 
# Created: Jun 2022, E. Botero
# Modified: 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

from MARC.Core import Units, Data, to_numpy
import numpy as np 
import Vehicles
import Procedure
from MARC.Optimization import Nexus
from MARC.Optimization.Package_Setups.pyoptsparse_setup import Pyoptsparse_Solve
#from MARC.Input_Output.OpenVSP import write
from MARC.Plots.Geometry import plot_propeller
import matplotlib.pyplot as plt
import jax

# ----------------------------------------------------------------------        
#   Run the whole thing
# ----------------------------------------------------------------------  
def main():

    # Matthew these settings are where you can change things
    problem = Nexus() 
    
    # Settings
    problem.rotor_settings = Data()    
    
    # Problem formulation settings
    rs = problem.rotor_settings
    rs.traditional_format = True
    rs.prop_rotor         = False    

    rs.orientation_euler_angles   = [0,np.pi/2,0]
    rs.tip_radius                 = 2.7/2
    rs.hub_radius                 = 0.15 * rs.tip_radius  
    rs.number_of_blades           = 3    
    rs.design_angular_velocity           = 1500
    rs.design_Cl                  = 0.7
    rs.temperature_deviation      = 0.
    rs.number_of_stations         = 200
    rs.airfoil_polar_stations     = list(np.zeros(rs.number_of_stations))
    rs.airfoil_geometry           = ['./Airfoils/NACA_4412.txt']
    rs.airfoil_polars             = [['./Airfoils/Polars/NACA_4412_polar_Re_50000.txt',
                                         './Airfoils/Polars/NACA_4412_polar_Re_100000.txt',
                                         './Airfoils/Polars/NACA_4412_polar_Re_200000.txt',
                                         './Airfoils/Polars/NACA_4412_polar_Re_500000.txt',
                                         './Airfoils/Polars/NACA_4412_polar_Re_1000000.txt',
                                         './Airfoils/Polars/NACA_4412_polar_Re_3500000.txt',
                                         './Airfoils/Polars/NACA_4412_polar_Re_5000000.txt',
                                         './Airfoils/Polars/NACA_4412_polar_Re_7500000.txt']]                               
    rs.static_keys                = ['airfoil_geometry','airfoil_polars','airfoil_polar_stations','temperature_deviation']
    
    # Hover Settings
    rs.hover  = Data()
    rs.hover.design_thrust         = 23544/8  # based on Stopped-Rotor V2, vehicle weight/(number of rotors - 1 )
    rs.hover.design_freestream_velocity   = np.sqrt(rs.hover.design_thrust/(2*1.2*np.pi*(rs.tip_radius**2))) # Ideal power  
    rs.hover.design_altitude       = 20 * Units.feet 
    rs.hover.design_tip_mach       = 0.5 
    rs.hover.static_keys           = ['design_altitude']
    
    # Hover Settings
    rs.OEI  = Data()
    rs.OEI.design_thrust         = 23544/6  # based on Stopped-Rotor V2, vehicle weight/(number of rotors - 1 )
    rs.OEI.design_freestream_velocity   = np.sqrt(rs.OEI.design_thrust/(2*1.2*np.pi*(rs.tip_radius**2))) # Ideal power  
    rs.OEI.design_altitude       = 20 * Units.feet 
    rs.OEI.design_tip_mach       = 0.5 
    rs.OEI.static_keys           = ['design_altitude']    

    # Cruise Settings
    rs.cruise = Data()
    rs.cruise.freeestream_velocity = 40.
    rs.cruise.design_thrust        = 23544/12
    rs.cruise.design_altitude      = 3000. * Units.feet
    rs.cruise.design_tip_mach      = 0.5 
    rs.cruise.static_keys          = ['design_altitude']
                         
    rs.tip_mach_range                    = [0.3,0.7]
    rs.microphone_evaluation_angle       = 135 * Units.degrees  # 45 degrees behind rotor plane
    rs.slack_constraint                  = 1E-3  
    rs.ideal_SPL_dBA                     = 45 
    
    # Optimization problem settings
    rs.alpha = 0.5
    rs.beta  = 0.5
    rs.gamma = 0.5
    
    # Ideally we won't need to touch below here
    
    # Setup the problem
    problem = setup(problem)    
    
    
    # Higher precision for finite differencing
    from jax.config import config
    config.update("jax_enable_x64", True)
    
    problem.use_jax_derivatives = True
    problem.jitable             = True
    problem.record_objective    = False
    
    
    ###### Derivative checks:
    
    ## Get the current values:
    #input_array = problem.optimization_problem.inputs.pack_array()
    #x = input_array[0::5]/input_array[3::5]       
        
    ## Doing a forward pass
    #print('Objective Value')
    #output = problem.objective(x)    
    #print(output)      
    
    #print('Constraint Value')
    #output = problem.all_constraints(x)    
    #print(output)         
    
    ## try taking a gradient
    #print('Gradient of Objective')
    #grad1 = problem.grad_objective(x)
    #print(grad1)
    
    ## try taking a jacobian of the constraints
    #print('Jacobian of Constraints')
    #jac1 = problem.jacobian_all_constraints(x)
    #print(jac1)
    
    ## Compare against FD gradients
    #print('FD Gradients')
    #fd_grad, fd_jac = problem.finite_difference(x)
    #print(fd_grad)
    #print(fd_jac)
    #print('FD Grad Error')
    #print(np.nan_to_num(np.abs((fd_grad-grad1)/fd_grad),posinf=0))
    #print('Max FD Jacobian Error')
    #print(np.nan_to_num(np.abs((fd_jac-jac1)/fd_jac),posinf=0))   
    
    
    ###### Derivative checks Done

    # now optimize the function
    outputs = Pyoptsparse_Solve(problem,solver="SLSQP")
    print(outputs)
            
    # Write the rotor back
    problem.vehicle_configurations.hover.networks = to_numpy(problem.vehicle_configurations.hover.networks)
    
    # Alias
    rotor = problem.vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor

    # Plot the rotor
    #plot_propeller(rotor)
    
    # Output a file
    write_vsp_BEM_out(rotor, problem.results.hover.full_results)

    return

# ----------------------------------------------------------------------        
#   Inputs, Objective, & Constraints
# ----------------------------------------------------------------------  

def setup(nexus):

    problem = nexus.optimization_problem
    rs      = nexus.rotor_settings

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #   [ tag       , initial,     lb , ub        , scaling , units ]
    problem.inputs = np.array([
        ['root_incidence'   , 0, -np.pi/6   , np.pi/6, 1., 1.],
        ['collective_cruise', 0, -np.pi/6   , np.pi/6, 1., 1.],
        ['collective_OEI'   , 0, -np.pi/6   , np.pi/6, 1., 1.],
        ['tip_mach_hover'   , rs.hover.design_tip_mach, rs.tip_mach_range[0],rs.tip_mach_range[1], 1., 1.],
        ['tip_mach_OEI'     , rs.hover.design_tip_mach, rs.tip_mach_range[0],0.85, 1., 1.],
    ],dtype=object)
    
    if rs.prop_rotor:
        new_inputs = np.array([
            ['tip_mach_cruise', rs.cruise.design_tip_mach, rs.tip_mach_range[0],rs.tip_mach_range[1], 1., 1.],
        ],dtype=object)        
        problem.inputs = np.vstack((problem.inputs,new_inputs))

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([
         [ 'objective', 1. , 1.*Units.less],
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    SL = rs.slack_constraint
    # [ tag, sense, edge, scaling, units ]
    problem.constraints = np.array([
         [ 'thrust_hover',          '>', rs.hover.design_thrust-SL,rs.hover.design_thrust, 1.*Units.newtons],
         [ 'thrust_OEI',            '>', rs.OEI.design_thrust,rs.OEI.design_thrust, 1.*Units.newtons],
         [ 'thrust_hover_2',        '<', rs.hover.design_thrust+SL,rs.hover.design_thrust, 1.*Units.newtons],
         [ 'taper',                 '>', 0.3, 1.0, 1.*Units.less],
         [ 'taper',                 '<', 0.9, 1.0, 1.*Units.less],
         [ 'max_cl_hover',          '<', 1.2, 1.0, 1.*Units.less],
         [ 'twist_pq',              '>', 0.5, 1.0, 1.*Units.less],
         [ 'chord_pq',              '>', 0.5, 1.0, 1.*Units.less],
    ],dtype=object)
    # need to add boolean to add extra residual constraints for prop rotor     
    if rs.prop_rotor:
        
        new_constraints = np.array([
             [ 'thrust_cruise',     '>', rs.cruise.design_thrust-SL,rs.hover.design_thrust, 1.*Units.newtons],
             [ 'thrust_cruise_2',   '<', rs.cruise.design_thrust+SL,rs.hover.design_thrust, 1.*Units.newtons],
            ],dtype=object)  
        
        problem.constraints = np.vstack((problem.constraints,new_constraints))
        
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] 
    problem.aliases = [
        [ 'objective'        , 'summary.objective'    ],
        [ 'thrust_hover'     , 'results.hover.thrust'       ],    
        [ 'thrust_hover_2'   , 'results.hover.thrust'       ],    
        [ 'thrust_cruise'    , 'results.cruise.thrust'       ],    
        [ 'thrust_cruise_2'  , 'results.cruise.thrust'       ],  
        [ 'thrust_OEI'       , 'results.OEI.thrust'       ],  
        [ 'collective_cruise', 'vehicle_configurations.cruise.networks.Battery_Electric_Rotor.rotors.prop_rotor.inputs.pitch_command'    ],
        [ 'collective_OEI'   , 'vehicle_configurations.oei.networks.Battery_Electric_Rotor.rotors.prop_rotor.inputs.pitch_command'    ],
        [ 'root_incidence'   , 'rotor_settings.root_incidence'       ],            
        [ 'tip_mach_hover'   , 'tip_mach_hover'       ],      
        [ 'tip_mach_OEI'     , 'tip_mach_OEI'       ],      
        [ 'tip_mach_cruise'  , 'tip_mach_cruise'       ],      
        [ 'taper'            , 'summary.blade_taper_constraint'   ],            
        [ 'max_cl_hover'     , 'summary.max_sectional_cl_hover'   ],               
        [ 'twist_pq'         , 'summary.twist_p_to_q_ratio'       ],            
        [ 'chord_pq'         , 'summary.chord_p_to_q_ratio'       ],               
    ]
    # need to add boolean to add extra aliases constraints for prop rotor   
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles.setup(nexus)
    R  = rs.tip_radius
    # Add extra variables
    if rs.traditional_format:
        nexus.add_array_inputs('vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.chord_distribution', 0.05*R, 0.2*R)
        nexus.add_array_inputs('vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.twist_distribution', -np.pi/4, np.pi/2)
    else:
        new_inputs = np.array([
                            ['chord_r',  0.1*R    , 0.05*R     , 0.2*R     , 1., 1.],
                            ['chord_p',  2        , 0.25       , 2.0       , 1., 1.],
                            ['chord_q',  1        , 0.25       , 1.5       , 1., 1.],
                            ['chord_t',  0.05*R   , 0.02*R     , 0.2*R     , 1., 1.],
                            ['twist_r',  np.pi/6  , -np.pi/4   , np.pi/2   , 1., 1.],
                            ['twist_p',  1        , 0.25       , 2.0       , 1., 1.],
                            ['twist_q',  0.5      , 0.25       , 1.5       , 1., 1.],
                            ['twist_t',  np.pi/10 , -np.pi/4   , np.pi/2   , 1., 1.],
        ],dtype=object)
        
        new_aliases = [
            ['chord_r', 'vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.chord_r'],
            ['chord_p', 'vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.chord_p'],
            ['chord_q', 'vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.chord_q'],
            ['chord_t', 'vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.chord_t'],
            ['twist_r', 'vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.twist_r'],
            ['twist_p', 'vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.twist_p'],
            ['twist_q', 'vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.twist_q'],
            ['twist_t', 'vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.twist_t'],
    ]
        problem.inputs  = np.vstack((problem.inputs,new_inputs))
        problem.aliases = problem.aliases + new_aliases
    
    # convert these to the new style
    nexus.convert_problem_arrays()    
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.setup(rs)
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------  
    nexus  = setup_dictionaries(nexus)

    return nexus



# ----------------------------------------------------------------------        
#   Setup Variables
# ----------------------------------------------------------------------   

def setup_dictionaries(nexus):
    
    nexus.summary = Data()    
    nexus.summary.nothing = 0.
    nexus.total_number_of_iterations = 0.
    
    nexus.results.hover  = Data()
    nexus.results.cruise = Data()
    nexus.results.OEI    = Data()
    
    nexus.results.cruise.full_results = Data()
    nexus.results.cruise.full_results.lift_coefficient = np.array([[0.]])
    nexus.results.cruise.efficiency                    = np.array([[0.]])
    nexus.results.cruise.mean_SPL                      = 0
    nexus.results.hover.mean_SPL                       = 0
    
    nexus.vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.chord_p = 1.
    nexus.vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.chord_q = 1.
    nexus.vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.twist_p = 1.
    nexus.vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor.twist_q = 1.    
    
    
    return nexus

# ----------------------------------------------------------------------
#   write_BEMT_out
# ----------------------------------------------------------------------


def write_vsp_BEM_out(rotor,rotor_outputs,filename='BEMT_out.bem'):
    
    #pull out the required data
    Nr    = int(rotor_outputs.number_radial_stations)
    B     = int(rotor.number_of_blades)
    R     = rotor.tip_radius/ Units.ft
    r     = rotor.radius_distribution / Units.ft
    r_R   = r /R
    c_R   = (rotor.chord_distribution/ Units.ft)/ R
    twist = rotor.twist_distribution / Units.degrees
    t_c   = rotor.thickness_to_chord
    origin = rotor.origin[0]
    
    rake  = np.zeros_like(r)
    skew  = np.zeros_like(r)
    sweep = np.zeros_like(r)
    CLi   = np.ones_like(r)
    tan   = np.zeros_like(r)
    axial = np.zeros_like(r)
    

    f = open(filename,'w')
    f.write('...'+rotor.tag+'...\n')
    f.write('Num_Sections:'+str(Nr)+'\n')
    f.write('Num_Blade: ' + str(B) +'\n')
    f.write('Diameter: ' + str(2*R) +'\n')
    f.write('Beta 3/4 (deg): 0.\n')
    f.write('Feather (deg): 0.00000000\n')
    f.write('Pre_Cone (deg): 0.00000000\n')
    f.write('Center: '+str(origin[0]/ Units.ft)+','+str(origin[1]/ Units.ft)+','+str(origin[2]/ Units.ft)+'\n')    
    f.write('Normal: 1.0, 0.0, 0.0\n')
    f.write('Radius/R, Chord/R, Twist (deg), Rake/R, Skew/R, Sweep, t/c, CLi, Axial, Tangential\n')
    
    for ii in range(Nr):
        string = ''
        string = string + '{0: >10}'.format(str('%7.6f' %r_R[ii])) + ', '
        string = string + '{0: >10}'.format(str('%7.6f' %c_R[ii])) + ', '
        string = string + '{0: >10}'.format(str('%7.6f' %twist[ii])) + ', '
        string = string + '{0: >10}'.format(str('%7.6f' %rake[ii])) + ', '
        string = string + '{0: >10}'.format(str('%7.6f' %skew[ii])) + ', '
        string = string + '{0: >10}'.format(str('%7.6f' %sweep[ii])) + ', '
        string = string + '{0: >10}'.format(str('%7.6f' %t_c[ii])) + ', '
        string = string + '{0: >10}'.format(str('%7.6f' %CLi[ii])) + ', '
        string = string + '{0: >10}'.format(str('%7.6f' %axial[ii])) + ', '        
        string = string + '{0: >12}'.format(str('%7.6f' %tan[ii])) +'\n'

        f.write(string)
        
    f.close()

    return

if __name__ == '__main__':
    #with jax.disable_jit():
    main()
    plt.show()
