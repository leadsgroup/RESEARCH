# Prop_Rotor_Optimization 

from RCAIDE.Core import Units, Data, to_numpy
from RCAIDE.Optimization import Nexus
from RCAIDE.Optimization.Package_Setups.pyoptsparse_setup import Pyoptsparse_Solve 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    
import matplotlib.pyplot as plt
import jax.numpy as jnp
import pickle
import os 
import numpy as np   
import sys 
sys.path.append('../Common')

from Optimize_example import setup
from Optimize_example import write_vsp_BEM_out

def main():
    #alphas = np.linspace(0,1,5)
    #betas  = np.linspace(0,1,5)
    alphas = np.array([0.8,0.8,0.9,0.2])
    betas  = np.array([0.1,0.9,0.5,0.5])
    traditional_format =  False  
              
    # Matthew these settings are where you can change things
    problem = Nexus() 
    
    # Settings
    problem.rotor_settings = Data()    
    
    # Problem formulation settings
    rs = problem.rotor_settings
    rs.traditional_format = traditional_format
    rs.prop_rotor         = True
 
    rs.orientation_euler_angles   = [0,np.pi/2,0]
    rs.tip_radius                 = 3/2
    rs.hub_radius                 = 0.15 * rs.tip_radius  
    rs.number_of_blades           = 3    
    rs.design_angular_velocity           = 1500
    rs.design_Cl                  = 0.7
    rs.temperature_deviation      = 0.
    rs.number_of_stations         = 200
    rs.airfoil_polar_stations     = list(np.zeros(rs.number_of_stations)) 
    ospath    = os.path.abspath(__file__)
    separator = os.path.sep
    rel_path  = os.path.dirname(ospath) + separator  
    
    rs.airfoil_geometry           = [ rel_path + '../Common/Airfoils/NACA_4412.txt']
    rs.airfoil_polars             = [[rel_path + '../Common/Airfoils/Polars/NACA_4412_polar_Re_50000.txt',
                                      rel_path + '../Common/Airfoils/Polars/NACA_4412_polar_Re_100000.txt',
                                      rel_path + '../Common/Airfoils/Polars/NACA_4412_polar_Re_200000.txt',
                                      rel_path + '../Common/Airfoils/Polars/NACA_4412_polar_Re_500000.txt',
                                      rel_path + '../Common/Airfoils/Polars/NACA_4412_polar_Re_1000000.txt',
                                      rel_path + '../Common/Airfoils/Polars/NACA_4412_polar_Re_3500000.txt',
                                      rel_path + '../Common/Airfoils/Polars/NACA_4412_polar_Re_5000000.txt',
                                      rel_path + '../Common/Airfoils/Polars/NACA_4412_polar_Re_7500000.txt']]                                
    rs.static_keys                = ['airfoil_geometry','airfoil_polars','airfoil_polar_stations','temperature_deviation']
    
    # Hover Settings
    rs.hover  = Data()
    rs.hover.design_thrust               = 23175.5364/6  
    rs.hover.design_freestream_velocity         = np.sqrt(rs.hover.design_thrust/(2*1.2*np.pi*(rs.tip_radius**2))) # Ideal power  
    rs.hover.design_altitude             = 20 * Units.feet 
    rs.hover.design_tip_mach             = 0.5 
    rs.hover.static_keys                 = ['design_altitude']
    
    # Hover Settings
    rs.OEI  = Data()
    rs.OEI.design_thrust                 = 23175.5364/5  
    rs.OEI.design_freestream_velocity           = np.sqrt(rs.OEI.design_thrust/(2*1.2*np.pi*(rs.tip_radius**2))) # Ideal power  
    rs.OEI.design_altitude               = 20 * Units.feet 
    rs.OEI.design_tip_mach               = 0.5 
    rs.OEI.static_keys                   = ['design_altitude']    
 
    # Cruise Settings
    rs.cruise = Data()
    rs.cruise.freeestream_velocity       = 175*Units.mph 
    rs.cruise.design_thrust              = 4000/6
    rs.cruise.design_altitude            = 3000. * Units.feet
    rs.cruise.design_tip_mach            = 0.5 
    rs.cruise.static_keys                =  ['design_altitude']
                         
    rs.tip_mach_range                    = [0.3,0.7]
    rs.microphone_evaluation_angle       = 135 * Units.degrees  # 45 degrees behind rotor plane
    rs.slack_constraint                  = 1E-3  
    rs.ideal_SPL_dBA                     = 45 
    
    rs.alpha = jnp.array(0.5)
    rs.beta  = jnp.array(0.5)
    rs.gamma = jnp.array(1.0)      
    
    # Setup the problem
    problem = setup(problem)    
    problem.use_jax_derivatives = True
    problem.jitable             = True 
    problem.record_objective    = True
    
    # Higher precision for finite differencing
    from jax.config import config
    config.update("jax_enable_x64", True)  
    
    sys.stdout.flush()
              
    for i in range(len(alphas)): 
            
        # Optimization problem settings
        problem.rotor_settings.alpha = jnp.array(alphas[i])
        problem.rotor_settings.beta  = jnp.array(betas[i])
        problem.rotor_settings.gamma = jnp.array(1.)
        problem.write_file = 'optimization_outputs' + str(alphas[i])+str(betas[i]) + '.txt'     

        # now optimize the function
        outputs = Pyoptsparse_Solve(problem,solver="SNOPT")
        print(outputs)
                
        # Write the rotor back
        problem.vehicle_configurations.hover.networks = to_numpy(problem.vehicle_configurations.hover.networks)
        
        # Alias
        rotor = problem.vehicle_configurations.hover.networks.Battery_Electric_Rotor.rotors.prop_rotor
     
        # Output a file
        write_vsp_BEM_out(rotor, problem.results.hover.full_results,filename = 'BEMT_out' + str(alphas[i]) + str(betas[i]) +'.bem') 
          
        # append final rotor performance to data structure 
        rotor.results     = problem.results 
         
        # save rotor geomtry 
        alpha_opt_weight = str(format(alphas[i],'.5f'))
        alpha_opt_weight = alpha_opt_weight.replace('.','_')    
        beta_opt_weight  = str(format(betas[i],'.5f'))
        beta_opt_weight  = beta_opt_weight.replace('.','_')  
        name       = 'PR_Alpha_' + alpha_opt_weight + '_Beta_' + beta_opt_weight  
        save_blade_geometry(rotor,name)  
        sys.stdout.flush()

    return 

def save_blade_geometry(rotor,filename): 
    ospath    = os.path.abspath(__file__)
    separator = os.path.sep
    rel_path  = os.path.dirname(ospath) + separator     
    pickle_file  = rel_path  + filename + '.pkl'
    with open(pickle_file, 'wb') as file:
        pickle.dump(rotor, file) 
    return     


if __name__ == '__main__':
    #with jax.disable_jit():
    main()
    plt.show()