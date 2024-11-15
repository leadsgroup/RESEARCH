# nonuniform_propeller_inflow.py

import RCAIDE 
from RCAIDE.Framework.Core import Units, Data 
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor         import design_propeller  
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor         import design_propeller   
from RCAIDE.Library.Methods.Aerodynamics.Common.Lift.compute_wing_wake import compute_wing_wake
from RCAIDE.Library.Methods.Aerodynamics.Common.Lift.compute_propeller_nonuniform_freestream import compute_propeller_nonuniform_freestream  
from RCAIDE.Library.Plots import * 

import os 
import numpy as np
import pylab as plt
 

def main():
    '''
    This example shows a propeller operating in three cases of nonuniform freestream flow:
    First, a propeller operates at a nonzero thrust angle relative to the freestream.
    Second, a propeller operates with an arbitrary upstream disturbance.
    Third, a propeller operates in the wake of an upstream wing

    '''
    # setup a simple vehicle
    Na=24
    Nr=101

    # setup the atmospheric conditions
    conditions = test_conditions()
     
    # set operating conditions for propeller test
    prop                                   = RCAIDE.Library.Components.Propulsors.Converters.Propeller() 
    prop.number_of_blades                  = 2
    prop.tip_radius                        = 38.    * Units.inches
    prop.hub_radius                        = 8.     * Units.inches
    prop.cruise.design_freestream_velocity = 135.   * Units['mph']
    prop.cruise.design_angular_velocity    = 1300.  * Units.rpm
    prop.cruise.design_Cl                  = 0.8
    prop.cruise.design_altitude            = 12000. * Units.feet
    prop.cruise.design_thrust              = 1200.
    prop.origin                            = [[0.,0.,0.]]
    prop.number_azimuthal_stations         = Na
    prop.rotation                          = 1
    prop.symmetry                          = True 
    airfoil                                = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil.tag                            = 'NACA_4412' 
    separator                             = os.path.sep 
    airfoil.coordinate_file                 =     '..' + separator + 'Airfoils' + separator  + 'NACA_4412.txt'   # absolute path   
    airfoil.polar_files                     = [   '..' + separator + 'Airfoils' + separator  + 'Polars' + separator  + 'NACA_4412_polar_Re_50000.txt',
                                                  '..' + separator + 'Airfoils' + separator  + 'Polars' + separator  + 'NACA_4412_polar_Re_100000.txt',
                                                  '..' + separator + 'Airfoils' + separator  + 'Polars' + separator  + 'NACA_4412_polar_Re_200000.txt',
                                                  '..' + separator + 'Airfoils' + separator  + 'Polars' + separator  + 'NACA_4412_polar_Re_500000.txt',
                                                  '..' + separator + 'Airfoils' + separator  + 'Polars' + separator  + 'NACA_4412_polar_Re_1000000.txt']   
     
    prop.append_airfoil(airfoil)      
 
    
    prop.airfoil_polar_stations    = list(np.zeros(Nr).astype(int))
    prop                           = design_propeller(prop,Nr)
    
    
    prop.inputs.omega              = np.ones_like(conditions.aerodynamics.angles.alpha)*prop.cruise.design_angular_velocity
    prop.orientation_euler_angles  = [0.,20.*Units.degrees,0]
    prop.use_2d_analysis           = True
    
    # spin propeller in nonuniform flow
    thrust, torque, power, Cp, outputs , etap = prop.spin(conditions)

    # plot velocities at propeller plane and resulting performance
    plot_rotor_disc_performance(prop,outputs,title='Case 1: Operating at Thrust Angle')
    
    thrust   = np.linalg.norm(thrust)
    thrust_r = 1743.0258191335301
    torque_r = 748.87304348
    power_r  = 101948.342247
    Cp_r     = 0.46942838
    etap_r   = 0.71821729
    print('\nCase 1 Errors: \n')
    print('Thrust difference = ', np.abs(thrust - thrust_r) / thrust_r )
    print('Torque difference = ', np.abs(torque - torque_r) / torque_r )
    print('Power difference = ', np.abs(power - power_r) / power_r )
    print('Cp difference = ', np.abs(Cp - Cp_r) / Cp_r )
    print('Etap difference = ', np.abs(etap - etap_r) / etap_r )
    assert (np.abs(thrust - thrust_r) / thrust_r < 1e-6), "Nonuniform Propeller Thrust Angle Regression Failed at Thrust Test"
    assert (np.abs(torque - torque_r) / torque_r < 1e-6), "Nonuniform Propeller Thrust Angle Regression Failed at Torque Test"
    assert (np.abs(power - power_r) / power_r < 1e-6), "Nonuniform Propeller Thrust Angle Regression Failed at Power Test"
    assert (np.abs(Cp - Cp_r) / Cp_r < 1e-6), "Nonuniform Propeller Thrust Angle Regression Failed at Power Coefficient Test"
    assert (np.abs(etap - etap_r) / etap_r < 1e-6), "Nonuniform Propeller Thrust Angle Regression Failed at Efficiency Test"

    return
 

def test_conditions():
    # --------------------------------------------------------------------------------------------------
    # Atmosphere Conditions:
    # --------------------------------------------------------------------------------------------------
    atmosphere = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmo_data  = atmosphere.compute_values(altitude=14000 * Units.ft)
    rho        = atmo_data.density
    mu         = atmo_data.dynamic_viscosity
    T          = atmo_data.temperature
    a          = atmo_data.speed_of_sound 

    # aerodynamics analyzed for a fixed angle of attack
    aoa   = np.array([[ 3 * Units.deg  ]])
    Vv    = np.array([[ 100 * Units.mph]])
    ones  = np.ones_like(aoa)

    mach  = Vv/a

    conditions                              = RCAIDE.Framework.Mission.Common.Results()
    conditions.freestream.density           = rho* ones
    conditions.freestream.dynamic_viscosity = mu* ones
    conditions.freestream.speed_of_sound    = a* ones
    conditions.freestream.temperature       = T* ones
    conditions.freestream.mach_number       = mach* ones
    conditions.freestream.velocity          = Vv * ones
    conditions.aerodynamics.angles.alpha    = aoa
    conditions.frames.body.transform_to_inertial = np.array( [[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]] ) 
    velocity_vector = np.zeros([len(aoa), 3])
    velocity_vector[:, 0] = Vv
    conditions.frames.inertial.velocity_vector = velocity_vector 

    return conditions
 
 
if __name__ == '__main__':
    main()
    plt.show()
