# motor_test.py
# 
 
#----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import RCAIDE
from RCAIDE.Core import Units
from RCAIDE.Core import Data, Container
from RCAIDE.Methods.Propulsion.electric_motor_sizing import size_from_mass , size_optimal_motor 
import numpy as np 
import pickle 
import matplotlib.pyplot as plt  
import matplotlib.cm as cm 

def main(): 
    '''
    Plot Parameters
    '''
    

    # Universal Plot Settings 
    plt.rcParams['axes.linewidth'] = 2.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 22,
                  'xtick.labelsize': 18,
                  'ytick.labelsize': 18,
                  'axes.titlesize': 22}
    plt.rcParams.update(parameters)
    plot_parameters                  = Data()
    plot_parameters.line_width       = 3 
    plot_parameters.line_style       = '-' 
    plot_parameters.figure_width     = 10 
    plot_parameters.figure_height    = 7 
    plot_parameters.marker_size      = 10 
    plot_parameters.legend_font_size = 20 
    plot_parameters.plot_grid        = True   
    plot_parameters.markers          = ['o','v','s','^','p','^','D','X','*']
    plot_parameters.colors           = cm.inferno(np.linspace(0,1,5))     
    plot_parameters.lw               = 2                              # line_width                
    plot_parameters.legend_font      = 20                             # legend_font_size       
    
    '''
    Design Old Motor 
    ''' 
    motor                              = RCAIDE.Components.Energy.Converters.Motor()
    motor.mass_properties.mass         = 9. * Units.kg 
    motor.efficiency                   = 0.935
    motor.gear_ratio                   = 1. 
    motor.gearbox_efficiency           = 1. # Gear box efficiency     
    motor.no_load_current              = 2.0 
    motor.nominal_voltage              = 400  
    motor.rotor_radius                 = 1.5  # prop.tip_radius
    motor.design_torque                = 1000 # prop.cruise.design_torque
    motor.angular_velocity             = 100  # prop.cruise.design_angular_velocity/motor.gear_ratio
    motor                              = size_optimal_motor(motor)   
    
    '''
    Design New Motor 
    ''' 
    motor_2                             = RCAIDE.Components.Energy.Converters.Motor_2()
    
    '''
    Operating conditions
    '''
 
    design_altitude = 1000 # m 
    flight_velocity = 100  # m/s   
    propeller_power_coefficient = 0.02
    
    atmosphere                                          = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmosphere_conditions                               = atmosphere.compute_values(design_altitude)  
    conditions                                          = Data()
    conditions.freestream                               = Data()
    conditions.propulsion                               = Data()
    conditions.frames                                   = Data()
    conditions.frames.body                              = Data()
    conditions.frames.inertial                          = Data()
    conditions.freestream.update(atmosphere_conditions)
    conditions.freestream.dynamic_viscosity             = atmosphere_conditions.dynamic_viscosity
    conditions.freestream.velocity                      = np.array([[flight_velocity,0,0]])
    conditions.propulsion.throttle                      = np.array([[1.0]])
    conditions.frames.body.transform_to_inertial        = np.array([np.eye(3)])  
    
    #motor_omega_function(motor,conditions)
    #motor_current_function(motor,conditions)
    
    
    motor_torque_function(motor_2,conditions) 
    
    return 

#------------------------------------
# Motor Omega Function  
#------------------------------------
def motor_omega_function(motor_0,conditions):
    
    # create copy of motor to test functions 
    motor = motor_0    
    
    # Define function specific inputs  
    motor.inputs.voltage  = np.atleast_2d(np.linspace(0,500,100)).T
    motor.inputs.rotor_CP = np.ones_like(motor_0.inputs.voltage)*motor.angular_velocity
    
    # Run Motor Omega Function 
    omega  = motor.omega(conditions)   
    torque = motor.outputs.torque  
    
    fig_1 = plt.figure('Omega_Function')
    fig_1.set_size_inches(10,6)
    axis_1_1 = fig_1.add_subplot(1,2,1)
    axis_1_2 = fig_1.add_subplot(1,2,2) 
    
    axis_1_1.plot(motor.inputs.voltage[:,0],omega[:,0],)
    axis_1_1.set_xlim([0,500])
    axis_1_1.set_xlabel('Voltage')
    axis_1_1.set_ylabel('Angular Velocity')

    axis_1_2.plot(motor.inputs.voltage[:,0],torque[:,0]) 
    axis_1_2.set_xlim([0,500])
    axis_1_2.set_xlabel('Voltage')
    axis_1_2.set_ylabel('Torque')
    
    fig_1.tight_layout()
    return 


#------------------------------------
# Motor Current Function 
#------------------------------------ 
def motor_current_function(motor_0,conditions):   
    # create copy of motor to test functions 
    motor = motor_0
    
    # Define function specific inputs   
    motor.inputs.voltage = np.atleast_2d(np.linspace(0,500,100)).T
    motor.outputs.omega  = np.ones_like(motor_0.inputs.voltage)*motor_0.angular_velocity
    
    # Run Motor Current Function 
    i, etam = motor.current(conditions) 

    fig_2 = plt.figure('Current_Function')
    fig_2.set_size_inches(10,6)
    axis_2_1 = fig_2.add_subplot(1,2,1)
    axis_2_2 = fig_2.add_subplot(1,2,2) 
    
    axis_2_1.plot(motor.inputs.voltage[:,0], i[:,0])
    axis_2_1.set_xlim([0,500])
    axis_2_1.set_xlabel('Voltage')
    axis_2_1.set_ylabel('Current')
    
    axis_2_2.plot(motor.inputs.voltage[:,0], etam[:,0])
    axis_2_2.set_xlim([0,500])
    axis_2_2.set_xlabel('Voltage')
    axis_2_2.set_ylabel('Efficiency')
    
    fig_2.tight_layout()
    return
def motor_torque_function(motor_0,conditions):

    # create copy of motor to test functions 
    motor = motor_0    
    
    # Define function specific inputs  
    motor.inputs.curent   = np.atleast_2d(np.linspace(0,500,100)).T 
    
    # Run Motor Omega Function 
    torque  = motor.torque(conditions)    
    
    fig_1 = plt.figure('Omega_Function')
    fig_1.set_size_inches(10,6)
    axis_1_1 = fig_1.add_subplot(1,2,1)
    axis_1_2 = fig_1.add_subplot(1,2,2) 
     

    axis_1_2.plot(motor.inputs.curent[:,0],torque[:,0]) 
    axis_1_2.set_xlim([0,500])
    axis_1_2.set_xlabel('Current')
    axis_1_2.set_ylabel('Torque')
    
    fig_1.tight_layout()    
    return 
    
# ----------------------------------------------------------------------        
#   Call Main
# ----------------------------------------------------------------------    

if __name__ == '__main__':
    main()
    plt.show()