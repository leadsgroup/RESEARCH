# motor_test.py
# 
 
#----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import MARC
from MARC.Core import Units
from MARC.Core import Data, Container
from MARC.Methods.Propulsion.electric_motor_sizing import size_from_mass , size_optimal_motor 
import numpy as np 
import pylab as plt 

def main(): 
    '''
    Design Motor 
    ''' 
    motor                              = MARC.Components.Energy.Converters.Motor()
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
    Operating conditions
    '''
 
    design_altitude = 1000 # m 
    flight_velocity = 100  # m/s   
    propeller_power_coefficient = 0.02
    
    atmosphere                                          = MARC.Analyses.Atmospheric.US_Standard_1976()
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
    conditions.propulsion.propeller_power_coefficient   = np.array([[propeller_power_coefficient]]) 
    
    return 

#------------------------------------
# Motor Omega Function  
#------------------------------------
def motor_omega_function(motor,conditions):
    
    # create copy of motor to test functions 
    motor_1 = motor    
    
    # Define function specific inputs 
    voltage_1 =  400
    motor_1.inputs.voltage  = np.array([[voltage_1]]) 
    motor_1.inputs.rotor_CP = conditions.propulsion.propeller_power_coefficient
    
    # Run Motor Omega Function 
    omega_1  = motor_1.omega(conditions)   
    torque_1 = motor_1.outputs.torque[0][0]  
    
    fig_1 = plt.figure('Omega_Function')
    axis_1_1 = fig_1.add_subplot(1,2,1)
    axis_1_2 = fig_1.add_subplot(1,2,2) 
    axis_1_1.plot(omega_1,voltage_1)
    axis_1_2.plot(torque_1,voltage_1) 
    
    return 

#------------------------------------
# Motor Current Function 
#------------------------------------ 
def motor_current_function(motor,conditions):   
    # create copy of motor to test functions 
    motor_2 = motor    
    
    # Define function specific inputs  
    voltage_1 =  400
    motor_2.inputs.voltage = np.array([[voltage_1]])
    motor_2.outputs.omega  = np.array([[motor.angular_velocity]])
    
    # Run Motor Current Function 
    i, etam = motor_2.current(conditions) 
    current_2 = i[0][0]
    

    fig_2 = plt.figure('Current_Function')
    axis_2_1 = fig_2.add_subplot(1,2,1)
    axis_2_2 = fig_2.add_subplot(1,2,2) 
    axis_2_1.plot(current_2,voltage_1)
    axis_2_2.plot(etam,voltage_1) 
        
    return
    

   
#------------------------------------
# Motor Torque Function 
#------------------------------------   
def motor_torque_function(motor,conditions):    
    # create copy of motor to test functions  
    motor_3  = motor      
    
    # Define function specific inputs  
    voltage_1 =  400
    motor_3.inputs.voltage = np.array([[voltage_1]]) 
    motor_3.inputs.omega  = np.array([[motor.angular_velocity]])   
    
    # Run Motor Torque Function 
    motor_3.torque(conditions)
    torque_3 = motor_3.outputs.torque[0][0]  


    fig_3 = plt.figure('Torque_Function')
    axis_3_1 = fig_3.add_subplot(1,2,1)
    axis_3_2 = fig_3.add_subplot(1,2,2) 
    axis_3_1.plot(torque_3,voltage_1)  
    return
    
#------------------------------------
# Motor Torque Function 
#------------------------------------   
def motor_torque_function(motor,conditions):  
    # create copy of motor to test functions   
    motor_4  = motor     
    
    # Define function specific inputs    
    voltage_1 =  400
    motor_4.inputs.voltage  = np.array([[voltage_1]]) 
    motor_4.inputs.rotor_CP = conditions.propulsion.propeller_power_coefficient
    
    # Run Motor Omega Function   
    torque_1 = motor_4.outputs.torque[0][0] 
    motor_4.inputs.torque = np.array([[torque_1]])
    
    # Run Motor Voltage-Current Function 
    motor_4.voltage_current(conditions) 
    voltage_4 = motor_4.outputs.voltage[0][0]
    current_4 = motor_4.outputs.current[0][0]  
    

    fig_4= plt.figure('CTorque_Function')
    axis_4_1 = fig_4.add_subplot(1,2,1)
    axis_4_2 = fig_4.add_subplot(1,2,2) 
    axis_4_1.plot(current_4,voltage_1)
    axis_4_2.plot(voltage_4,voltage_1) 
        
     
    return

# ----------------------------------------------------------------------        
#   Call Main
# ----------------------------------------------------------------------    

if __name__ == '__main__':
    main()
    plt.show()