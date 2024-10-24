# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core import Units ,Data  
from RCAIDE.Library.Plots  import * 
from RCAIDE.Library.Methods.Propulsors.Converters.Ducted_Fan.design_ducted_fan import design_ducted_fan
from RCAIDE.Library.Methods.Propulsors.Converters.Ducted_Fan.compute_ducted_fan_performance import compute_ducted_fan_performance
from RCAIDE.Framework.Mission.Common    import  Conditions , Results 
from RCAIDE.Library.Methods.Propulsors.Converters.DC_Motor          import design_motor 
from RCAIDE.Framework.Mission.Segments.Segment   import Segment  

# python imports  
import matplotlib.pyplot as plt                                      
import matplotlib.cm     as cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.axes3d import get_test_data
from scipy import interpolate 
import pickle    
import os
import pandas as pd 
import numpy as  np
import os                                                   
import time  
 
def main(): 
    
    ti                  = time.time()      

    tip_mach            = np.array([0.2, 0.3, 0.4, 0.5, 0.6 ])     
    mach_number         = np.array([0.01,0.1, 0.2, 0.3 , 0.4, 0.5]) 
    altitude            = np.array([0, 5000, 10000, 15000,  20000]) *Units.feet

    thrust              = np.zeros((len(altitude),len(mach_number),len(tip_mach)))
    torque              = np.zeros((len(altitude),len(mach_number),len(tip_mach)))
    power               = np.zeros((len(altitude),len(mach_number),len(tip_mach)))

    bus                                                = RCAIDE.Library.Components.Energy.Distributors.Electrical_Bus()
    bus.voltage                                        = 1000                         
    electric_ducted_fan                                = RCAIDE.Library.Components.Propulsors.Electric_Ducted_Fan()
    DF =  define_ducted_fan()
    electric_ducted_fan.ducted_fan                     = DF
    electric_ducted_fan.motor                          = define_ducted_fan_motor(bus,electric_ducted_fan.ducted_fan)
    electric_ducted_fan.electronic_speed_controller    = define_electronic_speed_controller()
    bus.propulsors.append(electric_ducted_fan)    
    
    for i in range(len(altitude)): 
        for j in range(len(mach_number)):
            for k in  range(len(tip_mach)):
                
                atmosphere                                        = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
                atmo_data                                         = atmosphere.compute_values(altitude[i])
                                                                   
                ctrl_pts =  1
                AoA =  0
                true_course      = 0
                fligth_path_angle =  0
                 
                
                segment                                                = RCAIDE.Framework.Mission.Segments.Segment() 
                
                conditions                                             = Results()  
                conditions.aerodynamics.angle_of_attack                = np.ones((ctrl_pts,1)) * AoA
                conditions.freestream.density                          = atmo_data.density  
                conditions.freestream.dynamic_viscosity                = atmo_data.dynamic_viscosity
                conditions.freestream.speed_of_sound                   = atmo_data.speed_of_sound
                conditions.freestream.temperature                      = atmo_data.temperature
                conditions.freestream.altitude                         = np.ones((ctrl_pts,1)) *altitude[i]
                conditions.frames.inertial.velocity_vector             = np.array([[mach_number[j] *atmo_data.speed_of_sound[0,0] , 0. ,0.]]) 
                conditions.frames.planet.true_course                   = np.zeros((ctrl_pts,3,3)) 
                conditions.frames.planet.true_course[:,0,0]            = np.cos(true_course),
                conditions.frames.planet.true_course[:,0,1]            = - np.sin(true_course)
                conditions.frames.planet.true_course[:,1,0]            = np.sin(true_course)
                conditions.frames.planet.true_course[:,1,1]            = np.cos(true_course) 
                conditions.frames.planet.true_course[:,2,2]            = 1 
                conditions.frames.wind.transform_to_inertial           = np.zeros((ctrl_pts,3,3))   
                conditions.frames.wind.transform_to_inertial[:,0,0]    = np.cos(fligth_path_angle) 
                conditions.frames.wind.transform_to_inertial[:,0,2]    = np.sin(fligth_path_angle) 
                conditions.frames.wind.transform_to_inertial[:,1,1]    = 1 
                conditions.frames.wind.transform_to_inertial[:,2,0]    = -np.sin(fligth_path_angle) 
                conditions.frames.wind.transform_to_inertial[:,2,2]    = np.cos(fligth_path_angle)  
                conditions.frames.body.transform_to_inertial           = np.zeros((ctrl_pts,3,3))
                conditions.frames.body.transform_to_inertial[:,0,0]    = np.cos(AoA)
                conditions.frames.body.transform_to_inertial[:,0,2]    = np.sin(AoA)
                conditions.frames.body.transform_to_inertial[:,1,1]    = 1
                conditions.frames.body.transform_to_inertial[:,2,0]    = -np.sin(AoA)
                conditions.frames.body.transform_to_inertial[:,2,2]    = np.cos(AoA)     
                segment.state.conditions                               = conditions  
                segment.state.conditions.energy[bus.tag]               = Conditions() 
                
                electric_ducted_fan.append_operating_conditions(segment,bus)     
                for tag, item in  electric_ducted_fan.items(): 
                    if issubclass(type(item), RCAIDE.Library.Components.Component):
                        item.append_operating_conditions(segment,bus,electric_ducted_fan) 
                
                # set throttle
                segment.state.conditions.energy[bus.tag][electric_ducted_fan.tag].throttle[:,0] = 1   
            
                # Run BEMT
                segment.state.conditions.expand_rows(ctrl_pts)
                ducted_fan_conditions             = segment.state.conditions.energy[bus.tag][electric_ducted_fan.tag][electric_ducted_fan.ducted_fan.tag]     
                ducted_fan_conditions.omega[:,0]  = ((tip_mach[k]*atmo_data.speed_of_sound[0,0]) /DF.tip_radius) 
                compute_ducted_fan_performance(electric_ducted_fan,segment.state,bus)
                
                thrust[i, j, k] = segment.state.conditions.energy[bus.tag][electric_ducted_fan.tag][electric_ducted_fan.ducted_fan.tag].thrust[0, 0] 
                torque[i, j, k] = segment.state.conditions.energy[bus.tag][electric_ducted_fan.tag][electric_ducted_fan.ducted_fan.tag].torque[0, 0]
                power[i, j, k]  = segment.state.conditions.energy[bus.tag][electric_ducted_fan.tag][electric_ducted_fan.ducted_fan.tag].power[0, 0]
        
    plot_results(altitude,mach_number,tip_mach,thrust,torque,power)
    
    tf                      = time.time()                                           # [s]       Define the final simulation time
    elapsed_time            = round((tf-ti),2)                                      # [s]       Compute the total simulation time

    print('Simulation Time: ' + str(elapsed_time) + ' seconds per timestep')        # [-]       Print the value of total simulation time    
    
    return

def plot_results(altitude,mach_number,tip_mach,thrust,torque,power): 
    #fig    =  plt.figure('Thrust')
    #fig.set_size_inches(7, 6)
    #axis_1 = fig.add_subplot(2,2,1) 
    #axis_2 = fig.add_subplot(2,2,2) 
    #axis_3 = fig.add_subplot(2,2,3) 
    #axis_4 = fig.add_subplot(2,2,4)
    
    #x, y = np.meshgrid(tip_mach, mach_number)
    


    #thrust_levels   = np.linspace(0,np.max(thrust),20)
    
    #thrust_cmap     = plt.get_cmap('turbo') 
    
    
    #contour_1 =  axis_1.contourf(x,y, thrust[0, :, :], thrust_levels,cmap = thrust_cmap)   
    #contour_2 =  axis_2.contourf(x,y, thrust[1, :, :], thrust_levels,cmap = thrust_cmap)   
    #contour_3 =  axis_3.contourf(x,y, thrust[2, :, :], thrust_levels,cmap = thrust_cmap)   
    #contour_4 =  axis_4.contourf(x,y, thrust[3, :, :], thrust_levels,cmap = thrust_cmap)   
    #axis_1.set_ylabel('Mach')
    #axis_2.set_ylabel('Mach')
    #axis_3.set_ylabel('Mach')
    #axis_4.set_ylabel('Mach')
    #axis_1.set_xlabel('Tip Mach')
    #axis_2.set_xlabel('Tip Mach')
    #axis_3.set_xlabel('Tip Mach')
    #axis_4.set_xlabel('Tip Mach')  
    #axis_1.set_title('Alt = 0 ft')
    #axis_2.set_title('Alt = 5000 ft')
    #axis_3.set_title('Alt = 10000 ft')
    #axis_4.set_title('Alt = 15000 ft') 
    
    #cbar  = fig.colorbar(contour_1, ax = axis_1) 
    #cbar  = fig.colorbar(contour_2, ax = axis_2) 
    #cbar  = fig.colorbar(contour_3, ax = axis_3) 
    #cbar  = fig.colorbar(contour_4, ax = axis_4) 
    #cbar.ax.set_ylabel('Thrust', rotation =  90)
    
    #fig.tight_layout()
    
    # ---------------------------------------------
    
    #fig    =  plt.figure('Power')
    #fig.set_size_inches(7, 6)
    #axis_1 = fig.add_subplot(2,2,1) 
    #axis_2 = fig.add_subplot(2,2,2) 
    #axis_3 = fig.add_subplot(2,2,3) 
    #axis_4 = fig.add_subplot(2,2,4)
    
    #x, y = np.meshgrid(tip_mach, mach_number)
    


    #thrust_levels   = np.linspace(0,np.max(power),20)
    
    #thrust_cmap     = plt.get_cmap('turbo') 
    
    
    #contour_1 =  axis_1.contourf(x,y, power[0, :, :], thrust_levels,cmap = thrust_cmap)   
    #contour_2 =  axis_2.contourf(x,y, power[1, :, :], thrust_levels,cmap = thrust_cmap)   
    #contour_3 =  axis_3.contourf(x,y, power[2, :, :], thrust_levels,cmap = thrust_cmap)   
    #contour_4 =  axis_4.contourf(x,y, power[3, :, :], thrust_levels,cmap = thrust_cmap)   
    #axis_1.set_ylabel('Mach')
    #axis_2.set_ylabel('Mach')
    #axis_3.set_ylabel('Mach')
    #axis_4.set_ylabel('Mach')
    #axis_1.set_xlabel('Tip Mach')
    #axis_2.set_xlabel('Tip Mach')
    #axis_3.set_xlabel('Tip Mach')
    #axis_4.set_xlabel('Tip Mach')  
    #axis_1.set_title('Alt = 0 ft')
    #axis_2.set_title('Alt = 5000 ft')
    #axis_3.set_title('Alt = 10000 ft')
    #axis_4.set_title('Alt = 15000 ft') 
    
    #cbar  = fig.colorbar(contour_1, ax = axis_1) 
    #cbar  = fig.colorbar(contour_2, ax = axis_2) 
    #cbar  = fig.colorbar(contour_3, ax = axis_3) 
    #cbar  = fig.colorbar(contour_4, ax = axis_4) 
    #cbar.ax.set_ylabel('Power', rotation =  90)
    
    #fig.tight_layout()    
    
    # ---------------------------------------------
    
    fig    =  plt.figure('Torque')
    fig.set_size_inches(7, 6)
    axis_1 = fig.add_subplot(2,2,1) 
    axis_2 = fig.add_subplot(2,2,2) 
    axis_3 = fig.add_subplot(2,2,3) 
    axis_4 = fig.add_subplot(2,2,4)
    
    x, y = np.meshgrid(tip_mach, mach_number)
    


    thrust_levels   = np.linspace(0,np.max(torque),20)
    
    thrust_cmap     = plt.get_cmap('turbo') 
    
    
    contour_1 =  axis_1.contourf(x,y, torque[0, :, :], thrust_levels,cmap = thrust_cmap)   
    contour_2 =  axis_2.contourf(x,y, torque[1, :, :], thrust_levels,cmap = thrust_cmap)   
    contour_3 =  axis_3.contourf(x,y, torque[2, :, :], thrust_levels,cmap = thrust_cmap)   
    contour_4 =  axis_4.contourf(x,y, torque[3, :, :], thrust_levels,cmap = thrust_cmap)   
    axis_1.set_ylabel('Mach')
    axis_2.set_ylabel('Mach')
    axis_3.set_ylabel('Mach')
    axis_4.set_ylabel('Mach')
    axis_1.set_xlabel('Tip Mach')
    axis_2.set_xlabel('Tip Mach')
    axis_3.set_xlabel('Tip Mach')
    axis_4.set_xlabel('Tip Mach')  
    axis_1.set_title('Alt = 0 ft')
    axis_2.set_title('Alt = 5000 ft')
    axis_3.set_title('Alt = 10000 ft')
    axis_4.set_title('Alt = 15000 ft') 
    
    cbar  = fig.colorbar(contour_1, ax = axis_1) 
    cbar  = fig.colorbar(contour_2, ax = axis_2) 
    cbar  = fig.colorbar(contour_3, ax = axis_3) 
    cbar  = fig.colorbar(contour_4, ax = axis_4) 
    cbar.ax.set_ylabel('Torque', rotation =  90)
    
    fig.tight_layout()        
    return

def plot_style(number_of_lines= 10): 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 20,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14,
                  'axes.titlesize': 18,
                  #figure.dpi': 1200
                  }

    # Universal Plot Settings  
    plt.rcParams.update(parameters)
    plot_parameters                        = Data()
    plot_parameters.line_width             = 1.5  
    plot_parameters.line_style             = ['-','--']
    plot_parameters.marker_size            = 4
    plot_parameters.legend_fontsize        = '12'
    plot_parameters.legend_title_font_size = 14
    plot_parameters.axis_font_size         = 16
    plot_parameters.title_font_size        = 16   
    plot_parameters.markers                =  ['o','x','o','v','P','p','^','D','*']
    plot_parameters.color                  = cm.inferno(np.linspace(0,0.9,number_of_lines)) 

    return plot_parameters
 
def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2
 

# ----------------------------------------------------------------------
#   Save Results
# ----------------------------------------------------------------------
def save_results(results,filename): 
    pickle_file  =  filename + '.pkl'
    with open(pickle_file, 'wb') as file:
        pickle.dump(results, file) 
    return   

# ------------------------------------------------------------------
#   Load Results
# ------------------------------------------------------------------   
def load_results(filename):  
    load_file = filename + '.pkl' 
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results  


def define_ducted_fan():
    ospath                                = os.path.abspath(__file__)
    separator                             = os.path.sep
    rel_path                              = ospath.split('Ducted_Fan_Verification.py')[0]   
    
    
    ducted_fan                                   = RCAIDE.Library.Components.Propulsors.Converters.Ducted_Fan()
    ducted_fan.tag                               = 'test_ducted_fan'
    ducted_fan.number_of_rotor_blades            = 12 #22 
    ducted_fan.number_of_radial_stations         = 20
    ducted_fan.tip_radius                        = 3.124 / 2
    ducted_fan.hub_radius                        = 3.124 /2 * 0.25
    ducted_fan.blade_clearance                   = 0.01
    ducted_fan.length                            = 2
    ducted_fan.rotor_percent_x_location          = 0.4
    ducted_fan.stator_percent_x_location         = 0.7
    ducted_fan.cruise.design_thrust              = 10000  
    ducted_fan.cruise.design_altitude            = 20000 *Units.ft  
    ducted_fan.cruise.design_tip_mach            = 0.7
    speed_of_sound                               =  316.032
    ducted_fan.cruise.design_angular_velocity    = (ducted_fan.cruise.design_tip_mach *speed_of_sound) /ducted_fan.tip_radius  # 1352 RPM
    ducted_fan.cruise.design_freestream_velocity = 0.45* speed_of_sound
    ducted_fan.cruise.design_reference_velocity  = 0.45*  speed_of_sound 
    airfoil                                      = RCAIDE.Library.Components.Airfoils.NACA_4_Series_Airfoil() 
    airfoil.NACA_4_Series_code                   = '2208'
    ducted_fan.append_duct_airfoil(airfoil)  
    airfoil                                      = RCAIDE.Library.Components.Airfoils.NACA_4_Series_Airfoil()
    airfoil.NACA_4_Series_code                   = '0021'    
    ducted_fan.append_hub_airfoil(airfoil)    
    design_ducted_fan(ducted_fan) 
    
    return ducted_fan 

def define_ducted_fan_motor(bus, ducted_fan):

    motor                                            = RCAIDE.Library.Components.Propulsors.Converters.DC_Motor()
    motor.efficiency                                 = 0.98
    motor.origin                                     = [[2.,  2.5, 0.95]]
    motor.nominal_voltage                            = bus.voltage
    motor.no_load_current                            = 1
    motor.rotor_radius                               = ducted_fan.tip_radius
    motor.design_torque                              = ducted_fan.cruise.design_torque
    motor.angular_velocity                           = ducted_fan.cruise.design_angular_velocity 
    design_motor(motor)      
    return motor

def define_electronic_speed_controller():

    # Electronic Speed Controller       
    esc                                              = RCAIDE.Library.Components.Energy.Modulators.Electronic_Speed_Controller()
    esc.tag                                          = 'esc'
    esc.efficiency                                   = 0.95        
    return esc


    
if __name__ == '__main__': 
    main()
    plt.show() 