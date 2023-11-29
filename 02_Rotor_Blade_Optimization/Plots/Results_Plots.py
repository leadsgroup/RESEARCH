import MARC 
from MARC.Core import Units , Data, to_numpy
# Package Imports 
import matplotlib.cm as cm 
import numpy as np 
import jax 
import jax.numpy as jnp  
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker 
import matplotlib.colors as colors  
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d import Axes3D  
from MARC.Analyses.Mission.Segments.Segment                           import Segment 
from MARC.Methods.Noise.Fidelity_Zero.Propeller.propeller_mid_fidelity import propeller_mid_fidelity
from MARC.Analyses.Mission.Segments.Conditions.Aerodynamics           import Aerodynamics 
from MARC.Components.Energy.Networks.Battery_Electric_Rotor                import Battery_Electric_Rotor  
from MARC.Components.Energy.Converters                                import Lift_Rotor , Prop_Rotor, Propeller
from MARC.Methods.Geometry.Two_Dimensional.Cross_Section.Airfoil.generate_interpolated_airfoils import generate_interpolated_airfoils 

from MARC.Components.Energy.Networks.Battery_Electric_Rotor                import Battery_Electric_Rotor 
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from MARC.Methods.Geometry.Two_Dimensional.Cross_Section.Airfoil.import_airfoil_geometry import import_airfoil_geometry
from MARC.Methods.Geometry.Two_Dimensional.Cross_Section.Airfoil.compute_naca_4series    import compute_naca_4series
import os
import pickle
from scipy.interpolate import interp1d



def main(): 
    lift_rotor_plots()
    
    prop_rotor_plots()
    
    #traditional_new_parameterization_comparison()
    
    #traditional_new_parameterization_3d_comparison()
    #prop_rotor_traditional_new_parameterization_3d_comparison()
    return 


def lift_rotor_plots():
    
 
    rotor                            = Lift_Rotor() 
    rotor.tag                        = 'rotor'
    rotor.orientation_euler_angles   = [0, 90*Units.degrees,0]
    rotor.tip_radius                 = 2.7/2
    rotor.hub_radius                 = 0.15 * rotor.tip_radius  
    rotor.number_of_blades           = 3    
    
    # hover 
    rotor.hover = Data()
    rotor.hover.design_thrust               = 23544/(8)  # based on Stopped-Rotor V2, vehicle weight/(number of rotors - 1 )
    rotor.hover.design_freestream_velocity  = np.sqrt(rotor.hover.design_thrust/(2*1.2*np.pi*(rotor.tip_radius**2))) # Ideal power  
    rotor.hover.design_altitude             = 20 * Units.feet   
    
    # OEI 
    rotor.OEI = Data()
    rotor.OEI.design_thrust               = 23544/(6)  # based on Stopped-Rotor V2, vehicle weight/(number of rotors - 1 )
    rotor.OEI.design_freestream_velocity  = np.sqrt(rotor.OEI.design_thrust/(2*1.2*np.pi*(rotor.tip_radius**2))) # Ideal power  
    rotor.OEI.design_altitude             = 20 * Units.feet   
       
    
    airfoil                          = MARC.Components.Airfoils.Airfoil()    
    airfoil.coordinate_file          =  '..' + separator + 'Airfoils' + separator + 'NACA_4412.txt'
    airfoil.polar_files              = ['..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_50000.txt',
                                         '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_100000.txt',
                                         '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_200000.txt',
                                         '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_500000.txt',
                                         '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_1000000.txt']
    #rotor.append_airfoil(airfoil)   
    #rotor.airfoil_polar_stations          = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]         
    
    
    '''lift rotor design using new method '''
    alpha_weights                      = np.linspace(0,1.0,1001)  
    alpha_weights[-1]                  = 0.99
    save_figures                       = True   
    folder_name                        = '../Lift_Rotor_Optmization/'
    prop_rotor_flag = False 
    plot_rotor_blade_comparisons(rotor,prop_rotor_flag,folder_name,alpha_weights =alpha_weights,beta_weights = None,gamma_weights = None,add_plot_legends = False , save_figures = save_figures)    
        
    return 


def prop_rotor_plots():
    

    
    # DEFINE ROTOR OPERATING CONDITIONS 

    prop_rotor                                    = MARC.Components.Energy.Converters.Prop_Rotor() 
    prop_rotor.tag                                = 'prop_rotor'     
    prop_rotor.tip_radius                         = 3/2
    prop_rotor.hub_radius                         = 0.15 * prop_rotor.tip_radius
    prop_rotor.number_of_blades                   = 3  
    
    # HOVER 
    prop_rotor.hover = Data()
    prop_rotor.hover.design_altitude              = 20 * Units.feet                  
    prop_rotor.hover.design_thrust                = 23175.5364/6 # weight of joby-like aircrft
    prop_rotor.hover.design_freestream_velocity   = np.sqrt(prop_rotor.hover.design_thrust/(2*1.2*np.pi*(prop_rotor.tip_radius**2)))  
  
    # OEI 
    prop_rotor.OEI = Data()
    prop_rotor.OEI.design_thrust                  = 23175.5364/4  # based on Stopped-Rotor V2, vehicle weight/(number of rotors - 1 )
    prop_rotor.OEI.design_freestream_velocity     = np.sqrt(prop_rotor.OEI.design_thrust/(2*1.2*np.pi*(prop_rotor.tip_radius**2))) # Ideal power  
    prop_rotor.OEI.design_altitude                = 20 * Units.feet   
    
    # CRUISE                 
    prop_rotor.cruise = Data()  
    prop_rotor.cruise.design_altitude             = 2500 * Units.feet                      
    prop_rotor.cruise.design_thrust               = 4000/6
    prop_rotor.cruise.design_freestream_velocity  = 175*Units.mph 


    airfoil                                       = MARC.Components.Airfoils.Airfoil()    
    airfoil.coordinate_file                       =  '..' + separator + 'Airfoils' + separator + 'NACA_4412.txt'
    airfoil.polar_files                           = ['..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_50000.txt',
                                                      '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_100000.txt',
                                                      '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_200000.txt',
                                                      '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_500000.txt',
                                                      '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_1000000.txt'] 
    #prop_rotor.append_airfoil(airfoil)   
    prop_rotor.airfoil_polar_stations             = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]   
     
    '''prop rotor design using new method ''' 
    ## ----------------------
    ##  PLOTS 
    ## ---------------------- 
    alpha_weights                      = np.array([0.8,0.8])# ,0.9,0.2
    beta_weights                       = np.array([0.1,0.9 ]) # ,0.5,0.5   
    save_figures                       = True  
    prop_rotor_flag = True 
    folder_name                        = '../Prop_Rotor_Optimization/'
    plot_rotor_blade_comparisons(prop_rotor,prop_rotor_flag,folder_name,alpha_weights =alpha_weights,beta_weights = beta_weights,gamma_weights = None,add_plot_legends = False , save_figures = save_figures)    
         
    
    return 

def traditional_new_parameterization_3d_comparison():
    

    rotor                            = Lift_Rotor() 
    rotor.tag                        = 'rotor'
    rotor.orientation_euler_angles   = [0, 90*Units.degrees,0]
    rotor.tip_radius                 = 2.7/2
    rotor.hub_radius                 = 0.15 * rotor.tip_radius  
    rotor.number_of_blades           = 3    
    
    # hover 
    rotor.hover = Data()
    rotor.hover.design_thrust               = 23544/(8)  # based on Stopped-Rotor V2, vehicle weight/(number of rotors - 1 )
    rotor.hover.design_freestream_velocity  = np.sqrt(rotor.hover.design_thrust/(2*1.2*np.pi*(rotor.tip_radius**2))) # Ideal power  
    rotor.hover.design_altitude             = 20 * Units.feet   
    
    # OEI 
    rotor.OEI = Data()
    rotor.OEI.design_thrust               = 23544/(6)  # based on Stopped-Rotor V2, vehicle weight/(number of rotors - 1 )
    rotor.OEI.design_freestream_velocity  = np.sqrt(rotor.OEI.design_thrust/(2*1.2*np.pi*(rotor.tip_radius**2))) # Ideal power  
    rotor.OEI.design_altitude             = 20 * Units.feet   
       
     
    
    
    '''lift rotor design using new method '''
    alpha_weights                      = np.array([ 0.5,0.5 ]) 
    save_figures                       = True     
    folder_names     = ['../Lift_Rotor_Optmization/','../Lift_Rotor_Optmization_trad/']
    rotor_labels     = ['N.P.','T.P.']  
    fig_name         = ['New_Paramterization','Traditional_Paramterization' ]
    add_plot_legends = True 
    save_figures     = True   
   
    PP              = define_plot_parameters() 
    PP.colors_1     = cm.viridis(np.linspace(0,1,len(alpha_weights)))  
                
          
    
    for i in range(len(alpha_weights)):    
        alpha_opt_weight = str(format(alpha_weights[i],'.5f'))
        alpha_opt_weight = alpha_opt_weight.replace('.','_')     
        rotor_file_name  = 'LR_Alpha_' + alpha_opt_weight   
        rotor       = load_blade_geometry(folder_names[i],rotor_file_name) 
        
    
        airfoil                          = MARC.Components.Airfoils.Airfoil()    
        airfoil.coordinate_file          =  'Airfoils/NACA_4412.txt'
        airfoil.polar_files              = ['Airfoils/Polars/NACA_4412_polar_Re_50000.txt',
                                             'Airfoils/Polars/NACA_4412_polar_Re_100000.txt',
                                             'Airfoils/Polars/NACA_4412_polar_Re_200000.txt',
                                             'Airfoils/Polars/NACA_4412_polar_Re_500000.txt',
                                             'Airfoils/Polars/NACA_4412_polar_Re_1000000.txt']
        rotor.append_airfoil(airfoil)   
        rotor.airfoil_polar_stations          = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]        
        
        plot_3d_rotor_geometry(rotor,fig_name[i],save_figures, elevation_angle = 45,  aximuth_angle = 0, cpt=0,rotor_face_color='dodgerblue',
                        rotor_edge_color='darkblue',rotor_alpha=1.0)  

    return 


def prop_rotor_traditional_new_parameterization_3d_comparison():
    
    

    # DEFINE ROTOR OPERATING CONDITIONS 

    prop_rotor                                    = MARC.Components.Energy.Converters.Prop_Rotor() 
    prop_rotor.tag                                = 'prop_rotor'     
    prop_rotor.tip_radius                         = 3/2
    prop_rotor.hub_radius                         = 0.15 * prop_rotor.tip_radius
    prop_rotor.number_of_blades                   = 3  
    
    # HOVER 
    prop_rotor.hover = Data()
    prop_rotor.hover.design_altitude              = 20 * Units.feet                  
    prop_rotor.hover.design_thrust                = 23175.5364/6 # weight of joby-like aircrft
    prop_rotor.hover.design_freestream_velocity   = np.sqrt(prop_rotor.hover.design_thrust/(2*1.2*np.pi*(prop_rotor.tip_radius**2)))  
  
    # OEI 
    prop_rotor.OEI = Data()
    prop_rotor.OEI.design_thrust                  = 23175.5364/4  # based on Stopped-Rotor V2, vehicle weight/(number of rotors - 1 )
    prop_rotor.OEI.design_freestream_velocity     = np.sqrt(prop_rotor.OEI.design_thrust/(2*1.2*np.pi*(prop_rotor.tip_radius**2))) # Ideal power  
    prop_rotor.OEI.design_altitude                = 20 * Units.feet   
    
    # CRUISE                 
    prop_rotor.cruise = Data()  
    prop_rotor.cruise.design_altitude             = 2500 * Units.feet                      
    prop_rotor.cruise.design_thrust               = 4000/6
    prop_rotor.cruise.design_freestream_velocity  = 175*Units.mph 


    airfoil                                       = MARC.Components.Airfoils.Airfoil()    
    airfoil.coordinate_file                       =  '..' + separator + 'Airfoils' + separator + 'NACA_4412.txt'
    airfoil.polar_files                           = ['..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_50000.txt',
                                                      '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_100000.txt',
                                                      '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_200000.txt',
                                                      '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_500000.txt',
                                                      '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_1000000.txt'] 
    #prop_rotor.append_airfoil(airfoil)   
    prop_rotor.airfoil_polar_stations             = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]   
    
 
    '''lift rotor design using new method '''

    alpha_weights   = np.array([0.8,0.8,0.9,0.2]) 
    beta_weights    = np.array([0.1,0.9,0.5,0.5])     
    save_figures    = True  
    prop_rotor_flag = True 
    folder_name     = '../Prop_Rotor_Optimization/' 
    save_figures    = True   
   
    PP              = define_plot_parameters() 
    PP.colors_1     = cm.viridis(np.linspace(0,1,len(alpha_weights)))  
                
          
    
    for i in range(len(alpha_weights)):    

        alpha_opt_weight_1 = str(format(alpha_weights[i],'.1f'))
        alpha_opt_weight_1 = alpha_opt_weight_1.replace('.','_')    
        beta_opt_weight_1  = str(format(beta_weights[i],'.1f'))
        beta_opt_weight_1  = beta_opt_weight_1.replace('.','_')   
        
        save_file_name     = 'PR_A_' + alpha_opt_weight_1 + '_B_' + beta_opt_weight_1 
        
        
        alpha_opt_weight = str(format(alpha_weights[i],'.5f'))
        alpha_opt_weight = alpha_opt_weight.replace('.','_')    
        beta_opt_weight  = str(format(beta_weights[i],'.5f'))
        beta_opt_weight  = beta_opt_weight.replace('.','_')    
        rotor_file_name  = 'PR_Alpha_' + alpha_opt_weight + '_Beta_' + beta_opt_weight          
        rotor              = load_blade_geometry(folder_name,rotor_file_name)  
    
        airfoil                          = MARC.Components.Airfoils.Airfoil()    
        airfoil.coordinate_file          =  'Airfoils/NACA_4412.txt'
        airfoil.polar_files              = ['Airfoils/Polars/NACA_4412_polar_Re_50000.txt',
                                             'Airfoils/Polars/NACA_4412_polar_Re_100000.txt',
                                             'Airfoils/Polars/NACA_4412_polar_Re_200000.txt',
                                             'Airfoils/Polars/NACA_4412_polar_Re_500000.txt',
                                             'Airfoils/Polars/NACA_4412_polar_Re_1000000.txt']
        rotor.append_airfoil(airfoil)   
        rotor.airfoil_polar_stations          = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]        
        
        plot_3d_rotor_geometry(rotor,save_file_name,save_figures, elevation_angle = 45,  aximuth_angle = 0, cpt=0,rotor_face_color='dodgerblue',
                        rotor_edge_color='blue',rotor_alpha=1.0)  

    return 

def traditional_new_parameterization_comparison():
    

    rotor                            = Lift_Rotor() 
    rotor.tag                        = 'rotor'
    rotor.orientation_euler_angles   = [0, 90*Units.degrees,0]
    rotor.tip_radius                 = 2.7/2
    rotor.hub_radius                 = 0.15 * rotor.tip_radius  
    rotor.number_of_blades           = 3    
    
    # hover 
    rotor.hover = Data()
    rotor.hover.design_thrust               = 23544/(8)  # based on Stopped-Rotor V2, vehicle weight/(number of rotors - 1 )
    rotor.hover.design_freestream_velocity  = np.sqrt(rotor.hover.design_thrust/(2*1.2*np.pi*(rotor.tip_radius**2))) # Ideal power  
    rotor.hover.design_altitude             = 20 * Units.feet   
    
    # OEI 
    rotor.OEI = Data()
    rotor.OEI.design_thrust               = 23544/(6)  # based on Stopped-Rotor V2, vehicle weight/(number of rotors - 1 )
    rotor.OEI.design_freestream_velocity  = np.sqrt(rotor.OEI.design_thrust/(2*1.2*np.pi*(rotor.tip_radius**2))) # Ideal power  
    rotor.OEI.design_altitude             = 20 * Units.feet   
       
    
    airfoil                          = MARC.Components.Airfoils.Airfoil()    
    ospath                           = os.path.abspath(__file__)
    separator                        = os.path.sep
    rel_path                         = os.path.dirname(ospath) + separator  
    airfoil.coordinate_file          =  '..' + separator + 'Airfoils' + separator + 'NACA_4412.txt'
    airfoil.polar_files              = ['..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_50000.txt',
                                         '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_100000.txt',
                                         '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_200000.txt',
                                         '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_500000.txt',
                                         '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_1000000.txt']
    #rotor.append_airfoil(airfoil)   
    #rotor.airfoil_polar_stations          = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]         
    
    
    '''lift rotor design using new method '''
    alpha_weights                      = np.array([ 0.5,0.5 ]) 
    save_figures                       = True     
    folder_names     = ['../Lift_Rotor_Optmization/','../Lift_Rotor_Optmization_trad/']
    rotor_labels     = ['N.P.','T.P.']  
    fig_name         = ['New_Paramterization','Traditional_Paramterization' ]
    add_plot_legends = True 
    save_figures     = True   
   
    PP              = define_plot_parameters() 
    PP.colors_1     = cm.viridis(np.linspace(0,1,len(alpha_weights)))  
                
         
    AXES , FIGURES  = set_up_axes(PP,'LR') 
    
    # file names 
    fig_1_name = "Method_Twist_Comparison" 
    fig_2_name = "Method_Chord_Comparison"  
    fig_3_name = "Method_Thickness_Comparison"   
    fig_4_name = 'Method_Thrust_Comparison' 
    fig_5_name = 'Method_Torque_Comparison'   
    fig_6_name = 'Method_Total_SPL_Comparison'   
    fig_7_name = "Method_Power_Noise_Pareto"  
    fig_8_name = 'Method_Power_RPM_Pareto'  
    
    for i in range(len(alpha_weights)):    
        alpha_opt_weight = str(format(alpha_weights[i],'.5f'))
        alpha_opt_weight = alpha_opt_weight.replace('.','_')     
        rotor_file_name  = 'LR_Alpha_' + alpha_opt_weight   
        rotor       = load_blade_geometry(folder_names[i],rotor_file_name)
        plot_rotor_2d_geoemtry_and_performance(rotor,AXES,PP.colors_1[i],PP.line_styles[0],
                                                   PP.line_styles[1],PP.markers[0],PP.markers[1],rotor_labels[i]) 
        
        plot_3d_rotor_geometry(rotor,fig_name[i],save_figures, elevation_angle = 45,  aximuth_angle = 0, cpt=0,rotor_face_color='dodgerblue',
                        rotor_edge_color='darkblue',rotor_alpha=1.0)  
         
     
    if add_plot_legends:        
        AXES[0].legend(loc='upper right')
        AXES[1].legend(loc='upper right')
        AXES[2].legend(loc='upper right')
        AXES[3].legend(loc='upper right')
        AXES[4].legend(loc='upper right')
        AXES[5].legend(loc='upper right')
        AXES[6].legend(loc='upper right')
        AXES[7].legend(loc='upper right')
        
    
    # axis limits 
    ymin0, ymax0 = AXES[0].get_ylim()
    ymin1, ymax1 = AXES[1].get_ylim()
    ymin2, ymax2 = AXES[2].get_ylim()
    ymin3, ymax3 = AXES[3].get_ylim()
    ymin4, ymax4 = AXES[4].get_ylim()
    ymin5, ymax5 = AXES[5].get_ylim() 
    ymin6, ymax6 = AXES[6].get_ylim()
    ymin7, ymax7 = AXES[7].get_ylim()
    
    xmin6, xmax6 = AXES[6].get_xlim()
    xmin7, xmax7 = AXES[7].get_xlim() 
    
    
    AXES[0].set_ylim([ymin0, ymax0])    
    AXES[1].set_ylim([ymin1, ymax1])    
    AXES[2].set_ylim([ymin2, ymax2])    
    AXES[3].set_ylim([ymin3, ymax3])  
    AXES[4].set_ylim([ymin4, ymax4])     
    AXES[5].set_ylim([ymin5, ymax5]) 
    AXES[6].set_ylim([ymin6, ymax6])     
    AXES[7].set_ylim([ymin7, ymax7]) 

    AXES[6].set_xlim([xmin6, xmax6])     
    AXES[7].set_xlim([xmin7, xmax7])    
    
    label_location = 0.5

    # axis labels 
    AXES[0].text(0.02,  (ymax0+ymin0)*label_location, 'Root', fontsize=25)    
    AXES[1].text(0.02,  (ymax1+ymin1)*label_location, 'Root', fontsize=25)  
    AXES[2].text(0.02,  (ymax2+ymin2)*label_location, 'Root', fontsize=25)   
    AXES[3].text(0.02,  (ymax3+ymin3)*label_location, 'Root', fontsize=25) 
    AXES[4].text(0.02 , (ymax4+ymin4)*label_location, 'Root', fontsize=25)     
                            
    # Assign Colormap and save plots 
    cmap     = plt.get_cmap("viridis")
    new_cmap = truncate_colormap(cmap, 0.0, 1.0)
    norm     = plt.Normalize(0,1) 
    sm       =  ScalarMappable(norm=norm, cmap=new_cmap)
    ax_ticks = np.linspace(0,1,11)
    sm.set_array([])  

    cmap_2      = plt.get_cmap("viridis")
    new_cmap_2  = truncate_colormap(cmap_2, 0.0, 1.0)
    norm_2      = plt.Normalize(0.45,0.7) 
    sm_2        = ScalarMappable(norm=norm_2, cmap=new_cmap_2 )
    ax_ticks_2  = np.linspace(0.45,0.7,6)
    sm_2.set_array([])     
            
    sfmt = ticker.ScalarFormatter(useMathText=True) 
    sfmt = ticker.FormatStrFormatter('%.1f')     
    sfmt2 = ticker.ScalarFormatter(useMathText=True) 
    sfmt2 = ticker.FormatStrFormatter('%.2f')     
    cbar_1 = FIGURES[0].colorbar(sm,   ax = AXES[0], ticks = list(ax_ticks),  format= sfmt)
    cbar_2 = FIGURES[1].colorbar(sm,   ax = AXES[1], ticks = list(ax_ticks),  format= sfmt)
    cbar_3 = FIGURES[2].colorbar(sm,   ax = AXES[2], ticks = list(ax_ticks),  format= sfmt)
    cbar_4 = FIGURES[3].colorbar(sm,   ax = AXES[3], ticks = list(ax_ticks),  format= sfmt)
    cbar_5 = FIGURES[4].colorbar(sm,   ax = AXES[4], ticks = list(ax_ticks),  format= sfmt) 
    cbar_6 = FIGURES[5].colorbar(sm,   ax = AXES[5], ticks = list(ax_ticks),  format= sfmt)
    cbar_7 = FIGURES[6].colorbar(sm,   ax = AXES[6], ticks = list(ax_ticks),  format= sfmt)
    cbar_8 = FIGURES[7].colorbar(sm_2, ax = AXES[7], ticks = list(ax_ticks_2),  format= sfmt2)
    
    
    cbar_1.set_label(r'$\alpha$')
    cbar_2.set_label(r'$\alpha$')
    cbar_3.set_label(r'$\alpha$')
    cbar_4.set_label(r'$\alpha$') 
    cbar_5.set_label(r'$\alpha$')  
    cbar_6.set_label(r'$\alpha$')
    cbar_7.set_label(r'$\alpha$') 
    cbar_8.set_label(r'Tip Mach')   
     
    FIGURES[0].tight_layout()
    FIGURES[1].tight_layout()
    FIGURES[2].tight_layout()
    FIGURES[3].tight_layout()
    FIGURES[4].tight_layout() 
    FIGURES[5].tight_layout()
    FIGURES[6].tight_layout(rect= (0.05,0,1,1))
    FIGURES[7].tight_layout(rect= (0.05,0,1,1))  
    
    if save_figures:
        FIGURES[0].savefig(fig_1_name  + '.png')               
        FIGURES[1].savefig(fig_2_name  + '.png')               
        FIGURES[2].savefig(fig_3_name  + '.png')               
        FIGURES[3].savefig(fig_4_name  + '.png')          
        FIGURES[4].savefig(fig_5_name  + '.png')                
        FIGURES[5].savefig(fig_6_name  + '.png')               
        FIGURES[6].savefig(fig_7_name  + '.png')            
        FIGURES[7].savefig(fig_8_name  + '.png')    
    
        
    return 
# ------------------------------------------------------------------ 
# Define plot parameters 
# ------------------------------------------------------------------  
def define_plot_parameters(): 

    plt.rcParams['axes.linewidth'] = 2.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 32,
                  'legend.fontsize': 22,
                  'xtick.labelsize': 28,
                  'ytick.labelsize': 28,
                  'axes.titlesize': 32}
    plt.rcParams.update(parameters)
    plot_parameters                  = Data()
    plot_parameters.line_width       = 2
    plot_parameters.line_styles      = ['-','--','--',':','--']
    plot_parameters.figure_width     = 10
    plot_parameters.figure_height    = 7
    plot_parameters.marker_size      = 10
    plot_parameters.legend_font_size = 20
    plot_parameters.alpha_val        = 0.25
    plot_parameters.root_color       = 'grey'
    plot_parameters.plot_grid        = True   

    plot_parameters.colors           = [['black','firebrick','darkblue'],
                                        ['dimgray','red','blue'], 
                                        ['darkgray','salmon','deepskyblue']]  

    plot_parameters.colors_1         = ['black','darkmagenta','mediumblue','darkgreen','darkgoldenrod','darkred']   
    plot_parameters.colors_2         = ['grey','orchid','darkcyan','green','orange','red']        
    plot_parameters.markers          = ['s','o','v','P','p','^','D','X','*']   
    
    return plot_parameters 


# ------------------------------------------------------------------ 
# Setup Axes 
# ------------------------------------------------------------------ 
def set_up_axes(PP,rotor_tag):
    
    # ------------------------------------------------------------------
    #   Twist Distribition
    # ------------------------------------------------------------------
    fig_1_name =  rotor_tag + "_Twist_Comparison"  
    fig_1 = plt.figure(fig_1_name)
    fig_1.set_size_inches(PP.figure_width,PP.figure_height)
    axis_1 = fig_1.add_subplot(1,1,1)
    axis_1.set_ylabel(r'$\beta$ ($\degree$)') 
    axis_1.set_xlabel('r')    
    axis_1.minorticks_on()   
    
    # ------------------------------------------------------------------
    #   Chord Distribution
    # ------------------------------------------------------------------ 
    fig_2_name =   rotor_tag + "_Chord_Comparison"  
    fig_2 = plt.figure(fig_2_name)     
    fig_2.set_size_inches(PP.figure_width,PP.figure_height) 
    axis_2 = fig_2.add_subplot(1,1,1)  
    axis_2.set_ylabel('c (m)') 
    axis_2.set_xlabel('r')    
    axis_2.minorticks_on()    

    # ------------------------------------------------------------------
    #  Thickness Distribution
    # ------------------------------------------------------------------ 
    fig_3_name =  rotor_tag + "_Thickness_Comparison"  
    fig_3 = plt.figure(fig_3_name)     
    fig_3.set_size_inches(PP.figure_width,PP.figure_height) 
    axis_3 = fig_3.add_subplot(1,1,1)  
    axis_3.set_ylabel('t (m)') 
    axis_3.set_xlabel('r')    
    axis_3.minorticks_on()   

    # ------------------------------------------------------------------
    # Thrust Comparison
    # ------------------------------------------------------------------      
    fig_4_name =  rotor_tag + "Thrust_Comparison" 
    fig_4 = plt.figure(fig_4_name)    
    fig_4.set_size_inches(PP.figure_width, PP.figure_height) 
    axis_4= fig_4.add_subplot(1,1,1)    
    axis_4.set_ylabel(r"T' (N)") 
    axis_4.set_xlabel('r')  

    # ------------------------------------------------------------------
    # Torque Comparison
    # ------------------------------------------------------------------   
    fig_5_name =   rotor_tag + "_Torque_Comparison"  
    fig_5 = plt.figure(fig_5_name)            
    fig_5.set_size_inches(PP.figure_width, PP.figure_height) 
    axis_5 = fig_5.add_subplot(1,1,1)    
    axis_5.set_ylabel(r"Q' (N-m)") 
    axis_5.set_xlabel('r')   

    # ------------------------------------------------------------------
    # Noise Requency Spectrum 
    # ------------------------------------------------------------------ 
    fig_6_name =  rotor_tag + "_Total_SPL_Comparison" 
    fig_6 = plt.figure(fig_6_name)            
    fig_6.set_size_inches(PP.figure_width, PP.figure_height) 
    axis_6 = fig_6.add_subplot(1,1,1)    
    axis_6.set_xscale('log') 
    axis_6.set_ylabel(r'SPL$_{1/3}$ (dBA)')
    axis_6.set_xlabel('Frequency (Hz)') 
    axis_6.set_ylim([0,100])  

    # ------------------------------------------------------------------
    # Performance Pareto 1
    # ------------------------------------------------------------------       
    fig_7_name =  rotor_tag + "_Power_Noise_Pareto"  
    fig_7 = plt.figure(fig_7_name)     
    fig_7.set_size_inches(PP.figure_width,PP.figure_height) 
    axis_7 = fig_7.add_subplot(1,1,1)  
    axis_7.set_xlabel('Power (kW)') 
    axis_7.set_ylabel('SPL (dBA)')    
    axis_7.minorticks_on()  
 

    # ------------------------------------------------------------------
    # Performance Pareto 2
    # ------------------------------------------------------------------    
    fig_8_name =  rotor_tag + "_Power_RPM_Pareto"  
    fig_8 = plt.figure(fig_8_name)     
    fig_8.set_size_inches(PP.figure_width,PP.figure_height) 
    axis_8 = fig_8.add_subplot(1,1,1)  
    axis_8.set_xlabel('Power (kW)') 
    axis_8.set_ylabel('RPM')    
    axis_8.minorticks_on()  
    
    
    AXES    = [axis_1,axis_2,axis_3,axis_4,axis_5,axis_6,axis_7,axis_8]
    FIGURES = [fig_1,fig_2,fig_3,fig_4,fig_5,fig_6,fig_7,fig_8]
    return AXES , FIGURES 


 
# ------------------------------------------------------------------ 
# Plot lift rotor pareto fronteir 
# ------------------------------------------------------------------ 
def plot_rotor_blade_comparisons(rotor,prop_rotor,folder_name,alpha_weights = None,beta_weights = None,gamma_weights = None,add_plot_legends = False , save_figures = True):     
    PP              = define_plot_parameters() 
     
    PP.colors_1     = cm.viridis(np.linspace(0,1,len(alpha_weights)))  
  
        
    if prop_rotor: 
    
        PP.colors_1     = ['green','darkblue' ]
        PP.colors_2     = ['yellowgreen','dodgerblue']
        
        AXES , FIGURES  = set_up_axes(PP,'PR') 
        

        alpha_opt_weight_1 = str(format(alpha_weights[0],'.1f'))
        alpha_opt_weight_1 = alpha_opt_weight_1.replace('.','_')    
        beta_opt_weight_1  = str(format(beta_weights[0],'.1f'))
        beta_opt_weight_1  = beta_opt_weight_1.replace('.','_')   

        alpha_opt_weight_2 = str(format(alpha_weights[1],'.1f'))
        alpha_opt_weight_2 = alpha_opt_weight_2.replace('.','_')    
        beta_opt_weight_2  = str(format(beta_weights[1],'.1f'))
        beta_opt_weight_2  = beta_opt_weight_2.replace('.','_')   
        save_file_name  = '_A_' + alpha_opt_weight_1 + '_B_' + beta_opt_weight_1 + '_A_' + alpha_opt_weight_2 + '_B_' + beta_opt_weight_2
        
        
        # file names 
        fig_1_name = "PR_Twist_Comparison" + save_file_name 
        fig_2_name = "PR_Chord_Comparison"  + save_file_name 
        fig_3_name = "PR_Thickness_Comparison"  + save_file_name   
        fig_4_name = 'PR_Thrust_Comparison' + save_file_name  
        fig_5_name = 'PR_Torque_Comparison'  + save_file_name    
        fig_6_name = 'PR_Total_SPL_Comparison' + save_file_name    
        fig_7_name = "PR_Power_Noise_Pareto"  + save_file_name   
        fig_8_name = 'PR_Power_RPM_Pareto'  + save_file_name  
        
        label_letter = ['A','B']
        for i in range(len(alpha_weights)):   
            alpha_opt_weight = str(format(alpha_weights[i],'.5f'))
            alpha_opt_weight = alpha_opt_weight.replace('.','_')    
            beta_opt_weight  = str(format(beta_weights[i],'.5f'))
            beta_opt_weight  = beta_opt_weight.replace('.','_')    
            rotor_file_name  = 'PR_Alpha_' + alpha_opt_weight + '_Beta_' + beta_opt_weight  
            rotor_label      = 'Rotor ' + label_letter[i]  +', ' + r'$\alpha$ = ' + str(alpha_weights[i])  + r' $\beta$ = ' + str(beta_weights[i])   

            rotor       = load_blade_geometry(folder_name,rotor_file_name)
            plot_rotor_2d_geoemtry_and_performance(rotor,prop_rotor,AXES,PP.colors_1[i],PP.colors_2[i],PP.line_styles[0],
                                                   PP.line_styles[1],PP.markers[0],PP.markers[1],rotor_label)
            
             
        
    else:
        
        AXES , FIGURES  = set_up_axes(PP,'LR') 
        
        # file names 
        fig_1_name = "LR_Twist_Comparison" 
        fig_2_name = "LR_Chord_Comparison"  
        fig_3_name = "LR_Thickness_Comparison"   
        fig_4_name = 'LR_Thrust_Comparison' 
        fig_5_name = 'LR_Torque_Comparison'   
        fig_6_name = 'LR_Total_SPL_Comparison'   
        fig_7_name = "LR_Power_Noise_Pareto"  
        fig_8_name = 'LR_Power_RPM_Pareto' 
        
        try:   
            # Plot Rotor designed using Adkins and Liebeck  
            rotor_file_name  = 'LR_AL'  
            rotor_label      = 'Adkins & Liebeck'
            rotor            = load_blade_geometry(folder_name,rotor_file_name)       
            plot_rotor_2d_geoemtry_and_performance(rotor,prop_rotor,AXES,'black','grey',PP.line_styles[0],PP.line_styles[1],PP.markers[0],PP.markers[1],rotor_label)      
        except: 
            pass        
        
        for i in range(len(alpha_weights)):    
            # save rotor geomtry
            alpha_opt_weight = str(format(alpha_weights[i],'.5f'))
            alpha_opt_weight = alpha_opt_weight.replace('.','_')     
            rotor_file_name  =  'LR_Alpha_' + alpha_opt_weight  
            rotor_label      = r'$\alpha$ = ' + str(alpha_weights[i])  
            try: 
                rotor       = load_blade_geometry(folder_name,rotor_file_name)
                plot_rotor_2d_geoemtry_and_performance(rotor,prop_rotor,AXES,PP.colors_1[i],PP.colors_1[i],PP.line_styles[0],
                                                       PP.line_styles[1],PP.markers[0],PP.markers[1],rotor_label) 
            except: 
                pass   
                   
    AXES[0].legend(loc='upper right')
    AXES[1].legend(loc='upper left')
    AXES[2].legend(loc='upper left')
    AXES[3].legend(loc='upper left',)
    AXES[4].legend(loc='upper left')
    AXES[5].legend(loc='upper right')
    AXES[6].legend(loc='upper right')
    AXES[7].legend(loc='upper right')
        
    
    # axis limits 
    ymin0, ymax0 = AXES[0].get_ylim()
    ymin1, ymax1 = AXES[1].get_ylim()
    ymin2, ymax2 = AXES[2].get_ylim()
    ymin3, ymax3 = AXES[3].get_ylim()
    ymin4, ymax4 = AXES[4].get_ylim()
    ymin5, ymax5 = AXES[5].get_ylim() 
    ymin6, ymax6 = AXES[6].get_ylim()
    ymin7, ymax7 = AXES[7].get_ylim()
    
    xmin6, xmax6 = AXES[6].get_xlim()
    xmin7, xmax7 = AXES[7].get_xlim() 
    
    
    AXES[0].set_ylim([ymin0, ymax0])    
    AXES[1].set_ylim([ymin1, ymax1*1.1])    
    AXES[2].set_ylim([ymin2, ymax2*1.1])    
    AXES[3].set_ylim([ymin3, ymax3*1.5])  
    AXES[4].set_ylim([ymin4, ymax4*1.5])     
    AXES[5].set_ylim([ymin5, ymax5]) 
    AXES[6].set_ylim([ymin6, ymax6])     
    AXES[7].set_ylim([ymin7, ymax7]) 

    AXES[6].set_xlim([xmin6, xmax6])     
    AXES[7].set_xlim([xmin7, xmax7])    
    
    label_location = 0.5

    # axis labels 
    AXES[0].text(0.02,  (ymax0+ymin0)*label_location, 'Root', fontsize=25)    
    AXES[1].text(0.02,  (ymax1+ymin1)*label_location, 'Root', fontsize=25)  
    AXES[2].text(0.02,  (ymax2+ymin2)*label_location, 'Root', fontsize=25)   
    AXES[3].text(0.02,  (ymax3+ymin3)*label_location, 'Root', fontsize=25) 
    AXES[4].text(0.02 , (ymax4+ymin4)*label_location, 'Root', fontsize=25)     
                            
    ## Assign Colormap and save plots 
    #cmap     = plt.get_cmap("viridis")
    #new_cmap = truncate_colormap(cmap, 0.0, 1.0)
    #norm     = plt.Normalize(0,1) 
    #sm       =  ScalarMappable(norm=norm, cmap=new_cmap)
    #ax_ticks = np.linspace(0,1,11)
    #sm.set_array([])  

    #cmap_2      = plt.get_cmap("viridis")
    #new_cmap_2  = truncate_colormap(cmap_2, 0.0, 1.0)
    #norm_2      = plt.Normalize(0.45,0.7) 
    #sm_2        = ScalarMappable(norm=norm_2, cmap=new_cmap_2 )
    #ax_ticks_2  = np.linspace(0.45,0.7,6)
    #sm_2.set_array([])     
            
    #sfmt = ticker.ScalarFormatter(useMathText=True) 
    #sfmt = ticker.FormatStrFormatter('%.1f')     
    #sfmt2 = ticker.ScalarFormatter(useMathText=True) 
    #sfmt2 = ticker.FormatStrFormatter('%.2f')     
    #cbar_1 = FIGURES[0].colorbar(sm,   ax = AXES[0], ticks = list(ax_ticks),  format= sfmt)
    #cbar_2 = FIGURES[1].colorbar(sm,   ax = AXES[1], ticks = list(ax_ticks),  format= sfmt)
    #cbar_3 = FIGURES[2].colorbar(sm,   ax = AXES[2], ticks = list(ax_ticks),  format= sfmt)
    #cbar_4 = FIGURES[3].colorbar(sm,   ax = AXES[3], ticks = list(ax_ticks),  format= sfmt)
    #cbar_5 = FIGURES[4].colorbar(sm,   ax = AXES[4], ticks = list(ax_ticks),  format= sfmt) 
    #cbar_6 = FIGURES[5].colorbar(sm,   ax = AXES[5], ticks = list(ax_ticks),  format= sfmt)
    #cbar_7 = FIGURES[6].colorbar(sm,   ax = AXES[6], ticks = list(ax_ticks),  format= sfmt)
    #cbar_8 = FIGURES[7].colorbar(sm_2, ax = AXES[7], ticks = list(ax_ticks_2),  format= sfmt2)
    
    
    #cbar_1.set_label(r'$\alpha$')
    #cbar_2.set_label(r'$\alpha$')
    #cbar_3.set_label(r'$\alpha$')
    #cbar_4.set_label(r'$\alpha$') 
    #cbar_5.set_label(r'$\alpha$')  
    #cbar_6.set_label(r'$\alpha$')
    #cbar_7.set_label(r'$\alpha$') 
    #cbar_8.set_label(r'Tip Mach')   
     
    FIGURES[0].tight_layout()
    FIGURES[1].tight_layout()
    FIGURES[2].tight_layout()
    FIGURES[3].tight_layout()
    FIGURES[4].tight_layout() 
    FIGURES[5].tight_layout()
    FIGURES[6].tight_layout(rect= (0.05,0,1,1))
    FIGURES[7].tight_layout(rect= (0.05,0,1,1))  
    
    if save_figures:
        FIGURES[0].savefig(fig_1_name  + '.png')               
        FIGURES[1].savefig(fig_2_name  + '.png')               
        FIGURES[2].savefig(fig_3_name  + '.png')               
        FIGURES[3].savefig(fig_4_name  + '.png')          
        FIGURES[4].savefig(fig_5_name  + '.png')                
        FIGURES[5].savefig(fig_6_name  + '.png')               
        FIGURES[6].savefig(fig_7_name  + '.png')            
        FIGURES[7].savefig(fig_8_name  + '.png')    
     
    return   

# ------------------------------------------------------------------ 
# Plot geoemtry and performance of single rotor 
# ------------------------------------------------------------------  
def plot_rotor_2d_geoemtry_and_performance(rotor,prop_rotor,AXES,color,color_2,line_style,line_style_2, marker,marker_2,label = 'prop'):   
    PP       = define_plot_parameters()

    n_root_sections   = 1 # 8
    rotor_modified    = rotor # add_rotor_stem(rotor,number_of_root_sections = n_root_sections)
    c                 = rotor_modified.chord_distribution
    beta              = rotor_modified.twist_distribution/Units.degrees 
    r                 = rotor_modified.radius_distribution/rotor.tip_radius
    t                 = rotor_modified.max_thickness_distribution 
    
    if prop_rotor: 
        prop_rotor_flag    = True 
        T_hover            = to_numpy(rotor_modified.results.hover.full_results.blade_thrust_distribution[0])
        Q_hover            = to_numpy(rotor_modified.results.hover.full_results.blade_torque_distribution[0])  
        SPL_dBA_1_3_hover  = to_numpy(rotor_modified.results.hover.noise_data.SPL_1_3_spectrum_dBA[0,0]) 
        frequency          = to_numpy(rotor_modified.results.hover.noise_data.one_third_frequency_spectrum) 
        RPM                = to_numpy(rotor_modified.results.hover.full_results.omega[0][0])/Units.rpm  
        SPL_max            = to_numpy(rotor_modified.results.hover.mean_SPL) 
        T_cruise           = to_numpy(rotor_modified.results.cruise.full_results.blade_thrust_distribution[0]) 
        Q_cruise           = to_numpy(rotor_modified.results.cruise.full_results.blade_torque_distribution[0])  
        PM_cruise          = to_numpy(rotor_modified.results.cruise.collective) 
        design_power       = to_numpy(rotor_modified.results.cruise.power )
        
        a  = 343
        R  = rotor_modified.tip_radius
        pitch_hover = to_numpy(rotor_modified.results.hover.collective* 180/np.pi ) 
        pitch_oei   = to_numpy(rotor_modified.results.OEI.collective * 180/jnp.pi) 
        pitch_cruise= to_numpy(rotor_modified.results.cruise.collective * 180/jnp.pi) 
        tm_hover    = to_numpy(rotor_modified.results.hover.full_results.omega[0][0])*(R/a)
        tm_oei      = to_numpy(rotor_modified.results.OEI.full_results.omega[0][0])*(R/a)
        tm_cruise   = to_numpy(rotor_modified.results.cruise.full_results.omega[0][0])*(R/a)
        rpm_hover   = to_numpy(rotor_modified.results.hover.full_results.omega[0][0])/Units.rpm 
        rpm_oei     = to_numpy(rotor_modified.results.OEI.full_results.omega[0][0])/Units.rpm 
        rpm_cruise  = to_numpy(rotor_modified.results.cruise.full_results.omega[0][0])/Units.rpm 
        power_hover = to_numpy(rotor_modified.results.hover.power[0][0] )
        power_oei   = to_numpy(rotor_modified.results.OEI.power[0][0] )
        power_cruise= to_numpy(rotor_modified.results.cruise.power[0][0] )  
        hover_spl   = to_numpy(rotor_modified.results.hover.mean_SPL) 
        
                
        print('pitch command hover : ')
        print(pitch_hover  )
        print('pitch command OEI : ')
        print( pitch_oei  )
        print('pitch command cruise: ')
        print(pitch_cruise  ) 
        print('tip Mach hover : ')
        print(tm_hover  )
        print('tip Mach OEI : ')
        print( tm_oei  )
        print('tip Mach cruise: ')
        print(tm_cruise  )    
        print('RPM hover : ')
        print(rpm_hover  )
        print('RPM OEI : ')
        print(rpm_oei  )
        print('RPM cruise: ')
        print(rpm_cruise  )        
        print('power hover : ')
        print( power_hover  )
        print('power OEI : ' )
        print( power_oei  )
        print('power cruise: ')
        print( power_cruise  ) 
        print('hover SPL : ')
        print( hover_spl  )       
        
        

    else:
        prop_rotor_flag   = False 
        T_hover           = to_numpy(rotor_modified.results.hover.full_results.blade_thrust_distribution[0])
        Q_hover           = to_numpy(rotor_modified.results.hover.full_results.blade_torque_distribution[0])  
        SPL_dBA_1_3_hover = to_numpy(rotor_modified.results.hover.noise_data.SPL_1_3_spectrum_dBA[0,0]) 
        frequency         = to_numpy(rotor_modified.results.hover.noise_data.one_third_frequency_spectrum)   
        RPM               = to_numpy(rotor_modified.results.hover.full_results.omega[0][0])/Units.rpm 
        SPL_max           = to_numpy(rotor_modified.results.hover.mean_SPL)
        design_power      = to_numpy(rotor_modified.results.hover.power) 
    
     
    AXES[0].axvspan(0, r[n_root_sections-1], alpha=PP.alpha_val, color= PP.root_color)
    AXES[0].plot(r[:(n_root_sections)], beta[:(n_root_sections)],color = 'grey', linestyle =  line_style,linewidth = PP.line_width)  
    AXES[0].plot(r[(n_root_sections-1):], beta[(n_root_sections-1):],color = color , markersize = PP.marker_size ,marker = marker, linestyle = line_style,linewidth = PP.line_width, label = label )  
    AXES[0].set_ylabel(r'$\beta$ ($\degree$)') 
    AXES[0].set_xlabel('r')    
    AXES[0].set_xlim([0,max(r)])   
     
    AXES[1].axvspan(0, r[n_root_sections-1], alpha=PP.alpha_val, color= PP.root_color)
    AXES[1].plot(r[:(n_root_sections)], c[:(n_root_sections)] ,color = 'grey', linestyle =  line_style,linewidth = PP.line_width)  
    AXES[1].plot(r[(n_root_sections-1):], c[(n_root_sections-1):] ,color = color,marker = marker , markersize = PP.marker_size, linestyle = line_style,linewidth = PP.line_width, label = label )    
    AXES[1].set_ylabel('c (m)') 
    AXES[1].set_xlim([0,max(r)]) 
    AXES[1].set_xlabel('r')    
 
    AXES[2].axvspan(0, r[n_root_sections-1], alpha=PP.alpha_val, color= PP.root_color)
    AXES[2].plot(r[:(n_root_sections)], t[:(n_root_sections)] ,color = 'grey', linestyle = line_style,linewidth = PP.line_width)  
    AXES[2].plot(r[(n_root_sections-1):] , t[(n_root_sections-1):],color = color,marker = marker, markersize = PP.marker_size, linestyle = line_style,linewidth = PP.line_width, label = label )     
    AXES[2].set_ylabel('t (m)')  
    AXES[2].set_xlim([0,max(r)])   
    AXES[2].set_xlabel('r')    
         
     
    AXES[3].axvspan(0, r[n_root_sections-1], alpha=PP.alpha_val, color= PP.root_color)
    if prop_rotor_flag: 
        AXES[3].plot(r[(n_root_sections-1):] , T_hover ,color = color , markersize = PP.marker_size,marker = marker, linestyle =line_style, linewidth = PP.line_width , label = 'Hover: ' + label  )   
        AXES[3].plot(r[(n_root_sections-1):] , T_cruise,color = color_2 , markersize = PP.marker_size,marker = marker_2, linestyle =line_style_2, linewidth = PP.line_width, label = 'Cruise: ' + label  )    
    else: 
        AXES[3].plot(r[(n_root_sections-1):] , T_hover ,color = color, markersize = PP.marker_size,marker = marker, linestyle = line_style, linewidth = PP.line_width , label = label  )    
         
    AXES[3].set_xlim([0,max(r)])   
    
     
    AXES[4].axvspan(0, r[n_root_sections-1], alpha=PP.alpha_val, color= PP.root_color) 
    if prop_rotor_flag: 
        AXES[4].plot(r[(n_root_sections-1):] , Q_hover ,color = color , markersize = PP.marker_size,marker = marker, linestyle =line_style, linewidth = PP.line_width , label =  'Hover: ' + label  )   
        AXES[4].plot(r[(n_root_sections-1):] , Q_cruise,color = color_2 , markersize = PP.marker_size,marker = marker_2, linestyle =line_style_2, linewidth = PP.line_width, label = 'Cruise: ' + label  )    

    else: 
        AXES[4].plot(r[(n_root_sections-1):] , Q_hover ,color = color , markersize = PP.marker_size,marker = marker, linestyle = line_style, linewidth = PP.line_width , label = label  )            
    AXES[4].set_xlim([0,max(r)])   
 
    if prop_rotor_flag: 
        AXES[5].semilogx(frequency , SPL_dBA_1_3_hover ,color = color , markersize = PP.marker_size,marker = marker, linestyle = line_style, linewidth = PP.line_width, label = label   )    
    else: 
        AXES[5].semilogx(frequency , SPL_dBA_1_3_hover,color = color  , markersize = PP.marker_size,marker = marker, linestyle =line_style, linewidth = PP.line_width , label = label  )    
    
    #AXES[6].scatter(design_power/1E3,SPL_max, color  = color, marker = 'o', s      = 150, label  = label)  
       
    #AXES[7].scatter(design_power/1E3, RPM, color  = color,  marker = 'o',  s      = 150,  label  = label)  
    return   


# ------------------------------------------------------------------ 
# Truncate colormaps
# ------------------------------------------------------------------  
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

  
# ------------------------------------------------------------------ 
# Plot rotor 
# ------------------------------------------------------------------ 
def plot_3d_rotor_geometry(rotor,fig_name,save_figure, elevation_angle = 45,  aximuth_angle = 0, cpt=0,rotor_face_color='dodgerblue',
                        rotor_edge_color='dodgerblue',rotor_alpha=1.0):  
 
    n_root_sections   = 8
    rotor             = add_rotor_stem(rotor,number_of_root_sections = n_root_sections)
    
    fig_name = "3D_" + fig_name
    fig = plt.figure(fig_name) 
    fig.set_size_inches(12,8) 
    axes = plt.axes(projection='3d')  
    axes.view_init(elev= elevation_angle, azim= aximuth_angle)    
    RADIUS = 1.75 # Control this value.
    axes.set_xlim3d(-RADIUS / 2, RADIUS / 2)
    axes.set_zlim3d(-RADIUS / 2, RADIUS / 2)
    axes.set_ylim3d(-RADIUS / 2, RADIUS / 2)    
    axes.grid(False)  
    plt.axis('off')
        
    num_B     = rotor.number_of_blades
    n_points  = 21
    af_pts    = n_points-1
    dim       = len(rotor.radius_distribution) 

    for i in range(num_B):
        G = get_blade_coordinates(rotor,n_points,dim,i)
        # ------------------------------------------------------------------------
        # Plot Propeller Blade
        # ------------------------------------------------------------------------
        for sec in range(dim-1):
            for loc in range(af_pts):
                X = [G.XA1[cpt,sec,loc],
                     G.XB1[cpt,sec,loc],
                     G.XB2[cpt,sec,loc],
                     G.XA2[cpt,sec,loc]]
                Y = [G.YA1[cpt,sec,loc],
                     G.YB1[cpt,sec,loc],
                     G.YB2[cpt,sec,loc],
                     G.YA2[cpt,sec,loc]]
                Z = [G.ZA1[cpt,sec,loc],
                     G.ZB1[cpt,sec,loc],
                     G.ZB2[cpt,sec,loc],
                     G.ZA2[cpt,sec,loc]]
                rotor_verts = [list(zip(X, Y, Z))]
                rotor_collection = Poly3DCollection(rotor_verts)
                rotor_collection.set_facecolor(rotor_face_color)
                rotor_collection.set_edgecolor(rotor_edge_color)
                rotor_collection.set_alpha(rotor_alpha)
                axes.add_collection3d(rotor_collection)  
    
    if save_figure:
        fig.savefig(fig_name  + '.png') 
        
    return

# ------------------------------------------------------------------ 
# get blade geometry 
# ------------------------------------------------------------------ 
def get_blade_coordinates(rotor,n_points,dim,i,aircraftRefFrame=True): 
    # unpack
    num_B        = rotor.number_of_blades
    airfoils     = rotor.Airfoils 
    beta         = rotor.twist_distribution + rotor.inputs.pitch_command
    a_o          = rotor.start_angle
    b            = rotor.chord_distribution
    r            = rotor.radius_distribution
    MCA          = rotor.mid_chord_alignment
    t            = rotor.max_thickness_distribution
    a_loc        = rotor.airfoil_polar_stations
    origin       = rotor.origin
    
    if rotor.rotation==1:
        # negative chord and twist to give opposite rotation direction
        b = -b    
        beta = -beta
    
    theta  = np.linspace(0,2*np.pi,num_B+1)[:-1]
    flip_1 =  (np.pi/2)
    flip_2 =  (np.pi/2)

    MCA_2d             = np.repeat(np.atleast_2d(MCA).T,n_points,axis=1)
    b_2d               = np.repeat(np.atleast_2d(b).T  ,n_points,axis=1)
    t_2d               = np.repeat(np.atleast_2d(t).T  ,n_points,axis=1)
    r_2d               = np.repeat(np.atleast_2d(r).T  ,n_points,axis=1)
    airfoil_le_offset  = np.repeat(b[:,None], n_points, axis=1)/2  

    # get airfoil coordinate geometry
    if len(airfoils.keys())>0:
        xpts  = np.zeros((dim,n_points))
        zpts  = np.zeros((dim,n_points))
        max_t = np.zeros(dim)
        for af_idx,airfoil in enumerate(airfoils):
            geometry     = import_airfoil_geometry(airfoil.coordinate_file,n_points)
            locs         = np.where(np.array(a_loc) == af_idx)
            xpts[locs]   = geometry.x_coordinates  
            zpts[locs]   = geometry.y_coordinates  
            max_t[locs]  = geometry.thickness_to_chord 

    else: 
        airfoil_data = compute_naca_4series('2410',n_points)
        xpts         = np.repeat(np.atleast_2d(airfoil_data.x_coordinates) ,dim,axis=0)
        zpts         = np.repeat(np.atleast_2d(airfoil_data.y_coordinates) ,dim,axis=0)
        max_t        = np.repeat(airfoil_data.thickness_to_chord,dim,axis=0)
            
    # store points of airfoil in similar format as Vortex Points (i.e. in vertices)
    max_t2d = np.repeat(np.atleast_2d(max_t).T ,n_points,axis=1)

    xp      = (- MCA_2d + xpts*b_2d - airfoil_le_offset)     # x-coord of airfoil
    yp      = r_2d*np.ones_like(xp)                          # radial location
    zp      = zpts*(t_2d/max_t2d)                            # former airfoil y coord
    
    rotor_vel_to_body = rotor.prop_vel_to_body()
    cpts              = len(rotor_vel_to_body[:,0,0])
    
    matrix        = np.zeros((len(zp),n_points,3)) # radial location, airfoil pts (same y)
    matrix[:,:,0] = xp
    matrix[:,:,1] = yp
    matrix[:,:,2] = zp
    matrix        = np.repeat(matrix[None,:,:,:], cpts, axis=0)

    
    # ROTATION MATRICES FOR INNER SECTION
    # rotation about y axis to create twist and position blade upright
    trans_1        = np.zeros((dim,3,3))
    trans_1[:,0,0] = np.cos(flip_1 - beta)
    trans_1[:,0,2] = -np.sin(flip_1 - beta)
    trans_1[:,1,1] = 1
    trans_1[:,2,0] = np.sin(flip_1 - beta)
    trans_1[:,2,2] = np.cos(flip_1 - beta)
    trans_1        = np.repeat(trans_1[None,:,:,:], cpts, axis=0)

    # rotation about x axis to create azimuth locations
    trans_2 = np.array([[1 , 0 , 0],
                   [0 , np.cos(theta[i] + a_o + flip_2 ), -np.sin(theta[i] +a_o +  flip_2)],
                   [0,np.sin(theta[i] + a_o + flip_2), np.cos(theta[i] + a_o + flip_2)]])
    trans_2 = np.repeat(trans_2[None,:,:], dim, axis=0)
    trans_2 = np.repeat(trans_2[None,:,:,:], cpts, axis=0)

    # rotation about y to orient propeller/rotor to thrust angle (from propeller frame to aircraft frame)
    trans_3 =  rotor_vel_to_body
    trans_3 =  np.repeat(trans_3[:, None,:,: ],dim,axis=1) 
    
    trans     = np.matmul(trans_2,trans_1)
    rot_mat   = np.repeat(trans[:,:, None,:,:],n_points,axis=2)    

    # ---------------------------------------------------------------------------------------------
    # ROTATE POINTS
    if aircraftRefFrame:
        # rotate all points to the thrust angle with trans_3
        mat  =  np.matmul(np.matmul(rot_mat,matrix[...,None]).squeeze(axis=-1), trans_3)
    else:
        # use the rotor frame
        mat  =  np.matmul(rot_mat,matrix[...,None]).squeeze(axis=-1)
    # ---------------------------------------------------------------------------------------------
    # create empty data structure for storing geometry
    G = Data()
    
    # store node points
    G.X  = mat[:,:,:,0] + origin[0][0]
    G.Y  = mat[:,:,:,1] + origin[0][1]
    G.Z  = mat[:,:,:,2] + origin[0][2]
    
    # store points
    G.XA1  = mat[:,:-1,:-1,0] + origin[0][0]
    G.YA1  = mat[:,:-1,:-1,1] + origin[0][1]
    G.ZA1  = mat[:,:-1,:-1,2] + origin[0][2]
    G.XA2  = mat[:,:-1,1:,0]  + origin[0][0]
    G.YA2  = mat[:,:-1,1:,1]  + origin[0][1]
    G.ZA2  = mat[:,:-1,1:,2]  + origin[0][2]

    G.XB1  = mat[:,1:,:-1,0] + origin[0][0]
    G.YB1  = mat[:,1:,:-1,1] + origin[0][1]
    G.ZB1  = mat[:,1:,:-1,2] + origin[0][2]
    G.XB2  = mat[:,1:,1:,0]  + origin[0][0]
    G.YB2  = mat[:,1:,1:,1]  + origin[0][1]
    G.ZB2  = mat[:,1:,1:,2]  + origin[0][2]
    
    return G

# ------------------------------------------------------------------ 
# Add rotor stem 
# ------------------------------------------------------------------  
def add_rotor_stem(rotor,number_of_root_sections= 5):
    
    # define airfoil sections   
    a1            = "Circle_Section.txt"
    a2            = rotor.Airfoils[list(rotor.Airfoils.keys())[0]].coordinate_file    # first airfoil on rotor 
    new_files     = generate_interpolated_airfoils(a1, a2, number_of_root_sections,save_filename="Root_Airfoil")  
    
    
    
    num_sec = 20         
    new_radius_distribution         = np.linspace(rotor.radius_distribution[0],rotor.radius_distribution[-1],num_sec)
    func_twist_distribution         = interp1d(rotor.radius_distribution , rotor.twist_distribution , kind='cubic')
    func_chord_distribution         = interp1d(rotor.radius_distribution , rotor.chord_distribution , kind='cubic') 
    func_mca_distribution           = interp1d(rotor.radius_distribution , rotor.mid_chord_alignment , kind='cubic') 
    func_max_thickness_distribution = interp1d(rotor.radius_distribution , rotor.max_thickness_distribution, kind='cubic')  
    
    rotor.twist_distribution         = func_twist_distribution(new_radius_distribution)     
    rotor.chord_distribution         = func_chord_distribution(new_radius_distribution)         
    rotor.mid_chord_alignment        = func_mca_distribution(new_radius_distribution) 
    rotor.radius_distribution        = new_radius_distribution      
    rotor.max_thickness_distribution = func_max_thickness_distribution(new_radius_distribution) 
    rotor.thickness_to_chord         = rotor.max_thickness_distribution/rotor.chord_distribution 
    
        
    for i in range(number_of_root_sections-1):
        # import geometry  
        airfoil                     = MARC.Components.Airfoils.Airfoil()
        airfoil.coordinate_file     = new_files[i]         
        airfoil.tag                 = 'Root_Section_' + str(i)
        airfoil.geometry            = import_airfoil_geometry(airfoil.coordinate_file )
        # append geometry
        rotor.Airfoils.append(airfoil) 
    
    # modify rotor 
    x      = np.linspace(0,4,number_of_root_sections)  
    func_1 = (np.tanh(x-2) + 2)/3
    func_2 = (np.tanh(x-2) + 1)/3 
    
    root_radius = np.linspace(0.1,rotor.radius_distribution[0],number_of_root_sections)[:-1]
    root_chord  = func_1[:-1]*rotor.chord_distribution[0]
    root_twist  = func_2[:-1]*rotor.twist_distribution[0]
    root_aloc   = list(np.arange(1,number_of_root_sections)) 
 
    # update rotor geoetry  
    rotor.airfoil_polar_stations     = root_aloc +  rotor.airfoil_polar_stations 
    rotor.chord_distribution         = np.hstack(( root_chord, rotor.chord_distribution   )) 
    rotor.twist_distribution         = np.hstack((root_twist , rotor.twist_distribution  ))   
    rotor.radius_distribution        = np.hstack((root_radius , rotor.radius_distribution       ))  
    rotor.mid_chord_alignment        = np.hstack((np.ones(number_of_root_sections - 1)*rotor.mid_chord_alignment[0] , rotor.mid_chord_alignment       ))  
    
    t_max_root  = np.zeros((number_of_root_sections - 1))    
    t_c_root    = np.zeros((number_of_root_sections - 1))    
    if len(rotor.Airfoils.keys())>0:
        for j,airfoil in enumerate(rotor.Airfoils): 
            a_geo              = airfoil.geometry
            if a_geo == None:
                a_geo = import_airfoil_geometry(airfoil.coordinate_file)
            locs               = np.where(np.array(root_aloc) == j )
            t_max_root[locs]   = a_geo.thickness_to_chord*rotor.chord_distribution[locs]
            t_c_root[locs]     = a_geo.thickness_to_chord      
     
    rotor.max_thickness_distribution = np.hstack(( t_max_root, rotor.max_thickness_distribution)) 
    rotor.thickness_to_chord         = np.hstack((t_c_root , rotor.thickness_to_chord        ))  
    
    return rotor   


# ------------------------------------------------------------------ 
# Load data  
# ------------------------------------------------------------------     
def load_blade_geometry(folder_name,filename):  
    ospath    = os.path.abspath(__file__)
    separator = os.path.sep
    rel_path  = os.path.dirname(ospath) + separator     
    load_file = rel_path  +  folder_name + filename + '.pkl'
    with open(load_file, 'rb') as file:
        rotor = pickle.load(file) 
    return rotor

if __name__ == '__main__': 
    main() 
    plt.show()