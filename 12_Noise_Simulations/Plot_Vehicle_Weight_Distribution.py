# MARC Imports 
from  MARC.Core import Data, Units  

# Python Imports  
import numpy as np   
import pickle 
import matplotlib.pyplot as plt  
import matplotlib.cm as cm  

# import propeller/rotors geometries  
import sys 
 
# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main():  

    save_figure = True
    
    # Universal Plot Settings 
    plt.rcParams['axes.linewidth'] = 2.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 32,
                  'xtick.labelsize': 28,
                  'ytick.labelsize': 28,
                  'axes.titlesize': 32}
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
     
     
     
    aircraft   = ['Stopped_Rotor_CRM','Hexacopter_CRM','Tiltrotor_CRM','CTOL_CRM']   
    
    for i in range(len(aircraft)):  
        vehicle = load_aircraft_geometry(aircraft[i])  
    
        plot_weight_buildup(vehicle,
                            save_filename = "Weight_Breakdown_"+ aircraft[i],
                            SI_Units   = True,
                            save_figure = save_figure)
    
    
    return  

def plot_weight_buildup(vehicle,
                        save_figure = False,
                        show_legend=True,
                        SI_Units   = True,
                        save_filename = "Weight_Breakdown",
                        aircraft_name  = None,
                        file_type = ".png",
                        width = 7.5, height = 7.2): 
  
    
    b         =  vehicle.weight_breakdown    
    weight    =  vehicle.mass_properties.max_takeoff    
    if aircraft_name == None:
        aircraft_name = vehicle.tag
        
    vals_unorm =  np.array([b.rotors,
                            b.hubs,
                            b.booms,
                            b.fuselage,
                            b.landing_gear,
                            b.wings_total,
                            b.seats,
                            b.avionics,
                            b.ECS,
                            b.motors,
                            b.servos,
                            b.wiring ,
                            b.BRS,
                            b.battery,
                            b.payload,
                            b.passengers])  
    

    fig= plt.figure(save_filename)
    fig.set_size_inches(width,height)
    ax = fig.add_subplot(111, projection='polar') 
    size  = 0.2
    vals  = vals_unorm/np.sum(vals_unorm)*2*np.pi
    
    inner_colors  = cm.gray(np.linspace(0.2,0.8,3))  
    middle_colors = cm.viridis(np.linspace(0.2,0.8,9))  
    outer_colors  = cm.plasma(np.linspace(0.2,0.8,6))  
    
    inner_vals   = np.array([np.sum(vals[0:14]),vals[14],vals[15]])
    middle_vals  = np.array([np.sum(vals[0:6]),vals[6],vals[7],vals[8],vals[9],vals[10],vals[11],vals[12],vals[13] ])
    outer_vals   = vals[0:6]   
    
    inner_v  = np.cumsum(np.append(0, inner_vals[:-1]))
    middle_v = np.cumsum(np.append(0, middle_vals[:-1])) 
    outer_v = np.cumsum(np.append(0, outer_vals[:-1]))     
    w1 = ax.bar(x=inner_v ,
           width=inner_vals, bottom=1-3*size, height=size,
           color=inner_colors, edgecolor='w', linewidth=1, align="edge",
           label = ['OEW = ' + str(round(100*np.sum(vals_unorm[0:14])/weight,1)) + ' %',
                    'Payload = ' + str(round(100*vals_unorm[14]/weight,1)) + ' %',
                    'Passengers = ' + str(round(100*vals_unorm[15]/weight,1)) + ' %',] )   
    
    w2 =  ax.bar(x=middle_v,
           width=middle_vals, bottom=1-2*size, height=size,
           color=middle_colors, edgecolor='w', linewidth=1, align="edge",
           label = ['Structural = ' + str(round(100*np.sum(vals_unorm[0:6])/weight,1)) + ' %',
                    'Seats = '    + str(round(100*vals_unorm[6]/weight ,1)) + ' %',
                    'Avionics = ' + str(round(100*vals_unorm[7]/weight ,1)) + ' %',
                    'E.C.S. = '   + str(round(100*vals_unorm[8]/weight ,1)) + ' %',
                    'Motors = '   + str(round(100*vals_unorm[9]/weight ,1)) + ' %',
                    'Servos = '   + str(round(100*vals_unorm[10]/weight,1)) + ' %',
                    'Wiring = '   + str(round(100*vals_unorm[11]/weight,1)) + ' %',
                    'B.R.S. = '   + str(round(100*vals_unorm[12]/weight,1)) + ' %',
                    'Battery = '  + str(round(100*vals_unorm[13]/weight,1)) + ' %',]) 
    
    w3 = ax.bar(x=outer_v,
           width=outer_vals, bottom=1-size, height=size,
           color=outer_colors, edgecolor='w', linewidth=1, align="edge",
           label = ['Rotors = '    + str(round(100*vals_unorm[0]/weight,1)) + ' %',
                    'Hubs = '    + str(round(100*vals_unorm[1]/weight,1)) + ' %',
                    'Booms = '    + str(round(100*vals_unorm[2]/weight,1)) + ' %',
                    'Fuselage = '    + str(round(100*vals_unorm[3]/weight,1)) + ' %',
                    'Landing Gear = '    + str(round(100*vals_unorm[4]/weight,1)) + ' %',
                    'Wings = '    + str(round(100*vals_unorm[5]/weight,1)) + ' %'] )
     
    # Add weight of aircraft 
    weight_text = str((round(weight,2))) + ' kg'
    if not SI_Units:    
        weight = weight/Units.lbs
        weight_text = str((round(weight,1))) + ' lbs.' 
    ax.annotate('MTOW', (np.pi*3/4,0.2), size= 20) 
    ax.annotate(weight_text, (np.pi,0.2), size= 20)   
    
    if show_legend:
        lns = w1+w2+w3
        labs = [l.get_label() for l in lns] 
        #ax.legend(lns, labs, loc='center right', ncol = 1, prop={'size': 14}  ,bbox_to_anchor= (0.5, 0.35, 1, 0.3)  ) 
        ax.legend(lns, labs, loc='lower center', ncol = 3, prop={'size': 14}  ,bbox_to_anchor= (0,-.25, 1, 0.3)  )
    ax.set_axis_off()  
    
    fig.tight_layout()
    
    if save_figure:
        plt.savefig(save_filename + '.pdf')   
        
    
    return 

def load_aircraft_geometry(filename):  
    load_file = filename + '.pkl' 
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results

if __name__ == '__main__': 
    main()    
    plt.show()   
 
