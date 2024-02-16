import numpy as np
from RCAIDE.Core import Units, Data 
import matplotlib.pyplot as plt 
from scipy.interpolate import griddata  
import matplotlib.colors
import matplotlib.colors as colors  
from geopy.distance import geodesic as GD
import matplotlib.ticker as mticker
import matplotlib.cm as cm 
import pickle


# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():
    
    # Universal Plot Settings
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"]    = "Times New Roman"
    parameters = {'axes.labelsize': 26,
                  'xtick.labelsize': 24,
                  'ytick.labelsize': 24,
                  'axes.titlesize': 26}
    plt.rcParams.update(parameters)
    plot_parameters                  = Data()
    plot_parameters.line_width       = 2
    plot_parameters.line_styles      = ['-','--',':'] 
    plot_parameters.line_colors      = cm.inferno(np.linspace(0.2,0.8,3))     
    plot_parameters.markers          = ['o','P','s','^','p','^','D','X','*']
    plot_parameters.figure_width     = 7
    plot_parameters.figure_height    = 6
    plot_parameters.marker_size      = 10 
    plot_parameters.legend_font_size = 20 
    
    # topography data  
    use_lat_long                = True
    save_figure                 = True 
    file_type                   = '.png'
     
    # tag for city  
    city                        = 'LA' 
    aircraft_models             = ['SR','TR','HC'] 
    altitudes                   = ['1000']
    flight_intervals            = ['10']  
    
    eval_coordinates = np.array([[241.942, 33.8493],[242.0820114,33.81151],[242.07653,34.10682]])
    eval_tags        = ['A','B','C']
    
    plot_time(city,aircraft_models,altitudes,
                                  flight_intervals, use_lat_long,save_figure,file_type,plot_parameters,eval_coordinates,eval_tags)    
                          
    return 

 

def plot_time(city,aircraft_models,altitudes,flight_intervals, use_lat_long,save_figure,file_type,PP,eval_coordinates,eval_tags):  

    

    time_step         = 2 # seconds   
    T                 = 2*Units.hours
    num_time_steps    = int(T/time_step)
    time              = (np.linspace(0,2*Units.hours,num_time_steps))/(60*60*24)   
    
    for c_i in range(len(eval_coordinates)):
        fig_name = 'time_varying_noise_' + eval_tags[c_i]  
        fig = plt.figure(fig_name)
        fig.set_size_inches(PP.figure_width,PP.figure_height)   
        axis = fig.add_subplot(1,1,1)  
        
        for m_i in range(len(aircraft_models)): 
            for f_i in range(len(flight_intervals)): 
                
                aircraft  = aircraft_models[m_i] 
                altitude  = '1000'
                frequency = flight_intervals[f_i]
            
                # load data  
                processed_results_filename = 'Raw_Data_' + aircraft + '/'+ aircraft + '_' + altitude + 'ft_' + city + '_' + frequency +   '2_hrs' 
                Results = load_results(processed_results_filename)
                
                x_idx = abs(Results.lat_deg[:,0]   - eval_coordinates[c_i,1]).argmin()
                y_idx = abs(Results.long_deg[0,:]   - eval_coordinates[c_i,0]).argmin()
      
                label_tag  =  aircraft  #+ ' '  + flight_intervals[f_i] + ' min interval flight'
                axis.plot(time[0:3600],Results.cumulative_noise_exposure[:,x_idx,y_idx][0:3600] , color = PP.line_colors[m_i], linestyle = PP.line_styles[f_i],label =label_tag , linewidth = PP.line_width )  
                axis.xaxis.set_major_formatter(mticker.FuncFormatter(timeformat))
        axis.set_xlim([0,1/24])   
        axis.set_ylim([40,80]) 
        axis.set_ylabel(r'L$_{A}$ [dBA]')
        axis.set_xlabel('Time [hh:mm]')            
        axis.legend(loc='upper right', prop={'size': PP.legend_font_size}, ncol= 3)  
        fig.tight_layout()  
        if save_figure: 
            fig_name = 'time_varying_noise_' + eval_tags[c_i]            
            figure_title  = '../Papers_and_Presentation/Images/'  + fig_name 
                    
            plt.savefig(figure_title + file_type)    
 
  
    return  


def timeformat(x,pos=None):
    h = int(x*24.)
    m = int((x*24.-h)*60)
    return "{:02d}:{:02d}".format(h,m)


# ------------------------------------------------------------------
# SUPPORTING FUNCTIONS 
# ------------------------------------------------------------------ 

def get_topography_data(topography_file,N_lat,N_long ): 
 
    data = np.loadtxt(topography_file)
    Long = data[:,0]
    Lat  = data[:,1]
    Elev = data[:,2]  
    
    x_min_coord = np.min(Lat)
    x_max_coord = np.max(Lat)
    y_min_coord = np.min(Long)
    y_max_coord = np.max(Long)
    if np.min(Long)>180: 
        y_min_coord = np.min(Long)-360
    if np.max(Long)>180:
        y_max_coord = np.max(Long)-360  
    
    top_left_map_coords      = np.array([x_max_coord,y_min_coord])
    bottom_left_map_coords   = np.array([x_min_coord,y_min_coord])  
    bottom_right_map_coords  = np.array([x_min_coord,y_max_coord]) 
    
    x_dist_max = GD(top_left_map_coords,bottom_left_map_coords).m 
    y_dist_max = GD(bottom_right_map_coords,bottom_left_map_coords).m  
     
    [long_dist,lat_dist]  = np.meshgrid(np.linspace(0,y_dist_max,N_long),np.linspace(0,x_dist_max,N_lat))
    [long_deg,lat_deg]    = np.meshgrid(np.linspace(np.min(Long),np.max(Long),N_long),np.linspace(np.min(Lat),np.max(Lat),N_lat)) 
    elevation             = griddata((Lat,Long), Elev, (lat_deg, long_deg), method='linear')     
     
    return long_dist,lat_dist,long_deg,lat_deg,elevation

 
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


# ------------------------------------------------------------------ 
# Truncate colormaps
# ------------------------------------------------------------------  
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

class FixPointNormalize(matplotlib.colors.Normalize):
    """ 
    Inspired by https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
    Subclassing Normalize to obtain a colormap with a fixpoint 
    somewhere in the middle of the colormap.
    This may be useful for a `terrain` map, to set the "sea level" 
    to a color in the blue/turquise range. 
    """
    def __init__(self, vmin=None, vmax=None, sealevel=0, col_val = 0.21875, clip=False):
        # sealevel is the fix point of the colormap (in data units)
        self.sealevel = sealevel
        # col_val is the color value in the range [0,1] that should represent the sealevel.
        self.col_val = col_val
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.sealevel, self.vmax], [0, self.col_val, 1]
        return np.ma.masked_array(np.interp(value, x, y)) 
    
if __name__ == '__main__': 
    main()   
    plt.show()
