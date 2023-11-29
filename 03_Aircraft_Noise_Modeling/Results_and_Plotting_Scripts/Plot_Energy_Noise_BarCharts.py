import numpy as np
import matplotlib.pyplot as plt
from MARC.Core import Units , Data 
import matplotlib.cm as cm 
import pickle

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():
    

    save_figure                 = True 
    file_type                   = '.png'

    
    # Universal Plot Settings
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"]    = "Times New Roman"
    parameters = {'axes.labelsize': 26,
                  'xtick.labelsize': 24,
                  'ytick.labelsize': 24,
                  'axes.titlesize': 26}
    plt.rcParams.update(parameters)
    plot_parameters                  = Data()
    plot_parameters.line_width       = 3
    plot_parameters.line_styles      = ['-','--',':'] 
    plot_parameters.line_colors      = cm.inferno(np.linspace(0.2,0.8,3))    
    plot_parameters.line_colors_1    = cm.inferno(np.linspace(0.1,0.9,9))          
    plot_parameters.markers          = ['o','P','s','^','p','^','D','X','*']
    plot_parameters.figure_width     = 11 
    plot_parameters.figure_height    = 6
    plot_parameters.marker_size      = 10 
    plot_parameters.legend_font_size = 20     
    
    
    plot_percent_noise_increase_charts(plot_parameters,save_figure,file_type) 
    
    plot_total_energy_consumtion_charts(plot_parameters,save_figure,file_type)
    plot_energy_consumtion_charts(plot_parameters,save_figure,file_type)
    
    return 

def plot_percent_noise_increase_charts(PP,save_figure,file_type) : 

    labels = ['45-50','50-55','55-60','60-70','> 70'] 
    

    
    # load data  
    SR_f1   = 'Raw_Data_SR/SR_1000ft_LA_10min_All' 
    SR_f2   = 'Raw_Data_SR/SR_1000ft_LA_30min_All' 
    SR_f3   = 'Raw_Data_SR/SR_1000ft_LA_60min_All'  
    SR_10min = load_results(SR_f1)
    SR_30min = load_results(SR_f2)
    SR_60min = load_results(SR_f3)
    

    TR_f1   = 'Raw_Data_TR/TR_1000ft_LA_10min_All' 
    TR_f2   = 'Raw_Data_TR/TR_1000ft_LA_30min_All' 
    TR_f3   = 'Raw_Data_TR/TR_1000ft_LA_60min_All'  
    TR_10min = load_results(TR_f1)
    TR_30min = load_results(TR_f2)
    TR_60min = load_results(TR_f3)
    

    HC_f1   = 'Raw_Data_HC/HC_1000ft_LA_10min_All' 
    HC_f2   = 'Raw_Data_HC/HC_1000ft_LA_30min_All' 
    HC_f3   = 'Raw_Data_HC/HC_1000ft_LA_60min_All'  
    HC_10min = load_results(HC_f1)
    HC_30min = load_results(HC_f2)
    HC_60min = load_results(HC_f3)    
    
    
    L_AeqT_old =  np.array([[[0.3671,  0.0936,  1.0341, 0.0,  0.0],[0.2273,  0.1872,  0.5170,  0.0,  0.0],[57.9646,  20.2247,  2.68872,  0.19342 ,  0.0 ]],
                        [[ 0.7343, 1.0767, 1.1375, 0.1934,  0.0 ],[ 0.5070,  0.7022,  0.8273,  0.0,  0.0 ],[76.551,  110.67,  9.9276, 1.7408, 0.5025 ]],
                        [[1.3638,  1.4044,  1.3443,  2.1276, 0.0 ],[0.577,  0.468, 0.723,  0.386, 0.0 ], [68.54345,  251.5917,  148.9141, 11.02514,  0.502512]]])
    
    L_AeqT      = np.zeros((3,3,5))
    L_AeqT[0,0] = np.array([SR_60min.increase_band_45_50,-SR_60min.increase_band_50_55,SR_60min.increase_band_55_60,SR_60min.increase_band_60_70,SR_60min.increase_band_70_80])
    L_AeqT[0,1] = np.array([TR_60min.increase_band_45_50,TR_60min.increase_band_50_55,TR_60min.increase_band_55_60,TR_60min.increase_band_60_70,TR_60min.increase_band_70_80]) 
    L_AeqT[0,2] = np.array([HC_60min.increase_band_45_50,HC_60min.increase_band_50_55,HC_60min.increase_band_55_60,HC_60min.increase_band_60_70,HC_60min.increase_band_70_80]) 
    L_AeqT[1,0] = np.array([SR_30min.increase_band_45_50,SR_30min.increase_band_50_55,SR_30min.increase_band_55_60,SR_30min.increase_band_60_70,SR_30min.increase_band_70_80])
    L_AeqT[1,1] = np.array([TR_30min.increase_band_45_50,TR_30min.increase_band_50_55,TR_30min.increase_band_55_60,TR_30min.increase_band_60_70,TR_30min.increase_band_70_80]) 
    L_AeqT[1,2] = np.array([HC_30min.increase_band_45_50,HC_30min.increase_band_50_55,HC_30min.increase_band_55_60,HC_30min.increase_band_60_70,HC_30min.increase_band_70_80]) 
    L_AeqT[2,0] = np.array([SR_10min.increase_band_45_50,SR_10min.increase_band_50_55,SR_10min.increase_band_55_60,SR_10min.increase_band_60_70,SR_10min.increase_band_70_80])
    L_AeqT[2,1] = np.array([TR_10min.increase_band_45_50,TR_10min.increase_band_50_55,TR_10min.increase_band_55_60,TR_10min.increase_band_60_70,TR_10min.increase_band_70_80]) 
    L_AeqT[2,2] = np.array([HC_10min.increase_band_45_50,HC_10min.increase_band_50_55,HC_10min.increase_band_55_60,HC_10min.increase_band_60_70,HC_10min.increase_band_70_80])    
    bar_width = 0.3
    opacity = 0.8
    
    aircraft = ['SR', 'TR', 'HC'] 
    for i in range(len(aircraft)): 
        X = np.arange(5)
        fig_name = aircraft[i] + '_noise_increase'
        fig = plt.figure(fig_name)
        fig.set_size_inches(7,6)    
        ax = plt.axes()
        plt.bar(X - 0.3, L_AeqT[0,i], bar_width , color = PP.line_colors[0], label = '60 min')
        plt.bar(X      , L_AeqT[1,i], bar_width, color = PP.line_colors[1] , label =  '30 min')
        plt.bar(X + 0.3, L_AeqT[2,i], bar_width, color = PP.line_colors[2] ,  label = '10 min')
        plt.ylabel('% L$_{AeqT}$ Increase')
        plt.xlabel('L$_{AeqT}$ Range')
        plt.xticks(X, ('45-50','50-55','55-60','60-70','> 70'))
        plt.legend( prop={'size': PP.legend_font_size})  
        plt.tight_layout()
        if save_figure:  
            figure_title  = '../Papers_and_Presentation/Images/'  + fig_name  
            plt.savefig(figure_title + file_type)    
    
    return 


def plot_total_energy_consumtion_charts(PP,save_figure,file_type) : 
   
    bar_width = 0.3 
    
    # load data  
    SR_f1   = 'Raw_Data_SR/SR_1000ft_LA_10min_All' 
    SR_f2   = 'Raw_Data_SR/SR_1000ft_LA_30min_All' 
    SR_f3   = 'Raw_Data_SR/SR_1000ft_LA_60min_All'  
    SR_10min = load_results(SR_f1)
    SR_30min = load_results(SR_f2)
    SR_60min = load_results(SR_f3)
    

    TR_f1   = 'Raw_Data_TR/TR_1000ft_LA_10min_All' 
    TR_f2   = 'Raw_Data_TR/TR_1000ft_LA_30min_All' 
    TR_f3   = 'Raw_Data_TR/TR_1000ft_LA_60min_All'  
    TR_10min = load_results(TR_f1)
    TR_30min = load_results(TR_f2)
    TR_60min = load_results(TR_f3)
    

    HC_f1   = 'Raw_Data_HC/HC_1000ft_LA_10min_All' 
    HC_f2   = 'Raw_Data_HC/HC_1000ft_LA_30min_All' 
    HC_f3   = 'Raw_Data_HC/HC_1000ft_LA_60min_All'  
    HC_10min = load_results(HC_f1)
    HC_30min = load_results(HC_f2)
    HC_60min = load_results(HC_f3) 
    
    fig = plt.figure('Total_Energy')
    fig.set_size_inches(12,6)    
    X = np.arange(3)
    ax = plt.axes()     
    old_SR_energy = np.array([0,0,0])
    old_TR_energy = np.array([0,0,0])
    old_HC_energy = np.array([0,0,0])
            
    SR_energy = np.array([SR_10min.total_energy_consumed, SR_30min.total_energy_consumed,SR_60min.total_energy_consumed])*2.77778E-10
    TR_energy = np.array([TR_10min.total_energy_consumed, TR_30min.total_energy_consumed,TR_60min.total_energy_consumed])*2.77778E-10
    HC_energy = np.array([HC_10min.total_energy_consumed, HC_30min.total_energy_consumed,HC_60min.total_energy_consumed])*2.77778E-10
        
        
    plt.bar(X - 0.3, SR_energy, bar_width ,bottom= old_SR_energy , color = PP.line_colors[0] , label = 'Stopped-rotor')
    plt.bar(X      , TR_energy, bar_width ,bottom= old_TR_energy , color = PP.line_colors[1] , label = 'Tilt-rotor')
    plt.bar(X + 0.3, HC_energy, bar_width ,bottom= old_HC_energy , color = PP.line_colors[2] , label = 'Hexacopter')
    
    plt.xticks(X, ('10 min','30 min','60 min'))   
            
    ax.set_ylabel('Total Energy [MW-hr]')
    ax.set_xlabel('Flight Frequency ')   
    plt.legend( prop={'size': PP.legend_font_size})  
    plt.tight_layout()
    if save_figure:  
        figure_title  = '../Papers_and_Presentation/Images/total_energy'  
        plt.savefig(figure_title + file_type)    
    
    
def plot_energy_consumtion_charts(PP,save_figure,file_type) : 
   
    bar_width = 0.3
    opacity   = 0.8
    
    # load data  
    SR_f1   = 'Raw_Data_SR/SR_1000ft_LA_10min_All' 
    SR_f2   = 'Raw_Data_SR/SR_1000ft_LA_30min_All' 
    SR_f3   = 'Raw_Data_SR/SR_1000ft_LA_60min_All'  
    SR_10min = load_results(SR_f1)
    SR_30min = load_results(SR_f2)
    SR_60min = load_results(SR_f3)
    

    TR_f1   = 'Raw_Data_TR/TR_1000ft_LA_10min_All' 
    TR_f2   = 'Raw_Data_TR/TR_1000ft_LA_30min_All' 
    TR_f3   = 'Raw_Data_TR/TR_1000ft_LA_60min_All'  
    TR_10min = load_results(TR_f1)
    TR_30min = load_results(TR_f2)
    TR_60min = load_results(TR_f3)
    

    HC_f1   = 'Raw_Data_HC/HC_1000ft_LA_10min_All' 
    HC_f2   = 'Raw_Data_HC/HC_1000ft_LA_30min_All' 
    HC_f3   = 'Raw_Data_HC/HC_1000ft_LA_60min_All'  
    HC_10min = load_results(HC_f1)
    HC_30min = load_results(HC_f2)
    HC_60min = load_results(HC_f3)     
    
    fig = plt.figure('Energy_per_Route')
    fig.set_size_inches(12,6)    
    X = np.arange(9)
    ax = plt.axes()      
        
    plt.plot(X, SR_60min.energy_per_flight_per_route[0:9]*2.77778E-7,linewidth = PP.line_width ,color = PP.line_colors[0] , label = 'Stopped-rotor')
    plt.plot(X, TR_60min.energy_per_flight_per_route[0:9]*2.77778E-7,linewidth = PP.line_width ,color = PP.line_colors[1] , label = 'Tilt-rotor')
    plt.plot(X, HC_60min.energy_per_flight_per_route[0:9]*2.77778E-7,linewidth = PP.line_width ,color = PP.line_colors[2] , label = 'Hexacopter') 
        
    
    plt.xticks(X, ( '1',
                    '2',
                    '3',
                    '4',
                    '5',
                    '6',
                    '7',
                    '8',
                    '9'))   
            
    ax.set_ylabel('Energy per Route [kW-hr]')
    ax.set_xlabel('Route No.')   
    
    plt.legend( prop={'size': PP.legend_font_size})  
    plt.tight_layout()
    if save_figure:  
        figure_title  = '../Papers_and_Presentation/Images/energy_per_route'  
        plt.savefig(figure_title + file_type)        
    return  


# ------------------------------------------------------------------
#   Load Results
# ------------------------------------------------------------------   
def load_results(filename):  
    load_file = filename + '.pkl' 
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results  


    
if __name__ == '__main__': 
    main()   
    plt.show()
    
