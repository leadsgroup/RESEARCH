import RCAIDE
from   RCAIDE.Framework.Core                                   import Units, Data
from   RCAIDE.Library.Components.Energy.Sources.Fuel_Cell_Stacks.PEM import PEM_Cell
from   RCAIDE.Library.Components.Energy.Sources.Fuel_Cell_Stacks.CEM import CEM_Module
import  numpy as  np
import  matplotlib.pyplot as  plt
from FC.pressure_dependence import CDs_1, CDs_1_5, CDs_2_5, Vs_1, Vs_1_5, Vs_2_5
from FC.HT_data import V_P1, V_O2, V_P0, i_O2, i_P0, i_P1

def main():
    Arjun_Model()
    
    RCAIDE_model()
    
    
def RCAIDE_model():
    
    
    
    return

def Arjun_Model():
    
    thermo_state          = Data()
        
    thermo_state.Tt       = 298.15
    thermo_state.Tt       = 101325
    
    T_fc                  = 353.15
    P_H2_input_1          = 1
    P_H2_input_2          = 1
    P_H2_input_3          = 1
    P_air_1               = 1
    P_air_2               = 1.5
    p_air_3_5             = 2.5
    RH                    = 1
    CD_end_1              = 1.6
    CD_end_2              = 1.8
    CD_end_3              = 2.175
    lambda_O2             = 2

    t_m                   = 0.0024 #* Units.cm
    t_m2                  = 0.005 #* Units.cm  
    t_m_wang              = 0.0183 #* Units.cm
    alpha                 = 0.375  # Transfer coefficient (dimensionless)
    alpha2                = 0.375
    A                     = 50 #* Units.cm ** 2
    A2                    = 9 #* Units.cm ** 2
    A_wang                = 25# * Units.cm ** 2
    compressor_efficiency = 0.71
    expander_efficiency   = 0.73
    motor_efficiency      = 0.895 * 0.895
    generator_efficiency  = 1
    specific_weight       = None
    
    type_1                = "LT"
    a_c_1                 = 98
    L_c_1                 = 0.1
    lambda_eff_1          = 24
    lambda_eff_2          = 10
    type_2                = "HT"
    a_c_2                 = 315.5
    L_c_2                 = 0.4
    gamma_1               = 0.45
    gamma_2               = 1
    i0ref_1               = 0.09 * 10 ** -4 
    i0ref_2               = 4 * 10 ** -8
    i0ref_P_ref           = 1 #* Units.bar
    i0ref_T_ref_1         = 353 #* Units.K
    i0ref_T_ref_2         = 369 #* Units.K
    
    cem                   = CEM_Module(compressor_efficiency, expander_efficiency, motor_efficiency, generator_efficiency, specific_weight)
    fuel_cell             = PEM_Cell(type_1, t_m,      a_c_1, L_c_1, A,      cem, lambda_eff = lambda_eff_1, gamma = gamma_1, i0ref = i0ref_1, i0ref_T_ref = i0ref_T_ref_1, alpha=alpha)
    fuel_cell_2           = PEM_Cell(type_2, t_m2,     a_c_2, L_c_2, A2,     cem,               gamma = gamma_2, i0ref = i0ref_2, i0ref_P_ref = i0ref_P_ref, i0ref_T_ref = i0ref_T_ref_2, alpha=alpha2)
    fuel_cell_wang        = PEM_Cell(type_1, t_m_wang, a_c_1, L_c_1, A_wang, cem, lambda_eff = lambda_eff_2, gamma = gamma_1, i0ref = i0ref_1,              i0ref_T_ref = i0ref_T_ref_1)    
           
    compute_fuel_cell_performance(fuel_cell, T_fc, P_H2_input_1,P_H2_input_2,P_H2_input_3, P_air_1,P_air_2,p_air_3_5, RH, CD_end_1, CD_end_2,CD_end_3, thermo_state, lambda_O2)
    
    return 
    
def compute_fuel_cell_performance(fuel_cell, T_fc, P_H2_input_1,P_H2_input_2,P_H2_input_3, P_air_1,P_air_2,p_air_3_5, RH, CD_end_1, CD_end_2,CD_end_3, thermo_state, lambda_O2):
    rated_CD_1, rated_PD1 = fuel_cell.evaluate_max_gross_power(T_fc, P_H2_input_1, P_air_1, RH, lambda_O2)
    rated_CD_2, rated_PD2 = fuel_cell.evaluate_max_gross_power(T_fc, P_H2_input_2, P_air_2, RH, lambda_O2)
    rated_CD_3, rated_PD3 = fuel_cell.evaluate_max_gross_power(T_fc, P_H2_input_3, p_air_3_5, RH, lambda_O2)
    
    fuel_cell.set_rated_cd(rated_CD_1, rated_PD1)
    CDs_1_model, Vs_1_model = fuel_cell.generate_gross_polarization_curve_data(T_fc, P_H2_input_1, P_air_1, RH, lambda_O2, CD_end_1, thermo_state)
    fuel_cell.set_rated_cd(rated_CD_2, rated_PD2)
    CDs_2_model, Vs_2_model = fuel_cell.generate_gross_polarization_curve_data(T_fc, P_H2_input_2, P_air_2, RH, lambda_O2, CD_end_2, thermo_state)
    fuel_cell.set_rated_cd(rated_CD_3, rated_PD2)
    CDs_3_5_model, Vs_3_5_model = fuel_cell.generate_gross_polarization_curve_data(T_fc, P_H2_input_3, p_air_3_5, RH, lambda_O2, CD_end_3, thermo_state)

    powers_1_model = CDs_1_model * Vs_1_model
    powers_2_model = CDs_2_model * Vs_2_model
    powers_3_5_model = CDs_3_5_model * Vs_3_5_model
    Ps_1 = np.array(CDs_1) * np.array(Vs_1)
    Ps_1_5 = np.array(CDs_1_5) * np.array(Vs_1_5)
    Ps_2_5 = np.array(CDs_2_5) * np.array(Vs_2_5)
    
    fig, ax = plt.subplots() 
    ax.set_xlabel(f"Current Density (A/cm$^2$)")
    ax.set_ylabel(f"Voltage (V)")
    color1 = 'black'
    color2 = 'gray'
    color3 = 'orange'
    ax.scatter(CDs_1, Vs_1, label = 'Experimental', color = color1)
    ax.plot(CDs_1_model, Vs_1_model, label = 'Model', color = color1, linestyle = 'solid')

    ax.scatter(CDs_1_5, Vs_1_5, color = color2)
    ax.plot(CDs_2_model, Vs_2_model, color = color2, linestyle = 'solid')

    ax.scatter(CDs_2_5, Vs_2_5, color = color3)
    ax.plot(CDs_3_5_model, Vs_3_5_model, color = color3, linestyle = 'solid')
    ax.set_ylim(0.4, 0.95)
    #ax.set_xlim(-0.05, 1.75)
    ax.annotate("2.5", xy=(2.125, 0.425), color = color3)
    ax.annotate("1.5", xy=(1.75, 0.425), color = color2)
    ax.annotate("$p_{stack}$ (bar):", xy=(1, 0.425), color=color1)
    ax.annotate("1.0", xy=(1.55, 0.425), color = color1)
    ax.legend()
    #ax.set_xlim(0, 0.3)
    
    #plt.show()
    
    fig, ax = plt.subplots() 
    ax.set_xlabel("Current Density (A/cm$^2$)")
    ax.set_ylabel(f"Power Density (W/cm$^2$)")

    ax.scatter(CDs_1, Ps_1, label = 'Experimental', color=color1)
    ax.plot(CDs_1_model, powers_1_model, label = 'Model', color = color1, linestyle = 'solid')

    ax.scatter(CDs_1_5, Ps_1_5, color=color2)
    ax.plot(CDs_2_model, powers_2_model, color = color2, linestyle = 'solid')
    
    ax.scatter(CDs_2_5, Ps_2_5, color=color3)
    ax.plot(CDs_3_5_model, powers_3_5_model, color = color3, linestyle = 'solid')
    
    ax.set_xlim(-0.05, 2.3)
    #ax.set_ylim(-0.05, 0.8)
    ax.annotate("2.5 bar", xy=(1.95, 0.94), color = color3)
    ax.annotate("1.5 bar", xy=(1.85, 0.8), color = color2)
    #ax.annotate("$p_{stack}$ (kPa):", xy=(, 0.37), color=color1)
    ax.annotate("1.0 bar", xy=(1.65, 0.7), color = color1)
    #ax.legend(loc=2)
    
    plt.show()
    


    
    return

if __name__ == '__main__': 
    main()