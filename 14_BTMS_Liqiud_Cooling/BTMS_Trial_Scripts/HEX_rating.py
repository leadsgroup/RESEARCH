import numpy as np

def main(): 
   
   # take in variables from the output of sizing problem
   
   H     = 1.0
   L_c   = 0.3
   L_h   = 0.3
   
   # Inital assumed efficiency of HEX 
   eff_hex     = 0.75
   
   # Input Temperarutes (K)
   T_i_h       = 1173 #hex_opt.coolant_temperature_of_hot_fluid
   T_i_c       = 473 #hex_opt.inlet_temperature_of_cold_fluid
   
   # Inlet Pressures (Pa)
   P_i_h       = 160e3 #hex_opt.coolant_inlet_pressure 
   P_i_c       = 200e3 #hex_opt.air_inlet_pressure
   
   # Universal Gas Constant 
   R           = 287.04 #J/kg.K 
   
   # Hydraulic Diameters
   d_h_c       = 1.54e-3  #m #hex_opt.coolant_hydraulic_diameter 
   d_h_h       = 1.54e-3  #m #hex_opt.air_hydraulic_diameter 

   # Fin Height/Spaceing 
   b_c         = 2.54e-3 #m #hex_opt.fin_spacing_cold                                       
   b_h         = 2.54e-3 #m #hex_opt.fin_spacing_hot

   # Fin metal thickness
   delta_h     = 1e-4   # m #hex_opt.fin_metal_thickness_hot  
   delta_c     = 1e-4   # m #hex_opt.fin_metal_thickness_cold 

   # Platethickness
   delta_w     = 5e-4   # m #hex_opt.t_w 

   # Strip edge exposed 
   l_s_h       = 3.175e-3 #m #hex_opt.fin_exposed_strip_edge_hot  
   l_s_c       = 3.175e-3 #m #hex_opt.fin_exposed_strip_edge_cold 

   #Fin and wall Conductivity 
   k_f         = 18    # W/m.K #hex_opt.k_f  
   k_w         = 18    # W/m.K #hex_opt.k_w 

   # Ratio of finned area to total area 
   Af_A_h      = 0.785 #hex_opt.finned_area_to_total_area_hot  
   Af_A_c      = 0.785 #hex_opt.finned_area_to_total_area_cold        

   #Finned area density 
   beta_h      = 2254 #m^2/m^3 #hex_opt.fin_area_density_hot  
   beta_c      = 2254 #m^2/m^3 #hex_opt.fin_area_density_cold 
   
   #Efficiency 
   pump_efficiency  = 0.7 #hex_opt.pump_efficiency
   fan_efficiency   = 0.7 #hex_opt.fan_efficiency   
   
   # Calculate Surface geometry properties 
   
   # Assumes N passages for hot air and N+1 fpr cold air 
   N_p   = (H-b_c+2*delta_w)/(b_h+b_c+2*delta_w)
   
   # Frontal areas on hot and cokd sides
   
   A_fr_h = L_c*H
   A_fr_c = L_h*H
   
   # Heat exchnager volume between plates on each fluid side. 
   
   V_p_h       = L_h*L_c*b_h*N_p
   V_p_c       = L_h*L_c*b_c*(N_p+1)
   
   # The heat transfer areas 
   
   A_h         = beta_h*V_p_h
   A_c         = beta_c*V_p_c
   
   #The minimum free flow area
   
   A_o_h       = d_h_h*A_h/(4*L_h)     
   A_o_c       = d_h_c*A_c/(4*L_c)     
   
   # Minimum Free flow area 
   
   sigma_h    = A_o_h/A_fr_h
   sigma_c    = A_o_c/A_fr_c
   
   
   # Thermal Properties Evaluation 
   
   # EValuate input density Needed?
   
   rho_i_h  = P_i_h/(R*T_i_h)
   rho_i_c  = P_i_c/(R*T_i_c)
   
   # Mass flow rate use from sizing problem as input; hardcode for now 
   m_dot_h     = 1.66 #hex_opt.coolant_flow_rate
   m_dot_c     = 2.00 #hex_opt.air_flow_rate
   
   # ---------------------------------------------------------------------------------------------------------- 
   # Thermal Performance Calculation 
   # ----------------------------------------------------------------------------------------------------------       

   
   T_o_h       = T_i_h-eff_hex*(T_i_h-T_i_c)
   T_o_c       = T_i_c+eff_hex*(m_dot_h/m_dot_c)*(T_i_h-T_i_c)
   
   check             = 0
   itetation_counter = 0
   
   while check==0:
      
      T_m_h       = (T_i_h+T_o_h)/2
      T_m_c       = (T_i_c+T_o_c)/2
      
      #Prandtl Number 
      Pr_h    =  0.721 #coolant.compute_prandtl_number(T_m_h)
      Pr_c    =  0.692 #air.compute_prandtl_number(T_m_c)    
      
      #Absolute viscosity 
      mu_h    = 40.1e-6 #coolant.compute_absolute_viscosity(T_m_h)
      mu_c    = 33.6e-6 #air.compute_absolute_viscosity(T_m_c)       
      
      #Specific heat 
      c_p_h   = 1.122e3 #J/kg-K #coolant.compute_cp(T_m_h)/1000 #KJ/kg-K
      c_p_c   = 1.073e3 #J/kg-K #air.compute_cp(T_m_c)/1000     #KJ/kg-K   
      
      # Core mass velcoity 
      
      G_h   = m_dot_h/A_o_h
      G_c   = m_dot_c/A_o_c
      
      # Calculate Reynolds Number
      Re_h       = G_h * d_h_h / mu_h
      Re_c       = G_c * d_h_c / mu_c
      
      # Calculate the colburn factor and friction factor using curve fitted values (What about for Turbulent regim, check in london and Kays )
      j_c            = 0.0131 * (Re_c / 1000)**(-0.415)
      j_h            = 0.0131 * (Re_h / 1000)**(-0.415)
   
      f_c            = 0.0514 * (Re_c / 1000)**(-0.471)
      f_h            = 0.0514 * (Re_h / 1000)**(-0.471)
   
      j_f_h          = j_h/f_h
      j_f_c          = j_c/f_c
      
      # Heat Transfer Coefficients
      h_h = j_h * G_h * c_p_h / (Pr_h**(2/3))
      h_c = j_c * G_c * c_p_c / (Pr_c**(2/3))
      
      m_f_h = (np.sqrt((2*h_h)/(k_f*delta_h)))*np.sqrt(1+(delta_h/l_s_h))
      m_f_c = (np.sqrt((2*h_c)/(k_f*delta_c)))*np.sqrt(1+(delta_c/l_s_c))
   
   
      l_f_h = b_h / 2 - delta_h
      l_f_c = b_c / 2 - delta_c
   
      # Fin Efficiency
      eta_f_h = np.tanh(m_f_h * l_f_h) / (m_f_h * l_f_h)
      eta_f_c = np.tanh(m_f_c * l_f_c) / (m_f_c * l_f_c)
      
      # Overall Efficiency 
      eta_o_h = 1 - (1 - eta_f_h) * Af_A_h
      eta_o_c = 1 - (1 - eta_f_c) * Af_A_c
      
      # Wall ressistance 
      A_w   = L_c*L_h*(2*N_p+2)
      R_w   = delta_w/(k_w*A_w)
      
      # Calculate overall heat transfer without fouling
      UA = 1 / ((1 / (eta_o_h * h_h*A_h)) +R_w+ (1 / (eta_o_c* h_c*A_c)))
      
      # Heat Capcity 
      C_h            = m_dot_h*c_p_h
      C_c            = m_dot_c*c_p_c
      
      C_min, C_max   = min(C_h, C_c), max(C_h, C_c)
      C_r            = C_min / C_max  
      
      # NTU 
      NTU            = UA/C_min
      
      # Updated effectiveness and we neglect longitudnal conduction for now 
      eff_hex_updated        = (1 - np.exp(((NTU**0.22)/C_r)*(np.exp(-C_r*(NTU**(0.78))) - 1 ))) 
      
      # Heat trannsfer rate
      q= eff_hex*(T_i_h-T_i_c)*C_min
      
      # Updated Outlet temperatures 
      T_o_h_updated          = T_i_h-(q/C_h) 
      T_o_c_updated          = T_i_c+(q/C_c)
        
      
      residual = [abs(T_o_c-T_o_c_updated),abs(T_o_h-T_o_h_updated),abs(eff_hex-eff_hex_updated)]
         
      if all(element < 0.01 for element in residual):
         check = 1 
      else:
         itetation_counter +=1 
         T_o_c   = T_o_c_updated
         T_o_h   = T_o_h_updated
         eff_hex = eff_hex_updated
      
           
    # ---------------------------------------------------------------------------------------------------------- 
    # Pressure Drop Calculation 
    # ----------------------------------------------------------------------------------------------------------       
     
   #inital assumption 
   P_o_h   = P_i_h 
   P_o_c   = P_i_c  
   
   iteraion_counter_1 = 0
   check              = 0
   delta_p_c          = []
   delta_p_h          = []
   
   while check == 0:
           
      # from the inlet and outlet pressures given the mean density is calcualted. 
      rho_h_i  =  0.4751 #coolant.compute_density(T_i_h,P_i_h)
      rho_c_i  =  1.4726 #air.compute_density(T_i_c,P_i_c)    
      rho_h_o  =  0.9424 #coolant.compute_density(T_o_h_,P_o_h)
      rho_c_o  =  0.7126 #air.compute_density(T_o_c,P_o_c)    
      rho_h_m  = (rho_h_i+rho_h_o)/2
      rho_c_m  = (rho_c_i+rho_c_o)/2               
      
      # Kc_c, Ke_c     = compute_heat_exhanger_factors(kc_vals,ke_vals,sigma_c, Re_c) 
      # Need to check if the values obtained from the function are close to what is obtained ere 
      k_c_c = 0.36
      k_c_h = 0.36
      k_e_c = 0.42
      k_e_h = 0.42
      
      # Thermal Resistance on the hot and cold fluid sides
      
      R_h = 1 / (eta_o_h * h_h * A_h)
      R_c = 1 / (eta_o_c * h_c * A_c)
      
      # Compute Wall temperature
      T_w = (T_m_h + (R_h / R_c) * T_m_c) / (1 + R_h / R_c)
      
      # Considering temperature at wall effecting f value of 0.81 changes
      f_h_wall = f_h * np.power(((T_w + 273) / (273 + T_m_h)), 0.81)
      f_c_wall = f_c * np.power(((T_w + 273) / (273 + T_m_c)), 1)
      
      # Calculate Pressure Drop
      
      delta_p_c.append = np.power(G_c, 2) / (2 * rho_c_i) * ((1 - np.power(sigma_c, 2) + k_c_c)
                                                                   + 2 * (rho_c_i / rho_c_o - 1) + f_c_wall * 4 * L_c / d_h_c *
                                                                 rho_c_i / rho_c_m
                                                         - (1 - np.power(sigma_c, 2) - k_e_c) * rho_c_i / rho_c_o)
      
      delta_p_h.append = np.power(G_h, 2) / (2 * rho_h_i) * ((1 - np.power(sigma_h, 2) + k_c_h)
                                                                   + 2 * (rho_h_i / rho_h_o - 1) + f_h_wall * 4 * L_h / d_h_h *
                                                                    rho_h_i / rho_h_m
                                                         - (1 - np.power(sigma_h, 2) - k_e_h) * rho_h_i / rho_h_o)   
      
      if iteraion_counter_1 >=1:
         residual_pressure    = [abs(delta_p_c[iteraion_counter_1]-delta_p_c[iteraion_counter_1-1]),abs(delta_p_h[iteraion_counter_1]-delta_p_h[iteraion_counter_1-1])]
         
         if all(element<0.01 for element in residual_pressure):
            check =1
         else:
            P_o_c                         = (delta_p_c[iteraion_counter_1]+P_i_c)
            P_o_h                         = (delta_p_h[iteraion_counter_1]+P_i_h)
            iteraion_counter_1           += 1
      else:
         P_o_c                         = (delta_p_c[iteraion_counter_1]+P_i_c)
         P_o_h                         = (delta_p_h[iteraion_counter_1]+P_i_h)
         iteraion_counter_1           += 1
         
   # Caluclate Inlet and outler velocity
   u_i      =  np.sqrt((2*P_i_c)/rho_c_i)
   u_o      =  np.sqrt((2*P_o_c)/rho_c_o)   

        
   #Calculate Power drawn by HEX 
   P_hex    = ((m_dot_h*delta_p_h[iteraion_counter_1])/(pump_efficiency*rho_h_m))+((m_dot_c*delta_p_c[iteraion_counter_1]/rho_c_m)+(abs(u_i**2-u_o**2)/2))/(fan_efficiency)
                                       
   print('Power consumed bt the HEX is ', P_hex)
 
   return 
      
      
      
      
if __name__ == "__main__":
    
   main()
