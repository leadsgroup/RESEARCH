# imports  
import numpy as np           
import pandas as pd  

# Conceptual Design: Identification of a general description of the aircraft. The global shape and configuration is determined 
#                    based on the mission reuirements. 
#                    
#                    The characterization process is based on a small set of parameters:
#                    - MTOW = Maximum Take Off Weight/Design Gross Weight
#                    - P_b or T = Power (for Propeller driven aircraft) or Thrust (for Jet driven aircraft) -> W/P_b or T/W
#                    - S = Wing reference area -> W/S
#                    
#                    The conceptual design process is based on Mission Specification and derives from Market studies, 
#                    Manufacturing studies, Technological studies, Specific requirements
#                    The Mission Specification includes:
#                    - Mission Requirements
#                    - Certification Framework
#                    - Mission Profiles: Schematic description of the mission
#                    - Payload Requirements
#                    - Pressurization Requirements 
#
#                    The conceptual design phase Sizing to Performance looks at:
#                    - Performance Requirements: 
#                                  - Point performance: Stall Speed, Max speed, Economic cruise speed, manoeuvring speed, max RoC
#                                  - Intergral performance: Range, Endurance, Time to climb
#                                  - Field performance: TO distance, LND distance, Balance field length (BFL)
           

def main():
    Weight_Estimation()
    Performance_Requirements()    
    return

def Weight_Estimation():
    
#   Based on historical regressions for W_E (Empty) and Fuel Fractions method for W_F (Fuel)
#   W_TO = Take-Off weight
#   W_ZF = Zero Fuel Weight
#   W_F = Fuel Weight
#   W_EO = Empty Operative Weight
#   W_L = Payload Weight -> From Mission Requirements
#   W_TF = Trip Fuel Weight    
#   W_RF = Reserve Fuel Weight -> From Certification Requirements    
#   W_E = Empty Weight 
#   W_C = Crew Weight -> From Mission Requirements
#
#   W_TO = W_ZF + W_F
#         = (W_EO + W_L) + (W_TF + W_RF)
#         = (W_E + W_C + W_L) + (W_TF + W_RF)
#
#   Unknowns: W_E and W_TF
#
#   w_E = W_E/MTOW = f_E(MTOW)
#   w_F = W_F/MTOW = f_F(MTOW)
#   So: (1 - w_E - w_F)*MTOW = W_C + W_L
#
#   f_E(MTOW) found from historical regression, using statistical distribution 
#   Roskam: log_10(MTOW) = A + B*log_10(W_E) from diagram: log_10(W_E) vs log_10(MTOW)
#   Raymer: w_E = D*MTOW**C from diagram: w_E vs log_10(MTOW)
#
#   f_F(MTOW) found from Fuel Fractions method
#   w_F = W_F/MTOW = (MTOW - W_final)/MTOW = 1 - W_final/W_initial
#   Each mission phase has an initial and final weight (ex. 1, 2)
#   W_final/W_initial = (W_1/MTOW)*(W_2/W_1)*(W_3/W_2)*...*(W_N/W_(N-1))
#   For short missions, use statistical approach: ex. Warmup and TO: W_i/W_(i-1) = 0.97
#                                                     Climb: as cruise at C_Y
#                                                     LND: 0.995 
#   For long phases (Cruise, Loiter): Use Breguet Formulation
    
def Performance_Requirements():
#   Once the MTOW has been determined, we can determine the Project Point in the Sizing Matrix Plot
#   - Performance Requirements: 
#                 - Point performance: Stall Speed, Max speed, Economic cruise speed, manoeuvring speed, max RoC
#                 - Intergral performance: Range, Endurance, Time to climb
#                 - Field performance: TO distance, LND distance, Balance field length (BFL)

#   Point Performance:
#   - V_stall:
#   V_S = np.sqrt(2*W/(rho_0*S*C_L_max)) = f(W/S) <= V_S*
#   So: W/S <= (V_S*)**2/(2/(rho_0*C_L_max)) = (W/S)*
#
#   - LND distance: 
#   Use statistical approach
#   S_LND = k*V_SO**2 = f(W/S) <= S_LND*
#   So: W/S <= (W/S)*
#   
#   - TO distance: 
#   Use statistical approach
#   S_TO = k0*P_TO = f(W/S, W/P) <= S_TO*
#
#   - Maximum speed:
#   From P_req = P_av we obtain an expression in V**4 and V with two meaningful solutions, we pick the highest one
#   V_H =  f(W/S, W/P) >= V_H*
#
#   Climb Performance:
#   - Steepest Climb: max gamma 
#     sin(gamma) = (T-D)/W so -> V_SC = f(W/S, W/P)
#   - Fastest Climb: max V_v 
#     V_v = V*sin(gamma) -> V_FC = f(W/S, W/P)
#   - Ceiling: ex. V_v = 0
#   - Time to Climb: t = int(dh/V_v) -> t = f(W/S, W/P)
#
#   Turning flight:
#   - Instant turn: As in stall, V_s_n = np.sqrt(n)*V_stal -> f(W/S)
#   - Sustained turn -> f(W/S, W/P and n)
    
    
     
if __name__ == '__main__':
    main()