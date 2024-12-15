from RCAIDE.Framework.Core import Data

def main():

    # Input parameters
    input = {}
    input['Q_F']          = 16.2     # [cm**3/hr]        Jet-A feed rate
    input['Q_S']          = 60       # [cm**3/hr]        Deionized water feed rate
    input['Q_A']          = 600      # [sccm]            Air feed rate

    input['x_H']          = 0.1348   # [-]               mass fraction of hydrogen content in Jet-A
    input['x_C']          = 0.8637   # [-]               mass fraction of carbon content in Jet-A
    input['y_H2']         = 0.7      # [mol]          mole fraction of hydrogen content in reformate
    input['y_CO']         = 0.2     # [mol]          mole fraxtion of carbon monoxide content in reformate
    input['rho_F']        = 0.813    # [g/cm**3]         Density of Jet-A
    input['rho_S']        = 1        # [g/cm**3]         Density of water
    input['rho_A']        = 0.001293 # [g/cm**3]         Density of air
    input['MW_F']         = 160      # [g/g-mol]         Average molecular weight of Jet-A    
    input['MW_S']         = 18.01    # [g/g-mol]         Average molecular weight of steam
    input['MW_C']         = 12.01    # [g/g-mol]         Average molecular weight of carbon
    input['MW_H2']        = 2.016    # [g/g-mol]         Average molecular weight of hydrogen
    input['A_F_st_Jet_A'] = 14.62    # [lb_Air/lb_Jet_A] Stoichiometric air-to-fuel mass ratio 
    input['theta']        = 0.074    # [sec**-1]         Contact time
    input['LHV_F']        = 43.435   # [kJ/g-mol]        Lower heating value of Jet-A
    input['LHV_H2']       = 240.2    # [kJ/g-mol]        Lower heating value of Hydrogen
    input['LHV_CO']       = 283.1    # [kJ/g-mol]        Lower heating value of Carbon Monoxide
    input['V_cat']        = 9.653    # [cm**3]           Catalyst bed volume

    AutothermalReformer(input)

    return 

def AutothermalReformer(input):

    # Molar Feed Rates
    F_F = input['Q_F'] * input['rho_F'] / input['MW_F']                # [g-mol/hr] molar flow rate of Jet-A
    F_S = input['Q_S'] * input['rho_S'] / input['MW_S']                # [g-mol/hr] molar flow rate of steam
    F_A = (input['Q_A'] * 60) / 22414                                  # [g-mol/hr] molar flow rate of air
    F_C = input['Q_F'] * input['rho_F'] * input['x_C'] / input['MW_C'] # [g-mol/hr] molar flow rate of carbon

    # Effluent Gas Molar Flow Rate
    Q_R = (input['Q_F']/60) + (input['Q_S']/60)  + input['Q_A'] # [sccm] Reformer effluent gas feed rate
    F_R = Q_R * 60 / 22414                                  # [g-mol/hr] reformate effluent gas molar flow rate

    # Space Velocity
    GHSV = ((F_F + F_S + F_A) / input['V_cat']) * 22410 # [hr**-1] gas hourly space velocity
    LHSV = input['Q_F'] / input['V_cat']                # [hr**-1] liquid hourly space velocity

    # Steam to Carbon, Oxygen to Carbon and Equivalence Ratio 
    S_C = F_S / F_C                                                                                        # [mol_H20/mol_C] Steam-to-Carbon feed ratio
    O_C = 2 * 0.21 * F_A / F_C                                                                             # [mol_O/mol_C] Oxygen-to-Carbon feed ratio
    phi = input['A_F_st_Jet_A'] * (input['Q_F'] * input['rho_F']) / ((input['Q_A'] * 60) * input['rho_A']) # [-] Fuel to Air ratio

    # Reformer efficiency
    eta_ref = ((input['y_H2'] * input['LHV_H2'] + input['y_CO'] * input['LHV_CO']) * F_R / (input['Q_F'] * input['rho_F'] * input['LHV_F'])) * 100 # [-] Reformer efficiency

    # Hydrogen conversion efficiency
    X_H2 = ((input['y_H2'] * F_R)/ (((input['Q_F'] * input['rho_F'] * input['x_H'])/(input['MW_H2'])) + F_S)) * 100 # [-] Hydrogen conversion efficiency  

    return    

if __name__ == "__main__":
    main()