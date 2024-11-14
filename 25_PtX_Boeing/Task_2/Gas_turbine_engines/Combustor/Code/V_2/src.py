from typing import Any
import numpy as np
import cantera as ct


class PSR(object):

    def __init__(
        self,
        TP: list,
        phi: float,
        length: float,
        area: float,        
        fuel_file: str,
        fuel_species: dict,
        oxidizer_species: dict,
        mean_phi: float,
        sigma_phi: float,
        delta_phi: float,
        tot_mdot_air: float,
        tot_mdot_fuel: float
    ) -> None:
        self.TP = list(TP)
        self.phi = float(phi),        
        self.length = float(length)
        self.area = float(area)        
        self.fuel_file = str(fuel_file)
        self.fuel_species = dict(fuel_species)
        self.oxidizer_species = dict(oxidizer_species)
        self.mean_phi = float(mean_phi)
        self.sigma_phi = float(sigma_phi)
        self.delta_phi = float(delta_phi)
        self.tot_mdot_air = float(tot_mdot_air)
        self.tot_mdot_fuel = float(tot_mdot_fuel)

        self._init_fuel()
        self._set_gaussian()
        self._set_mdot()

        self.reactor = ct.IdealGasReactor(self.fuel)
        self.reactor.volume = self.length*self.area

        self.upstream = ct.Reservoir(self.fuel)
        self.inlet = ct.MassFlowController(self.upstream, self.reactor)
        self.inlet.mass_flow_rate = self.mdot

        self.rho = self.fuel.density
        self.t_res = self.rho * self.reactor.volume / self.mdot

        self.reacnet = ct.ReactorNet([self.reactor])

    def _init_fuel(self):
        self.fuel = ct.Solution(self.fuel_file)
        self.fuel.TP = self.TP
        self.fuel.set_equivalence_ratio(
            self.phi,
            fuel=self.fuel_species,
            oxidizer=self.oxidizer_species
        )
        self.fuel.equilibrate('HP')

    def _set_gaussian(self):
        self.f = (1 / (np.sqrt(2 * np.pi) * self.sigma_phi)) \
               * np.exp((-(self.phi - self.mean_phi) ** 2) \
               / (2 * self.sigma_phi ** 2)) * self.delta_phi

    def _set_mdot(self):
        self.mdot_air = self.tot_mdot_air * self.f
        self.mdot_fuel = self.tot_mdot_fuel * self.f
        self.mdot = self.mdot_air + self.mdot_fuel

    def mix(self, mixer):
        self.outlet = ct.MassFlowController(self.reactor, mixer, mdot=self.mdot)

    def integrate(self):
        self.reacnet.advance(self.t_res)


class Combustor(object):

    def __init__(
        self
    ) -> None:
        pass

    
        self.psrs = []
        self.mixers = []
        for i in range(combustor_input['N_PZ']):
            self.psrs.append(PSR(
                TP=[combustor_input['T_stag_0'], combustor_input['P_stag_0']],
                phi=phi_PSR[i],
                length=combustor_input['L_PZ'],
                area=combustor_input['A_PZ'],
                fuel_file=fuel_file,
                fuel_species=dict_fuel,
                oxidizer_species=dict_oxy,
                mean_phi=phi_sign,
                sigma_phi=sigma_phi,
                delta_phi=Delta_phi,
                tot_mdot_air=m_dot_air_PSR,
                tot_mdot_fuel=m_dot_fuel_PSR
            ))
            if (i%2 == 1):
                mixer = ct.IdealGasReactor(psrs[-2].fuel)
                self.mixers.append(mixer)
                self.psrs[-2].mix(mixer)
                self.psrs[-1].mix(mixer)





        

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        pass




# def combustor(
#     high_fidelity_kin_mech,
#     dict_fuel,
#     dict_oxy,
#     T_stag_0, P_stag_0, FAR, FAR_TO,FAR_st, m_dot_fuel, m_dot_fuel_TO, m_dot_air_id, m_dot_air_TO, N_PZ,V_PZ, phi_PZ_des, S_PZ, phi_SZ_des_1, phi_SZ_des_2, phi_SZ_des_3,F_SC, A_SZ, L_SZ):

#     # ------------------------------------------------------------------------------
#     # ----------------------------- Initial Parameters -----------------------------
#     # ------------------------------------------------------------------------------

#     f_air_PZ                = (m_dot_fuel_TO*F_SC)/(phi_PZ_des*m_dot_air_TO*FAR_st) # [-]       Fraction of total air present in the combustor that enters the Primary Zone
#     f_air_SZ                = 1 - f_air_PZ                                          # [-]       Fraction of total air present in the combustor that enters the Secondary Zone
#     m_dot_air_PZ            = f_air_PZ*m_dot_air_id                                 # [kg/s]    Air mass flow going through the Primary Zone
#     m_dot_air_SZ            = (f_air_SZ*m_dot_air_id)/3                             # [kg/s]    Air mass flow going through each dilution air inlet (3 inlets)
#     phi_sign                = ((m_dot_fuel*F_SC)/m_dot_air_PZ)/(FAR_st)             # [-]       Mean Equivalence Ratio
#     sigma_phi               = S_PZ*phi_sign                                         # [-]       Standard deviation of the Equivalence Ratio
#     m_dot_air_PSR           = m_dot_air_PZ                                          # [kg/s]    Air mass flow going through each PSR
#     m_dot_fuel_PSR          = m_dot_fuel                                            # [kg/s]    Fuel mass flow going through each PSR
#     V_PZ_PSR                = V_PZ/N_PZ                                             # [m**3]    Volume of each PSR
#     phi_PSR                 = np.linspace(0.001, 2*phi_sign, N_PZ)                  # [-]       Distribution of Equivalence Ratio through the PSRs
#     Delta_phi               = np.abs(phi_PSR[0] - phi_PSR[1])                       # [-]       Difference between two subsequent Equivalence Ratios
#     comp_fuel               = list(dict_fuel.keys())                                # [-]       Fuel components

#     if high_fidelity_kin_mech:
#         fuel_file              = ct.Solution('Jet_A_High_Fidelity.yaml')               # [-]       Import full fuel kinematic mechanism
#     else:
#         fuel_file              = ct.Solution('Jet_A_Low_Fidelity.yaml')                # [-]       Import surrogate fuel kinematic mechanism
    

#     nb_psr = 8
#     if (nb_psr%2 == 1):
#         raise ValueError("The number of PSR must be even.")

def combustor(
    high_fidelity_kin_mech,
    dict_fuel,
    dict_oxy,
    T_stag_0, 
    P_stag_0, 
    FAR, 
    FAR_TO,
    FAR_st, 
    m_dot_fuel, 
    m_dot_fuel_TO, 
    m_dot_air_id, 
    m_dot_air_TO, 
    N_PZ,
    V_PZ, 
    phi_PZ_des, 
    S_PZ, 
    phi_SZ_des_1, 
    phi_SZ_des_2, 
    phi_SZ_des_3,
    F_SC, 
    A_SZ, 
    L_SZ
    ):

    # ------------------------------------------------------------------------------
    # ----------------------------- Initial Parameters -----------------------------
    # ------------------------------------------------------------------------------

    f_air_PZ                = (m_dot_fuel_TO*F_SC)/(phi_PZ_des*m_dot_air_TO*FAR_st) # [-]       Fraction of total air present in the combustor that enters the Primary Zone
    f_air_SZ                = 1 - f_air_PZ                                          # [-]       Fraction of total air present in the combustor that enters the Secondary Zone
    m_dot_air_PZ            = f_air_PZ*m_dot_air_id                                 # [kg/s]    Air mass flow going through the Primary Zone
    m_dot_air_SZ            = (f_air_SZ*m_dot_air_id)/3                             # [kg/s]    Air mass flow going through each dilution air inlet (3 inlets)
    phi_sign                = ((m_dot_fuel*F_SC)/m_dot_air_PZ)/(FAR_st)             # [-]       Mean Equivalence Ratio
    sigma_phi               = S_PZ*phi_sign                                         # [-]       Standard deviation of the Equivalence Ratio
    m_dot_air_PSR           = m_dot_air_PZ                                          # [kg/s]    Air mass flow going through each PSR
    m_dot_fuel_PSR          = m_dot_fuel                                            # [kg/s]    Fuel mass flow going through each PSR
    V_PZ_PSR                = V_PZ/N_PZ                                             # [m**3]    Volume of each PSR
    phi_PSR                 = np.linspace(0.001, 2*phi_sign, N_PZ)                  # [-]       Distribution of Equivalence Ratio through the PSRs
    Delta_phi               = np.abs(phi_PSR[0] - phi_PSR[1])                       # [-]       Difference between two subsequent Equivalence Ratios
    comp_fuel               = list(dict_fuel.keys())                                # [-]       Fuel components

    if high_fidelity_kin_mech:
        fuel_file              = ct.Solution('Jet_A_High_Fidelity.yaml')               # [-]       Import full fuel kinematic mechanism
    else:
        fuel_file              = ct.Solution('Jet_A_Low_Fidelity.yaml')                # [-]       Import surrogate fuel kinematic mechanism
    

    N_PZ = 8
    if (N_PZ%2 == 1):
        raise ValueError("The number of PSR must be even.")

    psrs = []
    mixers = []
    for i in range(N_PZ):
        psrs.append(PSR(
            TP=[T_stag_0, P_stag_0],
            phi=phi_PSR[i],
            volume=V_PZ_PSR,
            fuel_file=fuel_file,
            fuel_species=dict_fuel,
            oxidizer_species=dict_oxy,
            mean_phi=phi_sign,
            sigma_phi=sigma_phi,
            delta_phi=Delta_phi,
            tot_mdot_air=m_dot_air_PSR,
            tot_mdot_fuel=m_dot_fuel_PSR
        ))
        if (i%2 == 1):
            mixer = ct.IdealGasReactor(psrs[-2].fuel)
            mixers.append(mixer)
            psrs[-2].mix(mixer)
            psrs[-1].mix(mixer)
