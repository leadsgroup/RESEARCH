## @ingroup Components-Energy-Cooling
# Cryogenic_Heat_Exchanger.py
#
# Created:  Oct. 2022,  C.R. Zhao

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE

import numpy as np
from scipy.optimize import fsolve

from SUAVE.Core import Data
from SUAVE.Components.Energy.Energy_Component import Energy_Component
from SUAVE.Attributes.Coolants.Glycol_Water import Glycol_Water

# ----------------------------------------------------------------------
#  Cryogenic Heat Exchanger Component
# ----------------------------------------------------------------------
## @ingroup Components-Energy-Cooling

class WavyChan_GasLiqHEX(Energy_Component):
    """This provides output values for a heat exchanger used to cool components
    
    Assumptions:
    None
    
    Source:
    N/A
    """
    
    def __defaults__(self):
        """This sets the default values for the component to function.

        Assumptions:
        None

        Source:
        source: Zhao C, Sousa A C M, Jiang F. Minimization of thermal non-uniformity in lithium-ion battery pack
        cooled by channeled liquid flow[J]. International journal of heat and mass transfer, 2019, 129: 660-670.

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        None
        """         
        
        self.tag = 'WavyChan_GasLiqHEX'

        self.wavychan = Data()
        self.wavychan.cross_section = Data()
        self.hex      = Data()
        # -----setting the default values for the different components
        # wavy channel: thermophysical properties
        # parameters                                    values                              units
        self.wavychan.density                           = 2719                              # kg/m^3
        self.wavychan.thermal_conductivity              = 202.4                             # W/m.K
        self.wavychan.specific_heat_capacity            = 871                               # J/kg.K
        # wavy channel: geometric properties
        self.wavychan.cross_section.a                   = 1e-3                              # m
        self.wavychan.cross_section.b                   = 5e-4                              # m
        self.wavychan.cross_section.c                   = 6.3e-2                            # m
        self.wavychan.cross_section.d                   = 2e-3                              # m
        self.wavychan.contact_angle                     = 47.5                              # ^o

        # heat exchanger: thermophysical properties
        self.hex.density                                = 2780                              # kg/m^3
        self.hex.thermal_conductivity                   = 121                               # W/m.K
        self.hex.specific_heat_capacity                 = 871                               # J/kg.K
        # heat exchanger: geometric properties
        # other parameters
        self.hex.t_w                                    = 5e-4                              # m
        self.hex.t_f                                    = 1e-4                              # m
        self.hex.k_f                                    = 121                               # W/m.K
        self.hex.k_w                                    = 121                               # W/m.K

        self.pump_efficiency                            = 0.7
        self.fan_efficiency                             = 0.7

        # generic parameters required
        self.wavychan.Reynolds_number_of_coolant_in_chan            = 1.
        self.wavychan.velocity_of_coolant_in_chan                   = 1.
        self.wavychan.contact_area_of_chan_with_batteries           = 1.
        self.wavychan.hydraulic_diameter_of_chan                    = 1.
        self.wavychan.aspect_ratio_of_chan                          = 1.
        self.wavychan.mass_flow_rate_of_one_module                  = 1.
        self.wavychan.contact_area_of_one_module                    = 1.

        self.wavychan.mass_of_chan                                  = 1.
        self.wavychan.power_of_chan                                 = 1.

        self.hex.fin_height_of_hot_fluid                            = 1.
        self.hex.fin_height_of_cold_fluid                           = 1.
        self.hex.finned_area_by_overall_area_at_hot_fluid           = 1.
        self.hex.finned_area_by_overall_area_at_cold_fluid          = 1.
        self.hex.area_density_of_hot_fluid                          = 1.
        self.hex.area_density_of_cold_fluid                         = 1.
        self.hex.number_of_plates                                   = 1.
        self.hex.frontal_area_of_hot_fluid                          = 1.
        self.hex.frontal_area_of_cold_fluid                         = 1.
        self.hex.heat_transfer_area_of_hot_fluid                    = 1.
        self.hex.heat_transfer_area_of_cold_fluid                   = 1.
        self.hex.minimum_free_flow_area_of_hot_fluid                = 1.
        self.hex.minimum_free_flow_area_of_cold_fluid               = 1.
        self.hex.ratio_of_free_flow_area_to_frontal_area_of_hot_fluid   = 1.
        self.hex.ratio_of_free_flow_area_to_frontal_area_of_cold_fluid  = 1.

        self.hex.mass_of_hex                                        = 1.
        self.hex.power_of_hex_fraction_of_hot_fluid                 = 1.

        self.hex.length_of_hot_fluid                        = 1.
        self.hex.length_of_cold_fluid                       = 1.
        self.hex.stack_height                               = 1.
        self.hex.hydraulic_diameter_of_hot_fluid_channel    = 1.
        self.hex.aspect_ratio_of_hot_fluid_channel          = 1.
        self.hex.hydraulic_diameter_of_cold_fluid_channel   = 1.
        self.hex.aspect_ratio_of_cold_fluid_channel         = 1.

        self.wavychan.channel_numbers_of_module             = 1.

        self.mass_flow_rate_of_hot_fluid                    = 1.
        self.hex.mass_flow_rate_of_cold_fluid               = 1.
        self.hex.pressure_ratio_of_hot_fluid                = 1.
        self.hex.pressure_ratio_of_cold_fluid               = 1.
        self.hex.pressure_at_inlet_of_hot_fluid             = 1.
        self.hex.inlet_temperature_of_hot_fluid             = 1.

        # operating conditions
        self.operating_conditions = Data()
        self.operating_conditions.wavy_channel_inlet_temperature                    = 288.
        self.operating_conditions.wavy_channel_outlet_temperature                   = 1.
        self.operating_conditions.heat_exchanger_inlet_temperature_of_hot_fluid     = 1.
        self.operating_conditions.heat_exchanger_outlet_temperature_of_hot_fluid    = 1.
        self.operating_conditions.heat_exchanger_inlet_temperature_of_cold_fluid    = 1.
        self.operating_conditions.heat_exchanger_outlet_temperature_of_cold_fluid   = 1.
        self.operating_conditions.heat_extraction_from_battery_to_wavy_channel      = 1.
        self.operating_conditions.heat_dissipation_to_ambient_from_heat_exchanger   = 1.
        self.operating_conditions.heat_exchanger_effectiveness                      = 1.
        self.operating_conditions.wavy_channel_effectiveness                        = 1.
        self.operating_conditions.wavy_channel_power                                = 1.
        self.operating_conditions.heat_exchanger_power                              = 1.
        self.operating_conditions.battery_current_temperature                       = 1.

        return

    def energy_calc(self, battery, total_cells, single_side_contact=True, dry_mass=True):

        """ This calculates the mass of wavychannel and heat exchanger

        Assumptions:
        None.

        Source:
        N/A

        Inputs:
        self.inputs
            cooling_power      [W]

        Outputs:
        self.outputs.
            mdot                   [kg/s]

        Properties Used:
        self.
            cryogen_inlet_temperature       [K]
            cryogen_outlet_temperature      [K]
            cryogen_pressure                [Pa]

        
        """
        Coolant                             = Glycol_Water()
        density_coolant                     = Coolant.density
        kinematic_viscosity_coolant         = Coolant.dynamic_viscosity / Coolant.density

        # unpack the values from self
        D_cell                   = battery.cell.diameter
        H_cell                   = battery.cell.height
        total_cells              = total_cells
        number_of_modules        = battery.module_config.number_of_modules

        btms                     = self

        a                        = btms.wavychan.cross_section.a
        b                        = btms.wavychan.cross_section.b
        c                        = btms.wavychan.cross_section.c
        d                        = btms.wavychan.cross_section.d
        theta                    = btms.wavychan.contact_angle
        density_chan             = btms.wavychan.density

        density_hex             = btms.hex.density
        t_w                     = btms.hex.t_w
        t_f                     = btms.hex.t_f
        k_f                     = btms.hex.k_f
        k_w                     = btms.hex.k_w

        pump_efficiency          = btms.pump_efficiency
        fan_efficiency           = btms.fan_efficiency

        # ---------------------------------------------------------------------------------
        # Input parameters
        # ---------------------------------------------------------------------------------
        # heat exchanger: overall dimensions
        L_h                 = btms.hex.length_of_hot_fluid
        L_c                 = btms.hex.length_of_cold_fluid
        H                   = btms.hex.stack_height
        # heat exchanger: geometric dimensions
        d_H_h               = btms.hex.hydraulic_diameter_of_hot_fluid_channel
        AP_chan_h           = btms.hex.aspect_ratio_of_hot_fluid_channel
        d_H_c               = btms.hex.hydraulic_diameter_of_cold_fluid_channel
        AP_chan_c           = btms.hex.aspect_ratio_of_cold_fluid_channel
        # wavy channel
        channel_numbers_per_module = btms.wavychan.channel_numbers_of_module
        # operating conditions
        m_dot_h             = btms.mass_flow_rate_of_hot_fluid
        m_dot_c             = btms.hex.mass_flow_rate_of_cold_fluid
        R_p_h               = btms.hex.pressure_ratio_of_hot_fluid
        R_p_c               = btms.hex.pressure_ratio_of_cold_fluid

        p_h_1               = btms.hex.pressure_at_inlet_of_hot_fluid

        # ---------------------------------------------------------------------------------
        # Compute Mass and Power: wavy channel
        # ---------------------------------------------------------------------------------
        # contact area/length of a single battery cell
        if single_side_contact:
            A_surf          = theta / 360 * np.pi * D_cell * H_cell
            L_surf          = theta / 360 * np.pi * (D_cell + b + d / 2)
        else:
            A_surf          = 2 * theta / 360 * np.pi * D_cell * H_cell
            L_surf          = 2 * theta / 360 * np.pi * (D_cell + b + d / 2)

        # contact area/length of the battery pack
        if single_side_contact:
            A_chan          = total_cells * A_surf
            L_chan          = total_cells * L_surf + 4 * D_cell
        else:
            A_chan          = 2 * total_cells * A_surf
            L_chan          = 2 * (total_cells * L_surf + 4 * D_cell)

        # line density of the wavy channel
        if dry_mass:
            rho_line        = density_chan * (2 * a * (2 * b + d) + 2 * b * c)
        else:
            rho_line        = density_chan * (2 * a * (2 * b + d) + 2 * b * c) + density_coolant * c * d

        # MASS [kg]
        m_CHAN              = rho_line * L_chan

        # GEOMETRIC PARAMETERS
        d_H_chan            = 2 * (c * d) / (c + d)
        AP_chan_wavy        = d / c
        Ac_chan             = c * d

        # calculate coolant velocity
        number_of_channels  = number_of_modules * channel_numbers_per_module
        m_dot_h_unit        = m_dot_h / number_of_channels
        u_coolant           = m_dot_h_unit / (density_coolant * Ac_chan)

        # calculate Fanning friction factor
        Re_coolant = u_coolant * d_H_chan / kinematic_viscosity_coolant

        if Re_coolant <= 2300:
            # Fanning friction factor; c_geom constant: f * Re
            f_coolant = 24 * (1 - 1.3553 * AP_chan_wavy + 1.9467 * np.power(AP_chan_wavy, 2) - 1.7012 * np.power(AP_chan_wavy, 3)
                           + 0.9564 * np.power(AP_chan_wavy, 4) - 0.2537 * np.power(AP_chan_wavy, 5)) / Re_coolant
        else:
            # Fanning friction factor
            f_coolant = 1 / (4 * (1.8 * np.log10(Re_coolant / 7.7)) ** 2)

        # pressure drop
        delta_P_chan = 2 * f_coolant * (density_coolant * u_coolant * u_coolant) * (L_chan / d_H_chan)

        # POWER  [kW]
        P_CHAN              = m_dot_h * delta_P_chan / (1e3 * density_coolant * pump_efficiency)

        # ---------------------------------------------------------------------------------
        # Compute Mass and Power: heat exchanger
        # ---------------------------------------------------------------------------------
        """geometrical parameters"""
        # fin height
        b_h                             = d_H_h * (1 + AP_chan_h) / 2
        b_c                             = d_H_c * (1 + AP_chan_c) / 2

        # finned area by overall area
        A_f_by_A_h                      = AP_chan_h / (AP_chan_h + 1)
        A_f_by_A_c                      = AP_chan_c / (AP_chan_c + 1)

        # area density
        beta_h                          = 4 * (1 + AP_chan_h) / (d_H_h * (1 + AP_chan_h) + 2 * AP_chan_h * t_f)
        beta_c                          = 4 * (1 + AP_chan_c) / (d_H_c * (1 + AP_chan_c) + 2 * AP_chan_c * t_f)

        # assume  N_p passages for the liquid and (N_p + 1) passages for the air
        N_p = (H - b_c - 2 * t_w) / (b_h + b_c + 2 * t_w)
        # N_p = np.ceil(N_p)

        # the frontal areas on both fluid sides
        A_fr_h                          = L_c * H
        A_fr_c                          = L_h * H

        # the heat exchanger volume between plates on both fluid sides
        V_p_h                           = L_h * L_c * N_p * b_h
        V_p_c                           = L_h * L_c * (N_p + 1) * b_c

        # the heat transfer areas
        A_h                             = beta_h * V_p_h
        A_c                             = beta_c * V_p_c

        # the minimum free-flow areas [cross-section area of the fluid passages]
        A_o_h                           = d_H_h * A_h / (4 * L_h)
        A_o_c                           = d_H_c * A_c / (4 * L_c)

        # the ratio of minimum free-flow area to the frontal area
        sigma_h                         = A_o_h / A_fr_h
        sigma_c                         = A_o_c / A_fr_c

        # MASS [kg]
        if dry_mass:
            m_HEX                       = density_hex * L_h * A_fr_h * (1 - sigma_h - sigma_c)
        else:
            m_HEX                       = density_hex * L_h * A_fr_h * (1 - sigma_h - sigma_c) + density_coolant * A_o_h * L_h

        # POWER [kW]
        delta_p_h                       = p_h_1 * (1 - R_p_h)
        P_HEX_h                         = m_dot_h * (delta_p_h / density_coolant) / pump_efficiency / 1e3
        # ---------------------------------------------------------------------------------
        # Pack outputs
        # ---------------------------------------------------------------------------------
        btms.wavychan.Reynolds_number_of_coolant_in_chan            = Re_coolant
        btms.wavychan.velocity_of_coolant_in_chan                   = u_coolant
        btms.wavychan.contact_area_of_chan_with_batteries           = A_chan
        btms.wavychan.hydraulic_diameter_of_chan                    = d_H_chan
        btms.wavychan.aspect_ratio_of_chan                          = AP_chan_wavy
        btms.wavychan.mass_flow_rate_of_one_module                  = m_dot_h_unit * channel_numbers_per_module
        btms.wavychan.contact_area_of_one_module                    = A_chan / battery.module_config.number_of_modules

        btms.wavychan.mass_of_chan                                  = m_CHAN
        btms.wavychan.power_of_chan                                 = P_CHAN

        btms.hex.fin_height_of_hot_fluid                            = b_h
        btms.hex.fin_height_of_cold_fluid                           = b_c
        btms.hex.finned_area_by_overall_area_at_hot_fluid           = A_f_by_A_h
        btms.hex.finned_area_by_overall_area_at_cold_fluid          = A_f_by_A_c
        btms.hex.area_density_of_hot_fluid                          = beta_h
        btms.hex.area_density_of_cold_fluid                         = beta_c
        btms.hex.number_of_plates                                   = N_p
        btms.hex.frontal_area_of_hot_fluid                          = A_fr_h
        btms.hex.frontal_area_of_cold_fluid                         = A_fr_c
        btms.hex.heat_transfer_area_of_hot_fluid                    = A_h
        btms.hex.heat_transfer_area_of_cold_fluid                   = A_c
        btms.hex.minimum_free_flow_area_of_hot_fluid                = A_o_h
        btms.hex.minimum_free_flow_area_of_cold_fluid               = A_o_c
        btms.hex.ratio_of_free_flow_area_to_frontal_area_of_hot_fluid   = sigma_h
        btms.hex.ratio_of_free_flow_area_to_frontal_area_of_cold_fluid  = sigma_c

        btms.hex.mass_of_hex                                        = m_HEX
        btms.hex.power_of_hex_fraction_of_hot_fluid                 = P_HEX_h

        return btms
        
