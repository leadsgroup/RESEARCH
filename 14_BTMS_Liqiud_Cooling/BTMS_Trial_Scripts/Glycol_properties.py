## @ingroup Attributes-Coolants
# Glycol_Water
#
# Created:  March 2024,  S.S Shekar

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
from .Coolant import Coolant

# ----------------------------------------------------------------------
#  Liquid H2 Cryogen Class
# ----------------------------------------------------------------------
## @ingroup Attributes-Coolants
class Glycol_Water(Coolant):
    """     """

    def __defaults__(self):
        """This sets the default values.

        Assumptions:
        We assume the mixture is 50%

        Source:

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        None
        """

        self.tag                       = 'Glycol_Water'
        self.percent_glycol            = 0.5 
        self.density                   = 1075                       # kg/m^3
        self.specific_heat_capacity    = 3300                       # J/kg.K
        self.thermal_conductivity      = 0.387                      # W/m.K
        self.dynamic_viscosity         = 0.0019                     # Pa.s
        self.Prandtl_number            = self.specific_heat_capacity * self.dynamic_viscosity / self.thermal_conductivity
        self.kinematic_viscosity       = self.dynamic_viscosity / self.density

    def compute_cp(self,T=300):
        self.specific_heat_capacity = 1e-11 * T**3 + 7e-9 * T**2 + 0.0007 * T + 0.6151
        
        return self.specific_heat_capacity
    
    def compute_absolute_viscosity(self,T=300):
        self.dynamic_viscosity  = 9e20 * T**(- 9.451)

        return self.dynamic_viscosity
   
    def compute_density(self,T=300.): 
        self.density = 0.00001*T**3 - 0.0136*T**2 + 3.9013*T + 760.75
        return self.density  
    
    def compute_thermal_conductivity(self,T=300.): 
        self.thermal_conductivity = 7e-7 * T**2 + 0.0007 * T + 0.0593
        return self.thermal_conductivity
    
    
    def compute_prandtl_number(self,T=300.): 
        Cp = self.compute_cp(T)
        mu = self.compute_dynamic_viscosity(T)
        K  = self.compute_thermal_conductivity(T)
        return  mu*Cp/K      