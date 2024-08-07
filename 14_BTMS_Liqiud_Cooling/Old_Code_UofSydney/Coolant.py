## @ingroup Attributes-Coolants
# Coolant.py
#
# Created:  Dec. 2022  C.R. Zhao

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from RCAIDE.Core import Data

# ----------------------------------------------------------------------
#  Class
# ----------------------------------------------------------------------
## @ingroup Attributes-Coolants
class Coolant(Data):
    """Holds values for a coolant liquid

    Assumptions:
    None

    Source:
    None
    """

    def __defaults__(self):
        """This sets the default values.

        Assumptions:
        None

        Source:
        Values commonly available

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        None
        """
        self.tag                       = 'Coolant'
        self.density                   = 0.0                       # kg/m^3
        self.specific_heat_capacity    = 0.0                       # J/kg.K
        self.thermal_conductivity      = 0.0                       # W/m.K
        self.dynamic_viscosity         = 0.0                       # Pa.s
        self.temperatures              = Data()
