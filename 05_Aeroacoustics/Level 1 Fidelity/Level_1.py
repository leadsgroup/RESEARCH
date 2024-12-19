## @ingroup Methods-Noise-Multi_Fidelity
# RCAIDE/Methods/Noise/Multi_Fidelity/Level_1.py
# 
# 
# Created:  Apr 2024, Niranjan Nanjappa

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports    
from RCAIDE.Framework.Core    import Data 
from RCAIDE.Library.Methods.Aerodynamics.Airfoil_Panel_Method    import airfoil_analysis
from RCAIDE.Library.Methods.Geometry.Two_Dimensional.Airfoil     import compute_naca_4series
from RCAIDE.Library.Methods.Geometry.Two_Dimensional.Airfoil     import import_airfoil_geometry

# package imports  
import numpy as np
import scipy as sp
from scipy.special import jv

# ----------------------------------------------------------------------------------------------------------------------
# Level_1
# ----------------------------------------------------------------------------------------------------------------------   
## @ingroup Methods-Noise-Multi_Fidelity
def harmonic_noise_l1(harmonics,freestream,angle_of_attack,coordinates,velocity_vector,rotor,aeroacoustic_data,settings,res):
    '''This computes the  harmonic noise (i.e. thickness and loading noise) of a rotor or rotor
    in the frequency domain

    Assumptions:
    Compactness of thrust and torque along blade radius from root to tip
    Thin airfoil assumption

    Source:
    1) Hanson, Donald B. "Helicoidal surface theory for harmonic noise of rotors in the far field."
    AIAA Journal 18.10 (1980): 1213-1220.

    2) Hubbard, Harvey H., ed. Aeroacoustics of flight vehicles: theory and practice. Vol. 1.
    NASA Office of Management, Scientific and Technical Information Program, 1991.


    Inputs: 
        harmonics                     - harmomics                                                                  [Unitless]
        freestream                    - freestream data structure                                                  [m/s]
        angle_of_attack               - aircraft angle of attack                                                   [rad]
        position_vector               - position vector of aircraft                                                [m]
        velocity_vector               - velocity vector of aircraft                                                [m/s] 
        rotors                        - data structure of rotors                                                   [None]
        aeroacoustic_data             - data structure of acoustic data                                            [None]
        settings                      - accoustic settings                                                         [None] 
        res                           - results data structure                                                     [None] 

    Outputs 
        res.                                    *acoustic data is stored and passed in data structures*                                                                            
            SPL_prop_harmonic_bpf_spectrum       - harmonic noise in blade passing frequency spectrum              [dB]
            SPL_prop_harmonic_bpf_spectrum_dBA   - dBA-Weighted harmonic noise in blade passing frequency spectrum [dbA]                  
            SPL_prop_harmonic_1_3_spectrum       - harmonic noise in 1/3 octave spectrum                           [dB]
            SPL_prop_harmonic_1_3_spectrum_dBA   - dBA-Weighted harmonic noise in 1/3 octave spectrum              [dBA] 
            p_pref_harmonic                      - pressure ratio of harmonic noise                                [Unitless]
            p_pref_harmonic_dBA                  - pressure ratio of dBA-weighted harmonic noise                   [Unitless]

    '''
    
    num_h        = len(harmonics)     
    num_cpt      = len(angle_of_attack) 
    num_mic      = len(coordinates.X_hub[0,:,0,0,0,0])
    num_rot      = len(coordinates.X_hub[0,0,:,0,0,0]) 
    phi_0        = np.array([rotor.phase_offset_angle])  # phase angle offset  
    num_sec      = len(rotor.radius_distribution) 
    orientation  = np.array(rotor.orientation_euler_angles) * 1 
    body2thrust  = sp.spatial.transform.Rotation.from_rotvec(orientation).as_matrix()
    
    
    # Omega
    Omega = 2*np.pi*rpm
    
    # span coordinate
    z   = np.linspace(0, rT, Nr)[1:]
    BD  = chord/D
    
    # distances
    R  = np.sqrt(np.sum((S_L-R_L)**2))
    r  = np.sqrt(np.sum(R_L**2))
    
    
    # source speed
    S_Speed = np.sqrt(np.sum(S_V**2))
    
    
    # angle between source velocity and reciever location (retarded)
    cosTheta  = ((S_L - R_L)*S_V)/(R*S_Speed)
    Theta     = np.arccos(cosTheta)
    sinTheta  = np.sin(Theta)
    
    
    # mach numbers
    # source mach number
    Mx = S_Speed/c0

    # tip rotational mach number            
    Mt = Omega*rT/c0
    
    # section relative mach number
    Mr = np.sqrt(Mx**2 + (z**2)*Mt**2)
    
    # OmegaD
    OmegaD = Omega/(1 - Mx*cosTheta)
    
    
    # phase shifts
    
    
    # dimensionless chordwise wavenumber
    kx = (2*m*B*BD*Mt)/(Mr*(1 - Mx*cosTheta))
    
    # ratio of maximum thickness to chord
    
    
    # Lift and Drag coefficients
    
    
    # wavenumber normal to chord
    ky = (2*m*B*BD/(z*Mr))*((Mr**2*cosTheta - Mx)/(1 - Mx*cosTheta))
    
    
    
# the code is under development
# https://www.grants.gov/search-results-detail/352892
# https://www.grants.gov/search-results-detail/352583
# https://www.grants.gov/search-results-detail/353777
# https://www.grants.gov/search-results-detail/353777 - important
