## @ingroup Library-Methods-Stability-Dynamic_Stability 
# RCAIDE/Library/Methods/Stability/Dynamic_Stability/compute_dynamic_flight_modes.py
# 
# 
# Created:  Apr 2024, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------

# RCAIDE imports
import RCAIDE
from RCAIDE.Framework.Core                            import Units 
from RCAIDE.Library.Components.Wings.Control_Surfaces import Aileron , Elevator  

# python imports 
import numpy as np 


def main():

    b_ref  = 0
    c_ref  = 0
    S_ref  = 0 
    Ixx    = 0
    Iyy    = 0
    Izz    = 0
    Ixy    = 0
    m      = 0
    theta0 =  0
    u0      =  0


    num_cases = 4  # two references, RCAIDE-AVL and RCAIDE-VLM
    ALon = np.zeros((num_cases,4,4))                   

    '''TO BE COMPLETED USING Table 6.1State-Space Modeling of the Rigid-Body Dynamics
   of a Navion Airplane From Flight Data, Using Frequency-Domain Identification Techniquesin and OUR CODE'''

    Cw         = np.array([0,0,0,0 ]) # Should be four values 
    Xu         = np.array([0,0,0,0 ]) # Should be four values 
    Xw         = np.array([0,0,0,0 ]) # Should be four values 
    Xq         = np.array([0,0,0,0 ]) # Should be four values 
    Zu         = np.array([0,0,0,0 ]) # Should be four values 
    Zw         = np.array([0,0,0,0 ]) # Should be four values 
    Zq         = np.array([0,0,0,0 ]) # Should be four values 
    Mu         = np.array([0,0,0,0 ]) # Should be four values 
    Mw         = np.array([0,0,0,0 ]) # Should be four values 
    Mq         = np.array([0,0,0,0 ]) # Should be four values 
    ZwDot      = np.array([0,0,0,0 ]) # Should be four values 
    MwDot      = np.array([0,0,0,0 ]) # Should be four values 


    ALon[:,0,0] = (Xu / m).T[0]
    ALon[:,0,1] = (Xw / m).T[0]
    #ALon[:,0,2] =  Xq.T[0] / m 
    ALon[:,0,3] = (-g * np.cos(theta0)).T[0]
    ALon[:,1,0] = (Zu / (m - ZwDot)).T[0]
    ALon[:,1,1] = (Zw / (m - ZwDot)).T[0]
    ALon[:,1,2] = ((Zq + (m * u0)) / (m - ZwDot) ).T[0]
    ALon[:,1,3] = (-m * g * np.sin(theta0) / (m - ZwDot)).T[0]
    ALon[:,2,0] = ((MwDot * Zu / (m - ZwDot)) / Iyy).T[0]  # ((Mu + MwDot * Zu / (m - ZwDot)) / Iyy).T[0] 
    ALon[:,2,1] = ((Mw + MwDot * Zw / (m - ZwDot)) / Iyy).T[0] 
    ALon[:,2,2] = ((Mq + MwDot * (Zq + m * u0) / (m - ZwDot)) / Iyy ).T[0] 
    ALon[:,2,3] = (-MwDot * m * g * np.sin(theta0) / (Iyy * (m - ZwDot))).T[0] 
    ALon[:,3,0] = 0
    ALon[:,3,1] = 0
    ALon[:,3,2] = 1
    ALon[:,3,3] = 0


    # Look at eigenvalues and eigenvectors
    LonModes                  = np.zeros((num_cases,4), dtype = complex)
    phugoidFreqHz             = np.zeros((num_cases,1))
    phugoidDamping            = np.zeros((num_cases,1))
    phugoidTimeDoubleHalf     = np.zeros((num_cases,1))
    shortPeriodFreqHz         = np.zeros((num_cases,1))
    shortPeriodDamping        = np.zeros((num_cases,1))
    shortPeriodTimeDoubleHalf = np.zeros((num_cases,1))

    for i in range(num_cases):
        D  , V = np.linalg.eig(ALon) # State order: u, w, q, theta
        LonModes[i,:] = D

        # Find phugoid
        phugoidInd               = np.argmax(D)  
        phugoidFreqHz[i]         = abs(D[phugoidInd]) / (2 * np.pi)
        phugoidDamping[i]        = np.sqrt(1/ (1 + ( D[phugoidInd].imag/ D[phugoidInd].real )**2 ))  
        phugoidTimeDoubleHalf[i] = np.log(2) / abs(2 * np.pi * phugoidFreqHz[i] * phugoidDamping[i])

        # Find short period
        shortPeriodInd               = np.argmin(D)  
        shortPeriodFreqHz[i]         = abs(D[shortPeriodInd]) / (2 * np.pi)
        shortPeriodDamping[i]        = np.sqrt(1/ (1 + (D[shortPeriodInd].imag/D[shortPeriodInd].real)**2 ))  
        shortPeriodTimeDoubleHalf[i] = np.log(2) / abs(2 * np.pi * shortPeriodFreqHz[i] * shortPeriodDamping[i]) 
 
        Ixp  = 0 
        Izp  = 0 
        Ixzp = 0 

        for i in range(num_cases): 
            R       = np.array( [[np.cos(AoA[i][0]*Units.degrees) ,  - np.sin(AoA[i][0]*Units.degrees) ], [ np.sin( AoA[i][0]*Units.degrees) , np.cos(AoA[i][0]*Units.degrees)]])
            modI    = np.array([[moments_of_inertia[0][0],moments_of_inertia[0][2]],[moments_of_inertia[2][0],moments_of_inertia[2][2]]] ) 
            INew    = R * modI  * np.transpose(R)
            IxxStab =  INew[0,0]
            IxzStab = -INew[0,1]
            IzzStab =  INew[1,1]
            Ixp[i]  = (IxxStab * IzzStab - IxzStab**2) / IzzStab
            Izp[i]  = (IxxStab * IzzStab - IxzStab**2) / IxxStab
            Ixzp[i] = IxzStab / (IxxStab * IzzStab - IxzStab**2) 

        Yv = np.array([0,0,0,0 ]) # Should be four values 
        Yp = np.array([0,0,0,0 ]) # Should be four values 
        Yr = np.array([0,0,0,0 ]) # Should be four values 
        Lv = np.array([0,0,0,0 ]) # Should be four values 
        Lp = np.array([0,0,0,0 ]) # Should be four values 
        Lr = np.array([0,0,0,0 ]) # Should be four values 
        Nv = np.array([0,0,0,0 ]) # Should be four values 
        Np = np.array([0,0,0,0 ]) # Should be four values 
        Nr = np.array([0,0,0,0 ]) # Should be four values   

        ALat[:,0,0] = (Yv / m).T[0] 
        ALat[:,0,1] = (Yp / m).T[0] 
        ALat[:,0,2] = (Yr/m - u0).T[0] 
        ALat[:,0,3] = (g * np.cos(theta0)).T[0] 
        ALat[:,1,0] = (Lv / Ixp + Ixzp * Nv).T[0] 
        ALat[:,1,1] = (Lp / Ixp + Ixzp * Np).T[0] 
        ALat[:,1,2] = (Lr / Ixp + Ixzp * Nr).T[0] 
        ALat[:,1,3] = 0
        ALat[:,2,0] = (Ixzp * Lv + Nv / Izp).T[0] 
        ALat[:,2,1] = (Ixzp * Lp + Np / Izp).T[0] 
        ALat[:,2,2] = (Ixzp * Lr + Nr / Izp).T[0] 
        ALat[:,2,3] = 0
        ALat[:,3,0] = 0
        ALat[:,3,1] = 1
        ALat[:,3,2] = (np.tan(theta0)).T[0] 
        ALat[:,3,3] = 0

        LatModes                    = np.zeros((num_cases,4),dtype=complex)
        dutchRollFreqHz             = np.zeros((num_cases,1))
        dutchRollDamping            = np.zeros((num_cases,1))
        dutchRollTimeDoubleHalf     = np.zeros((num_cases,1))
        rollSubsistenceFreqHz       = np.zeros((num_cases,1))
        rollSubsistenceTimeConstant = np.zeros((num_cases,1))
        rollSubsistenceDamping      = np.zeros((num_cases,1))
        spiralFreqHz                = np.zeros((num_cases,1))
        spiralTimeDoubleHalf        = np.zeros((num_cases,1))
        spiralDamping               = np.zeros((num_cases,1))
        dutchRoll_mode_real         = np.zeros((num_cases,1))

        for i in range(num_cases):        
            D  , V = np.linalg.eig(ALat[i,:,:]) # State order: u, w, q, theta
            LatModes[i,:] = D  

            # Find dutch roll (complex pair)
            done = 0
            for j in range(3):
                for k in range(j+1,4):
                    if LatModes[i,j].real ==  LatModes[i,k].real:
                        dutchRollFreqHz[i]         = abs(LatModes[i,j]) / (2 * np.pi)
                        dutchRollDamping[i]        = np.sqrt(1/ (1 + ( D[i].imag/ D[i].real )**2 ))  
                        dutchRollTimeDoubleHalf[i] = np.log(2) / abs(2 * np.pi * dutchRollFreqHz[i] * dutchRollDamping[i])
                        dutchRoll_mode_real[i]     = LatModes[i,j].real / (2 * np.pi)
                        done = 1
                        break  
                if done:
                    break  

            # Find roll mode
            diff_vec = np.arange(0,4)
            tmpInd   = np.setdiff1d(diff_vec , [j,k])
            rollInd  = np.argmax(abs(LatModes[i,tmpInd])) # higher frequency than spiral
            rollInd  = tmpInd[rollInd]
            rollSubsistenceFreqHz[i]       = abs(LatModes[i,rollInd]) / 2 / np.pi
            rollSubsistenceDamping[i]      = - np.sign(LatModes[i,rollInd].real)
            rollSubsistenceTimeConstant[i] = 1 / (2 * np.pi * rollSubsistenceFreqHz[i] * rollSubsistenceDamping[i])

            # Find spiral mode
            spiralInd               = np.setdiff1d(diff_vec,[j,k,rollInd])
            spiralFreqHz[i]         = abs(LatModes[i,spiralInd]) / 2 / np.pi
            spiralDamping[i]        = - np.sign(LatModes[i,spiralInd].real)
            spiralTimeDoubleHalf[i] = np.log(2) / abs(2 * np.pi * spiralFreqHz[i] * spiralDamping[i])



        return 


    if __name__ == '__main__': 
        main()    
        plt.show()