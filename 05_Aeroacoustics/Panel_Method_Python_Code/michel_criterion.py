
import numpy as np  
def michel_criterion(Re_theta,Re_x): 

    transition_criteria  = Re_theta - 1.174*(1 + 22400/Re_x)*Re_x**0.46  
    try: 
        transition_index     = np.where(transition_criteria>0)[0][0]  
        transition_flag = True 
    except:
        transition_index     = len(Re_theta) -1 
        print('flow does not transition')
        transition_flag = False 
        
    return  transition_index , transition_flag