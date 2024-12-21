import  numpy as  np
import  scipy.integrate as  integrate


def Vh_linear(t, v_climb, Vf, Vi, height):
    t_total = height / v_climb
    a = (Vf - Vi) / t_total
    Vh = np.sqrt(a**2*t**2+2*a*t*Vi+Vi**2-v_climb**2)
    
    return Vh
def Linear_Accel_Climb_Distance(v_start,  v_end, v_climb, height):    
    t_final = height / v_climb
    
    dist = integrate.quad(lambda t: Vh_linear(t, v_climb, v_end, v_start, height), 0, t_final)[0]
    return(dist)