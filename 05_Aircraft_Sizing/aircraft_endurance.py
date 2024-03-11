# imports  
import numpy as np         
import matplotlib.pyplot as plt  

def main():  
    
    # symmetry 
    n_pts           = 20
    span            = 1  
    semi_points     = int(np.ceil(n_pts/2))
    n               = np.linspace(semi_points,0,semi_points)    
    thetan          = n*(np.pi/2)/(semi_points+1)               
    x_points_half   = (span/2)*np.cos(thetan)                   
    x_points        = np.zeros(n_pts+1)
    x_points[0]     = 0 
    x_points[-1]    = span  
    x_points[1:-1]  = np.concatenate((x_points_half[:-1], (span - x_points_half[::-1])))
    
    

    # symmetry 
    n_pts           = 20
    span            = 1  
                
    semi_points     = int(np.ceil(n_pts/2))
    n               = np.linspace(n_pts,0,n_pts)    
    thetan          = n*(np.pi/2)/(n_pts)     
    x2_points       = span-(span)*np.cos(thetan)  
    
    
    y  = np.ones(len(x_points))
    
    fig = plt.figure()
    axis = fig.add_subplot(2,1,1) 
    axis.plot(x_points, y, marker = 'o')

    y2  = np.ones(len(x2_points))*2
        

    axis = fig.add_subplot(2,1,2) 
    axis.plot(x2_points, y2, marker = 'o')    
    return 
     
if __name__ == '__main__':
    main()
    plt.show()