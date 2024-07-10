
# python imports 
import numpy as np
import  scipy
import matplotlib.pyplot as plt  
import os   
from scipy.interpolate import CubicSpline

def main():
    
    x_new   = np.array([0,0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1,3, 5, 7, 8, 9, 10, 11.5, 55, 57, 58, 62, 73]) 
    
    fuselage_width =  np.array([
        [-0.0006225766209242423, 0],
        [0.1400797397080149, 0.409775784753364],
        [0.5584512289692842, 0.819551569506728],
        [1.2544918911628864, 1.229327354260085],
        [2.228201726288817, 1.63910313901345],
        [3.8967070703664985, 2.185470852017936],
        [5.704047000910347, 2.7318385650224215],
        [7.232472605280017, 3.0050224215246644],
        [8.343771873630264, 3.141614349775786],  
        [56.102869617991175, 3.141614349775786],
        [58.46243501129506, 3.0050224215246644],
        [61.37671617384268, 2.7318385650224215],
        [64.15091759668229, 2.185470852017936],
        [69.42165126942918, 1.0927354260089643],
        [73.16582706766918, 0.1365919282511214]])
    
    fit_1 = scipy.interpolate.interp1d(fuselage_width[:,0], fuselage_width[:,1], kind='linear')
    fuselage_width_new =  fit_1(x_new)
    
    
    fuselage_height_top =  np.array([
    [-0.080312384099261, -47.943766816143],
    [-0.07844465423648828, -47.5339910313],
    [0.3405494116457062, -46.987623318385],
    [0.7576757476651252, -46.851031390134],
    [1.7332533126538303, -46.031479820627],
    [2.985254897333019, -45.4851121076233],
    [4.376091068478371, -44.9387443946188],
    [5.765682086381875, -44.6655605381165],
    [7.2941076907515425, -44.392376681614],
    [8.821910718500284, -44.2557847533632],
    [11.46039043797835, -44.1191928251121],
    [56.7204656259483, -44.1191928251121], 
    [59.357700192184495, -44.255784753363], 
    [70.73777824606361, -45.2119282511210],
    [73.51446997538692, -45.2119282511210]])
    
    fuselage_height_top[:,1] =  fuselage_height_top[:,1] +  47.943766816143
    fit_2 = scipy.interpolate.interp1d(fuselage_height_top[:,0], fuselage_height_top[:,1], kind='linear')  
    fuselage_height_top_new =  fit_2(x_new) 
    
    
    fuselage_height_bottom =  np.array([     
    [73.51446997538692, -45.21192825112107],
    [58.63800161839576, -49.85605381165919],
    [57.24841060049228, -50.12923766816143],
    [54.74876546748037, -50.26582959641255],
    [7.545006068984121, -50.26582959641255],
    [5.740779021544895, -50.12923766816143],
    [4.214843723658924, -49.85605381165919],
    [2.827743012239117, -49.58286995515694],
    [1.7189340503725665, -49.1730941704035],
    [0.8871716848174254, -48.8999103139013],
    [0.3330784921946117, -48.6267264573990]])
    

    fuselage_height_bottom[:,1] =  fuselage_height_bottom[:,1] +  47.943766816143
    fuselage_height_bottom[:,1][0] = 2.21
    fit_3 = scipy.interpolate.interp1d(fuselage_height_bottom[:,0][::-1], fuselage_height_bottom[:,1][::-1], kind='linear', fill_value="extrapolate")   
    fuselage_height_bottom_new = fit_3(x_new)
    
    main_wing = np.array([   
    [24.98960299,0], 
    [30.5429,8],
    [32.95,11.44],
    [46.09183755,30.32340807],
    [47.89668718,30.32340807],
    [40.75,11.44],
    [40.1591,8],
    [40.12319549,0.136591928]])  
     
    

    fig = plt.figure()
    fig.set_size_inches(10, 10) 

    axes = plt.subplot(2,2,1)
    axes.plot( fuselage_width[:,0], fuselage_width[:,1], 'ro-') 
    axes.plot( x_new  , fuselage_width_new , 'bo-')
    axes.grid(True)

    axes = plt.subplot(2,2,2)
    axes.plot(fuselage_height_bottom[:,0], fuselage_height_bottom[:,1], 'ro-' ) 
    axes.plot(fuselage_height_top[:,0], fuselage_height_top[:,1], 'ro-' ) 
    axes.plot( x_new  , fuselage_height_top_new, 'bo-' ) 
    axes.grid(True)

    axes = plt.subplot(2,2,3)
    axes.plot(fuselage_height_bottom[:,0], fuselage_height_bottom[:,1], 'ro-' ) 
    axes.plot( x_new , fuselage_height_bottom_new , 'bo-' )  
    axes.grid(True)

    axes = plt.subplot(2,2,4) 
    axes.plot(main_wing[:,0] ,main_wing[:,1], 'bo-' )      
    axes.grid(True)
    
    f_width =  fuselage_width_new * 2
    f_height =  fuselage_height_top_new - fuselage_height_bottom_new
    f_center_z =  (fuselage_height_top_new + fuselage_height_bottom_new) / 2
    
    
    for i in  range(len(f_width)):
        print('\n')
        print('Segment ' + str(i+1) )
        print(x_new[i]/x_new[-1]) 
        print(f_center_z[i]/x_new[-1] ) 
        print(f_height[i]) 
        print(f_width[i]) 
    return

if __name__ == '__main__': 
    main()
    plt.show()
       
