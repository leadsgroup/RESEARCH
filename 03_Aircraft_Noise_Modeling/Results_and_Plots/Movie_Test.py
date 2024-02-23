import pylab as plt
import numpy
import matplotlib.animation as animation
#plt.rcParams['animation.ffmpeg_path'] = r"C:\some_path\ffmpeg.exe"   # if necessary
import matplotlib
matplotlib.use("Agg") 
from matplotlib.animation import FFMpegWriter

def main(): 
    # Generate data for plotting
    Lx = Ly = 3
    Nx = Ny = 11
    Nt = 20
    x = numpy.linspace(0, Lx, Nx)
    y = numpy.linspace(0, Ly, Ny)
    x,y = numpy.meshgrid(x,y)
    z0 = numpy.exp(-(x-Lx/2)**2-(y-Ly/2)**2)   # 2 dimensional Gaussian
    
    
    fig = plt.figure()
    fig.set_size_inches(10,5)
    ax = plt.axes(xlim=(0, Lx), ylim=(0, Ly), xlabel='x', ylabel='y')
    
    cvals = numpy.linspace(0,1,Nt+1)      # set contour values 
    cont = plt.contourf(x, y, some_data(0,z0,Nt), cvals)    # first image on screen
    plt.colorbar() 

    
    writer = FFMpegWriter(fps=1) 
    with writer.saving(fig, 'animation_3.mp4', 1000):
        for i in range(Nt):   
            z = some_data(i,z0,Nt)
            for c in cont.collections:
                c.remove()  # removes only the contours, leaves the rest intact
            cont = plt.contourf(x, y, z, cvals)
            plt.title('t = %i:  %.2f' % (i,z[5,5]))    

            writer.grab_frame()      
    
    return 

def some_data(i,z0,Nt):   # function returns a 2D data array
    return z0 * (i/Nt)

# animation function
def animate(i):
    global cont
    z = some_data(i)
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
    cont = plt.contourf(x, y, z, cvals)
    plt.title('t = %i:  %.2f' % (i,z[5,5]))
    return cont


if __name__ == '__main__': 
    main()     
    #plt.show() 