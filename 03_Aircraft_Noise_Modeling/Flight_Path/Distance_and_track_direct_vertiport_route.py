# ---------------------------------------------------------------------------------------------------------------------- 
# This script calculates the distance of a flight from one vertiport to another considering a circular pattern around the airport and flying point-to-point.
#
# Y is local approximation for North, X is local approximation for East
# In following convention a heading of 0 degrees is due north or along the Y-axis and goes clockwise
# Assumes a right hand pattern
#
# Departure pattern: circuit taken after takeoff around the first vertiport
# Path: route taken between the departure pattern and approach pattern between the two vertiports
# Approach pattern: circuit taken before the landing sequence once the eVTOL arrives near vertiport 2
#
# Last modified: 9/4/2024
# ---------------------------------------------------------------------------------------------------------------------- 

#import RCAIDE
#from RCAIDE.Framework.Core import Units   
import numpy as np
from matplotlib import pyplot as plt

# ---------------------------------------------------------------------------------------------------------------------- 
# Inputs
# ---------------------------------------------------------------------------------------------------------------------- 

vert_Loc1 = np.array([0, 0]) #[0 * Units.km, 0 * Units.km]) # X,Y coordinates of vertiport 1
vert_Loc2 = np.array([-10, -10]) #[10 * Units.km, 4 * Units.km]) # X,Y coordinates of vertiport 2
radius_Vert1 = 1 #* Units.km # circular pattern radius around vertiport 1
radius_Vert2 = 1 #* Units.km # circular pattern radius around vertiport 2
dep_Heading = 200 # Heading [degrees] of the departure from vertiport 1
app_heading = 90 # Heading [degrees] of the approach to vertiport 2

dep_distance = radius_Vert1
app_distance = radius_Vert2

# ---------------------------------------------------------------------------------------------------------------------- 
# Path heading:
# ---------------------------------------------------------------------------------------------------------------------- 
path_Slope = (vert_Loc2[1] -vert_Loc1[1])/(vert_Loc2[0] -vert_Loc1[0])
if (vert_Loc2[1] - vert_Loc1[1]) < 0:
    path_heading = (np.arctan(1/path_Slope) *180 / np.pi + 180) % 360
else:
    path_heading = (np.arctan(1/path_Slope) *180 / np.pi) % 360

# ---------------------------------------------------------------------------------------------------------------------- 
# Departure sector: follows a right-hand pattern
# ---------------------------------------------------------------------------------------------------------------------- 
if (path_heading > dep_Heading):
    dep_sector = path_heading - dep_Heading # negative means counter clockwise, positive means clockwise
else:
    dep_sector = 360 + path_heading - dep_Heading

# ---------------------------------------------------------------------------------------------------------------------- 
# Departure pattern distance:
# ---------------------------------------------------------------------------------------------------------------------- 
dep_pattern_distance = abs(dep_sector*np.pi/180) *radius_Vert1

# ---------------------------------------------------------------------------------------------------------------------- 
# Path distance:
# ---------------------------------------------------------------------------------------------------------------------- 
path_distance = np.sqrt((vert_Loc1[1]-vert_Loc2[1])**2 + (vert_Loc1[0]-vert_Loc2[0])**2) - radius_Vert1 - radius_Vert2

# ---------------------------------------------------------------------------------------------------------------------- 
# Path approach heading: The angle measures on the circle that is the traffic apttern are 180 degrees offset from the heading 
# ---------------------------------------------------------------------------------------------------------------------- 
path_heading_angle = (path_heading + 180) % 360 # in
app_angle = (app_heading + 180) % 360 #out 

# ---------------------------------------------------------------------------------------------------------------------- 
# Approach sector: assumes a right hand pattern 
# ---------------------------------------------------------------------------------------------------------------------- 
if path_heading_angle > app_angle:
    app_sector = 360 - (path_heading_angle - app_angle)
else :
    app_sector = app_angle - path_heading_angle

# ---------------------------------------------------------------------------------------------------------------------- 
# Approach pattern distance:
# ---------------------------------------------------------------------------------------------------------------------- 
app_pattern_distance = abs(app_sector*np.pi/180) * radius_Vert2

# ---------------------------------------------------------------------------------------------------------------------- 
# Total distance:
# ---------------------------------------------------------------------------------------------------------------------- 
total_distance = dep_distance + dep_pattern_distance + path_distance + app_pattern_distance + app_distance

# ---------------------------------------------------------------------------------------------------------------------- 
# Point-to-point distance:
# ---------------------------------------------------------------------------------------------------------------------- 
p2p_distance = np.sqrt((vert_Loc1[1]-vert_Loc2[1])**2 + (vert_Loc1[0]-vert_Loc2[0])**2)

# ---------------------------------------------------------------------------------------------------------------------- 
# Print distances
# ---------------------------------------------------------------------------------------------------------------------- 
print("Point-to-point Distance: "+str(p2p_distance))
print("Total Distance: "+str(total_distance))
print("Departure Distance: "+str(dep_distance))
print("Departure Pattern Distance: "+str(dep_pattern_distance))
print("Path Distance: "+str(path_distance))
print("Approach Distance: "+str(app_distance))
print("Approach Pattern Distance: "+str(app_pattern_distance))

# ---------------------------------------------------------------------------------------------------------------------- 
# Graphs:
# ---------------------------------------------------------------------------------------------------------------------- 
x_traj = [vert_Loc1[0]] # Array to hold the x coordinates of the trajectory
y_traj = [vert_Loc1[1]] # Array to hold the y coordinates of the trajectory
n = 300 # number of points to plot each pattern with

if path_heading < dep_Heading: # To correctly plot the pattern this ensures a clockwise pattern
    path_heading = path_heading + 360

# segment from departure end to the path 
for i in range(0, n):
    theta = (path_heading - dep_Heading) / n * i + dep_Heading
    x_traj.append(np.sin(np.pi/180*theta)+vert_Loc1[0])
    y_traj.append(np.cos(np.pi/180*theta)+vert_Loc1[1])

if app_angle < path_heading_angle: # THis ensures a right hand pattern
    app_angle = app_angle + 360
for i in range(0, n): # going from path to the approach end
    theta = (app_angle - path_heading_angle) / n * i + path_heading_angle
    x_traj.append(np.sin(np.pi/180*theta)+vert_Loc2[0])
    y_traj.append(np.cos(np.pi/180*theta)+vert_Loc2[1])
        
x_traj.append(vert_Loc2[0])
y_traj.append(vert_Loc2[1])

# Plot circles around vertiports
x_vert1 = []
y_vert1 = []
for i in range(100):
    theta = i * 2 * np.pi / 100
    x_vert1.append(vert_Loc1[0]+radius_Vert1*np.sin(theta))
    y_vert1.append(vert_Loc1[1]+radius_Vert1*np.cos(theta))

x_vert2 = []
y_vert2 = []
for i in range(100):
    theta = i * 2 * np.pi / 100
    x_vert2.append(vert_Loc2[0]+radius_Vert2*np.sin(theta))
    y_vert2.append(vert_Loc2[1]+radius_Vert2*np.cos(theta))

# Show plots   

plt.plot(x_vert1, y_vert1, linestyle='--')
plt.plot(x_vert2, y_vert2, linestyle='--')
plt.plot([vert_Loc1[0], vert_Loc2[0]], [vert_Loc1[1], vert_Loc2[1]], linestyle='--')
plt.text(vert_Loc1[0], vert_Loc1[1], "Vertiport_1")
plt.text(vert_Loc2[0], vert_Loc2[1], "Vertiport_2")
plt.plot(x_traj, y_traj)
plt.title("Direct Route With Traffic Patterns From Vertiport 1 to 2")
plt.show()


