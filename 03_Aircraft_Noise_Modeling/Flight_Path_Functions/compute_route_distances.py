from  RCAIDE.Framework.Core import  Units 
import  numpy as  np

def compute_route_distances(x1, y1, x2, y2, radius_Vert1, radius_Vert2, dep_heading_rad, app_heading_rad,high_speed_climb_distance):
    
    dep_heading  =  dep_heading_rad / Units.degree
    app_heading  =  app_heading_rad / Units.degree
    dep_distance = radius_Vert1
    app_distance = radius_Vert2
    
    # ---------------------------------------------------------------------------------------------------------------------- 
    # Path heading:
    # ---------------------------------------------------------------------------------------------------------------------- 
    path_heading = (np.arctan2((x2 - x1), (y2 - y1)) *180 /np.pi) % 360
    
    # ---------------------------------------------------------------------------------------------------------------------- 
    # Departure sector: follows a right-hand pattern
    # ---------------------------------------------------------------------------------------------------------------------- 
    if (path_heading > dep_heading):
        dep_sector = path_heading - dep_heading # negative means counter clockwise, positive means clockwise
    else:
        dep_sector = 360 + path_heading - dep_heading
    
    # ---------------------------------------------------------------------------------------------------------------------- 
    # Departure pattern distance:
    # ---------------------------------------------------------------------------------------------------------------------- 
    dep_pattern_distance = abs(dep_sector*np.pi/180) *radius_Vert1
    
    # ---------------------------------------------------------------------------------------------------------------------- 
    # Path distance:
    # ---------------------------------------------------------------------------------------------------------------------- 
    path_distance = np.sqrt((y1-y2)**2 + (x1-x2)**2) - radius_Vert1 - radius_Vert2 -  high_speed_climb_distance  
    
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
    p2p_distance = np.sqrt((y1-y2)**2 + (x1-x2)**2)
    
    return path_distance, path_heading, dep_sector, app_sector