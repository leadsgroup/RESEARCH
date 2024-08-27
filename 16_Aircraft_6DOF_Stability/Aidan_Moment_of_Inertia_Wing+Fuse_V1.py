def compute_wing_moment_of_inertia(self, center_of_gravity):
    '''AIDAN WILL COMPUTE MOMENT OF INERTIA OF WING'''
    
    # note that self a wing
    wing =  self
    
    # ADD CODE
      
    m = mass_of_component
    tr = wing.thickness_to_chord #root thickness
    tt = wing.thickness_to_chord #tip thickness
    ct = wing.chords.tip # tip chord
    cr = wing.chords.root # root chord
    
    # below a0-a4 values are defined for a NACA 4-digit airfoil. This holds for all NACA airfoils
    a0 = 2.969
    a1 = -1.260
    a2 = -3.516
    a3 = 2.843
    a4 = -1.015
    
    
    b = wing.spans.total /2 # half-span
    A = wing.sweeps.quarter_chord * math.pi / 180 # sweep angle in radians (quarter chord)
    
    # Calculate constants
    kd = (tr * (ct + cr) * (2 * cr ** 2 + cr * ct + 2 * ct ** 2)
          + tt * (cr ** 3 + 3 * cr ** 2 * ct + 6 * cr * ct ** 2 + 10 * ct ** 3))
    ke = (tr * (5 * cr ** 4 + 4 * cr ** 3 * ct + 3 * cr ** 2 * ct ** 2 + 2 * cr * ct ** 3 + ct ** 4)
          + tt * (cr ** 4 + 2 * cr ** 3 * ct + 3 * cr ** 2 * ct ** 2 + 4 * cr * ct ** 3 + 5 * ct ** 4))    
    kf = (tr * (cr ** 2 + 2 * cr * ct + 2 * ct ** 2) + tt * (cr ** 2 + 4 * cr * ct
                                                             + 10 * ct ** 2))
    kg = (tr ** 3 * (35 * cr ** 4 + 20 * cr ** 3 * ct + 10 * cr ** 2 * ct ** 2 + 4 * cr * ct ** 3 + ct ** 4)
          + tr ** 2 * tt * (15 * cr ** 4 + 20 * cr ** 3 * ct + 18 * cr ** 2 * ct ** 2 + 12 * cr * ct ** 3 + 5 * ct ** 4)
          + tr * tt ** 2 * (5 * cr ** 4 + 12 * cr ** 3 * ct + 18 * cr ** 2 * ct ** 2 + 20 * cr * ct ** 3 + 15 * ct ** 4)
          + tt ** 3 * (cr ** 4 + 4 * cr ** 3 * ct + 10 * cr ** 2 * ct ** 2 + 20 * cr * ct ** 3 + 35 * ct ** 4))
    vo = 1 / 60 * (40 * a0 + 30 * a1 + 20 * a2 + 15 * a3 + 12 * a4) # NACA 4 digit integral of thickness distribution.
    v1 = 1 / 60 * (56 * a0 + 50 * a1 + 40 * a2 + 33 * a3 + 28 * a4)
    v2 = 1 / 980 * (856 * a0 + 770 * a1 + 644 * a2 + 553 * a3 + 484 * a4)
    v3 = (2 / 5 * a0 ** 3 + a0 ** 2 * a1 + 3 / 4 * a0 ** 2 * a2 + 3 / 5 * a0 ** 2 * a3 + 1 / 2 * a0 ** 2 * a4 + 6 / 7 * a0 * a1 ** 2
          + 4 / 3 * a0 * a1 * a2 + 12 / 11 * a0 * a1 * a3 + 12 / 13 * a0 * a1 * a4 + 6 / 11 * a0 * a2 ** 2 + 12 / 13 * a0 * a2 * a3
          + 4 / 5 * a0 * a2 * a4 + 2 / 5 * a0 * a3 ** 2 + 12 / 17 * a0 * a3 * a4 + 6 / 19 * a0 * a4 ** 2 + 1 / 4 * a1 ** 3
          + 3 / 5 * a1 ** 2 * a2 + 1 / 2 * a1 ** 2 * a3 + 3 / 7 * a1 ** 2 * a4 + 1 / 2 * a1 * a2 ** 2 + 6 / 7 * a1 * a2 * a3
          + 3 / 4 * a1 * a2 * a4 + 3 / 8 * a1 * a3 ** 2 + 2 / 3 * a1 * a3 * a4 + 3 / 10 * a1 * a4 ** 2 + 1 / 7 * a2 ** 3
          + 3 / 8 * a2 ** 2 * a3 + 1 / 3 * a2 ** 2 * a4 + 1 / 3 * a2 * a3 ** 2 + 3 / 5 * a2 * a3 * a4 + 3 / 11 * a2 * a4 ** 2
          + 1 / 10 * a3 ** 3 + 3 / 11 * a3 ** 2 * a4 + 1 / 4 * a3 * a4 ** 2 + 1 / 13 * a4 ** 3)
    
    delta = 1 # 1 for right wing, -1 for left wing
    
    # Moment of inertia in local system
    Ixx = m * (56 * b ** 2 * kf * vo + kg * v3) / (280 * ka * vo)
    Iyy = m * (84 * b * (2 * b * kf * vo * math.tan(A) ** 2 + kd * v1 * math.tan(A)) + 49 * ke * v2 + 3 * kg * v3) / (840 * ka * vo)
    Izz = m * (12 * b * (2 * b * (math.tan(A) ** 2 + 1) * kf * vo + kd * v1 * math.tan(A)) + 7 * ke * v2) / (120 * ka * vo)
    Ixy = -1 * delta * b * m * (4 * b * kf * vo * math.tan(A) + kd * v1) / (20 * ka * vo)
    ## Ixz, Iyz are 0
    Ixz = 0
    Iyz = 0
    I = [[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]]
    
    # global system
    s = wing.origin - np.array(center_of_gravity) # Vector for the parallel axis theorem
    
    I_global = np.array(I) + m * np.array(np.dot(s, s)) * np.array(np.identity(3)) - np.array(np.dot(s, np.transpose(s)))
    
    return I_global


## By segments:

def compute_fuselage_moment_of_inertia(self, center_of_gravity):
    '''AIDAN WILL COMPUTE MOMENT OF INERTIA OF FUSELAGE'''
    import numpy as np
    # ADD CODE
    I_total = np.zeros((3, 3))
    
    # note that self a fuselage
    fuselage =  self
    
    for segment in fuselage.Segments:
        outer_radius = segment.height # segment.width # Figure out ellipse to circle conversion or moment of inertia of an elliptical rectangular solid
        #inner_radius = segment. STUFF
        h = fuselage.lengths * (segment.)
        m = mass_of_component
        
        # vector for parallel axis theorem
        s =  self.origin[0] - center_of_gravity # check type and order/vector
        
        I =  np.zeros((3, 3))
        
        # Moment of inertia in local system
        I[0][0] =  m / 2 *  (outer_radius ** 2 + inner_radius ** 2) # Ixx
        I[1][1] =  m / 12 * (3 * (outer_radius ** 2 + inner_radius ** 2) + h ** 2) # Iyy
        I[2][2] =  m / 12 * (3 * (outer_radius ** 2 + inner_radius ** 2) + h ** 2) # Izz
        
        # transform moment of inertia to global system
        I_global = np.array(I) + m * np.array(np.dot(s, s)) * np.array(np.identity(3)) - np.array(np.dot(s, np.transpose(s)))
        I_total = np.array(I_total) + np.array(I_global)
    
    return I_total

## By hemisphere, cylnder, and cone
## as of 8/10 needs mass and then to be tested. Unsure about numpy array addition.
## will need to add math and numpy library to Fuselage class

def compute_fuselage_moment_of_inertia(self, center_of_gravity):
    '''AIDAN WILL COMPUTE MOMENT OF INERTIA OF FUSELAGE'''
    import numpy as np
    # ADD CODE
    I_total = np.zeros((3, 3))
    
    # note that self a fuselage
    fuselage =  self
    
    ## Hemisphere
    
    origin_hemisphere =  np.array([fuselage.lengths.nose, 0, 0]) + np.array(fuselage.origin)
    m = mass_of_component
    I =  np.zeros((3, 3))
    outer_radius = fuselage.effective_diameter / 2
    inner_radius = 0.75 * fuselage.effective_diameter / 2 # Assume the inner radius is 75 % of the outer radius
    
    # Moment of inertia in local system
    I[0][0] =  2 * m / 5 *  (outer_radius ** 5 - inner_radius ** 5) /(outer_radius **3 -inner_radius **3) # Ixx
    I[1][1] =  2 * m / 5 *  (outer_radius ** 5 - inner_radius ** 5) /(outer_radius **3 -inner_radius **3)# Iyy
    I[2][2] =  2 * m / 5 *  (outer_radius ** 5 - inner_radius ** 5) /(outer_radius **3 -inner_radius **3) # Izz

    # global system
    s = np.array(origin_hemisphere) - np.array(center_of_gravity)
    
    I_global = np.array(I) + m * np.array(np.dot(s, s)) * np.array(np.identity(3)) - np.array(np.dot(s, np.transpose(s)))
    I_total = np.array(I_total) + np.array(I_global)
    
    ## cylinder
    
    h = fuselage.lengths.total - fuselage.lengths.nose - fuselage.lengths.tail
    m = mass_of_component
    origin_cylinder =  np.array([fuselage.lengths.nose + h / 2,0, 0]) + np.array(fuselage.origin) # origin of the cylinder is located a tthe middle of the cylinder
    
    I =  np.zeros((3, 3))
    
    # Moment of inertia in local system
    I[0][0] =  m / 2 *  (outer_radius ** 2 + inner_radius ** 2) # Ixx
    I[1][1] =  m / 12 * (3 * (outer_radius ** 2 + inner_radius ** 2) + h ** 2) # Iyy
    I[2][2] =  m / 12 * (3 * (outer_radius ** 2 + inner_radius ** 2) + h ** 2) # Izz
    
    # transform moment of inertia to global system
    s = np.array(origin_cylinder) - np.array(center_of_gravity)
    
    I_global = np.array(I) + m * np.array(np.dot(s, s)) * np.array(np.identity(3)) - np.array(np.dot(s, np.transpose(s)))
    I_total = np.array(I_total) + np.array(I_global)

    ## cone
    
    h = fuselage.lengths.tail
    m = mass_of_component
    origin_cone =  np.array([fuselage.lengths.total - h,0, 0]) + np.array(self.origin)
    
    I =  np.zeros((3, 3))
    
    # Moment of inertia in local system
    rho = (m / (1 / 3 * math.pi * (outer_radius ** 2 * h - inner_radius ** 2 * (h * inner_radius / outer_radius))))
    I[0][0] =  rho * (1 / 3 * math.pi * outer_radius ** 2 * h ** 3 + math.pi / 20 * outer_radius ** 4 *h - 1 / 3 * math.pi * inner_radius ** 2 * (h * inner_radius / outer_radius) ** 3 + math.pi / 20 * inner_radius ** 4 *(h * inner_radius / outer_radius))
    I[1][1] =  rho * (1 / 3 * math.pi * outer_radius ** 2 * h ** 3 + math.pi / 20 * outer_radius ** 4 *h - 1 / 3 * math.pi * inner_radius ** 2 * (h * inner_radius / outer_radius) ** 3 + math.pi / 20 * inner_radius ** 4 *(h * inner_radius / outer_radius))
    I[2][2] =  rho * (math.pi /10 *outer_radius **4 *h -math.pi /10 *inner_radius **4 *(h * inner_radius / outer_radius)) # Izz
    
    # transform moment of inertia to global system
    
    s = np.array(origin_cone) - np.array(center_of_gravity)
    
    I_global = np.array(I) + m * np.array(np.dot(s, s)) * np.array(np.identity(3)) - np.array(np.dot(s, np.transpose(s)))
    I_total = np.array(I_total) + np.array(I_global)
    
    
    return I_total