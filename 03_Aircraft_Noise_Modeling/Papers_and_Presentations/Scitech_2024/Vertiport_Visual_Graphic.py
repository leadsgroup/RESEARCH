import  numpy as  np
from  lxml import  etree
from  pykml.parser import  Schema
from pykml.factory import  KML_ElementMaker as  KML

def app_circle(path_heading, app_sector, path_heading):
    app_path_angle = path_heading + np.pi + 0.01 * range(100) *app_sector
    x = radius_vert * np.cos(app_path_angle)
    y = radius_vert * np.sin(app_path_angle)
    
    