
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from RCAIDE.Framework.Core import Units, Data  
from RCAIDE.Library.Plots        import *  

import pickle
import time  
import numpy as np
import pylab as plt
import pandas as pd
import sys
import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

sys.path.append('Common')  
import Vehicle
import Configurations
import Analyses 
import Missions
import Plots  
 

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():  
    # start simulation clock
    ti                         = time.time()
    RUN_NEW_MODEL_FLAG         = True
    
    # -------------------------------------------------------------------------------------------    
    # SET UP SIMULATION PARAMETERS   
    # -------------------------------------------------------------------------------------------  
    simulated_days             = 1               # number of days simulated
    flights_per_day            = 1               # number of flights per day  
    recharge_battery           = True            # flag to simulate battery recharge  
    plot_mission               = True            # plot mission flag  
    resize_aircraft            = False
   
    if RUN_NEW_MODEL_FLAG:    
         
        # -------------------------------------------------------------------------------------------    
        # SET UP VEHICLE
        # -------------------------------------------------------------------------------------------  
        vehicle = Vehicle.vehicle_setup(resize_aircraft,'e_twin_otter_vehicle')
        configs =  Configurations(vehicle)
    
       
        # -------------------------------------------------------------------------------------------    
        # SET UP MISSION PROFILE  
        # -------------------------------------------------------------------------------------------    
        for day in range(simulated_days):          
            
            print(' *********** '+ str(day) +' ***********  ')        
            base_mission[day]      = Missions.repeated_flight_operation_setup(vehicle,day,flights_per_day, recharge_battery)
            missions_analyses = Missions.missions_setup(base_mission[day])
       
    
        # -------------------------------------------------------------------------------------------    
        # RUN SIMULATION !!
        # -------------------------------------------------------------------------------------------
        results     = missions_analyses.base.evaluate() 
        
        # -------------------------------------------------------------------------------------------    
        # SAVE RESULTS
        # -------------------------------------------------------------------------------------------
        filename          = 'e_Twin_Otter_'
        save_results(results,filename)   
    
    else:
        filename          = 'e_Twin_Otter_' 
        results = load_results(filename) 
        
    if plot_mission: 
        Plots.plot_results(results,save_figure_flag = False)       
    
    
    tf = time.time() 
    print ('time taken: '+ str(round(((tf-ti)/60),3)) + ' mins')
    
    
    for i in range(len(flight_no)):
        starting_range =  results.segments[1+(i*14)].conditions.frames.inertial.aircraft_range[-1,0] /Units.nmi
        final_range = results.segments[12+(i*14)].conditions.frames.inertial.aircraft_range[-1,0]/Units.nmi
        print(final_range-starting_range)
    
    return 


# ----------------------------------------------------------------------
#   Save Results
# ----------------------------------------------------------------------
def save_results(results,filename): 
    pickle_file  =  filename + '.pkl'
    with open(pickle_file, 'wb') as file:
        pickle.dump(results, file) 
    return   

# ------------------------------------------------------------------
#   Load Results
# ------------------------------------------------------------------   
def load_results(filename):  
    load_file = filename + '.pkl' 
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results

  
    

def show_notification():
    os.system('osascript -e \'display notification "The simulation has completed successfully." with title "Simulation Complete"\'')

def play_sound():
    # Play a default system sound on macOS
    os.system('afplay /System/Library/Sounds/Glass.aiff')

def simulation_complete(airport, month):
    play_sound()
    show_notification()
    send_email()
    
def send_email():
    # Outlook credentials
    yahoo_user = 'sai.sankalp@yahoo.com'
    yahoo_password = 'zlofcfakwwmutfwl'
    
    # Email content
    msg = MIMEMultipart()
    msg['From'] = yahoo_user
    msg['To'] = '4085921097@tmomail.net'  
    body = 'The simulation is complete.'
    msg.attach(MIMEText(body, 'plain'))

    try:
        # Setup the server and send email
        server = smtplib.SMTP('smtp.mail.yahoo.com', 587)
        server.starttls()
        server.login(yahoo_user, yahoo_password)
        text = msg.as_string()
        server.sendmail(yahoo_user, msg['To'], text)
        server.quit()
        print('Email sent successfully!')
    except Exception as e:
        print(f'Failed to send email: {e}')

    

if __name__ == '__main__':
    main()
    simulation_complete
    plt.show()