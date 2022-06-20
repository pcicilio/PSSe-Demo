from __future__ import division
from collections import defaultdict
import os,sys
import io

#Change to your PSS/e Location and set up paths
sys.path.append(r"C:\Program Files (x86)\PTI\PSSE34\PSSPY37") #Give the path to PSSBIN to imoport psspy
sys.path.append(r"C:\Program Files (x86)\PTI\PSSE34\PSSBIN")
sys.path.append(r"C:\Program Files (x86)\PTI\PSSE34\PSSLIB")
sys.path.append(r"C:\Program Files (x86)\PTI\PSSE34\EXAMPLE")
os.environ['PATH'] = (r"C:\Program Files (x86)\PTI\PSSE34\PSSPY37;" + r"C:\Program Files (x86)\PTI\PSSE34\PSSBIN;" + r"C:\Program Files (x86)\PTI\PSSE34\EXAMPLE;" + os.environ['PATH'])


import psse34 #addition necessary for new version 34
import psspy
#import pssarrays #loads and saves data from PSSe functions (short circuit)
import redirect
import dyntools
import pssplot
#import pssexcel #loads and saves data from PSSe functions (PV, QV, ACCC)
import random
import copy
import math
import multiprocessing
import time
import sys
import csv
import numpy as np
import pandas as pd
import re
import redirect
import silence #This prevents the output bar of PSSe from outputting in the console, which reduces simulation time
from silence import silence
import matplotlib.pyplot as plt
redirect.psse2py();
#redirect.psse2py() 
#import matplotlib
#_i=psspy.getdefaultint()
#_f=psspy.getdefaultreal()
#_s=psspy.getdefaultchar()
#redirect.psse2py()

# For DLL import Error: Following checks path, path is good
import os
PSSE_PATH = r"C:\Program Files (x86)\PTI\PSSE34\PSSPY37" 
if not os.path.exists(PSSE_PATH):    
    print ("Bad path")
#Tries to find dll in path, is also good
import glob
fns = glob.glob(os.path.join(PSSE_PATH, 'psspy.*'))
if len(fns) == 0:
   print ("Who moved my libraries")


## Run PSSe Dyanmic Simulation using PSSe Python API =====================================================================================================
def Run_SIM(sav_file,snp_file,dyr_file,out_file,disturbance_type,channel_option,runtime): #inputs are strings
# This function takes the file names of the power system model and desired PSSe
# output file name. The location of the files will need to be updated if this
# code is being run for the first time. 
# Then the program runs a dynamic simulation, with a branch trip and branch close.
# The channels are also selected (which data is recorded). 
# If a different disturbance, or different simulation runtime, or different data
# channels are desired this code will need to be updated.
# The commands in this code are from the PSSe python API, and there is further
# instruction on what these commands are and their error codes in the API
# documentation, which is located in the PSSe DOCS folder where PSSe is installed.

    # Write .sav case name and .outx file name
    sav = r"""C:\Users\lebredeson\Documents\PSSe_Demo\Demo_Models\%s""" %sav_file    #find .sav file
    out = r"""C:\Users\lebredeson\Documents\PSSe_Demo\Outputs\%s""" %out_file       #name .outx file

    ierr = [1]*30 #check and record for error codes
    output = io.StringIO()
    with silence(output):    
        ierr[0] = psspy.psseinit(200000) #initialize PSSe. This number needs to be high, otherwise there will be not enough output channels available for recording outputs
        ierr[1] = psspy.case(sav) #load case information (.sav file)
        if snp_file != "None":
            snp = r"""C:\Users\lebredeson\Documents\PSSe_Demo\Demo_Models\%s""" %snp_file   #find .snp file
            ierr[2] = psspy.rstr(snp) #load dynamic snapshot information (.snp file)
        ierr[3] = psspy.fnsl([0,0,0,1,1,0,99,0]) #Solves power flow using fixed slope decoupled Newton-Raphson
        #ierr[3] = psspy.mslv([1,0,0,1,1,0]) #solves power flow using modified gauss-seidel method
        ierr[4] = psspy.cong(0)
        ierr[5] = psspy.conl(0,1,1,[0,0],[ 50.0, 50.0, 50.0, 50.0]) #initialize for load conversion
        ierr[6] = psspy.conl(0,1,2,[0,0],[ 50.0, 50.0, 50.0, 50.0]) # convert loads
        ierr[7] = psspy.conl(0,1,3,[0,0],[ 50.0, 50.0, 50.0, 50.0]) # postprocessing housekeeping
        ierr[8] = psspy.ordr(1)
        ierr[9] = psspy.fact()
        ierr[10] = psspy.tysl(0)
        if dyr_file != "None":
            dyre = r"""C:\Users\lebredeson\Documents\PSSe_Demo\Demo_Models\%s""" %dyr_file   #find .dyr file
            ierr[11] = psspy.dyre_new([1,1,1,1],dyre,"","","")
        
        # For when not using a subsystem for channels use chsb(0,1..., with subsystem do chsb(0,0...
        ierr[12] = psspy.delete_all_plot_channels() #Delete channels stored from any previous runs that might be in memory
        # Pick Channels: These are channels the user needs to define, including if 
        #                all buses/machines should be recorded or just a subset (aka subsystem)
        #                and which outputs (frequency, electric or mechanical power, line flows,etc.)
        #                Here are some examples. 
        if channel_option == 'All':
            ierr[12] = psspy.chsb(0,1,[-1,-1,-1,1,2,0]) #Machine electrical power
            ierr[13] = psspy.chsb(0,1,[-1,-1,-1,1,12,0]) #Bus Frequency Deviations (pu)
            ierr[14] = psspy.chsb(0,1,[-1,-1,-1,1,13,0]) #Bus Voltage and angle (pu)    
        elif channel_option == 'GVEA':
            ierr[12] = psspy.bsys(1,0,[0.0,0.0],0,[],9,[325,33601,33500,30110,30120,31001,32001,32300,32100],0,[],0,[])  
            ierr[13] = psspy.chsb(1,0,[-1,-1,-1,1,2,0]) #Machine electrical power
            ierr[14] = psspy.chsb(1,0,[-1,-1,-1,1,12,0]) #Bus Frequency Deviations (pu)
            ierr[15] = psspy.chsb(1,0,[-1,-1,-1,1,13,0]) #Bus Voltage and angle (pu)
        ierr[21] = psspy.strt_2([0,1],out) #Initialize dynamic simulation
        ierr[22] = psspy.run(0, 1.0,1,1,1) #Run simulation for 1 second to verify steady state
        if disturbance_type == "line_fault":
            ierr[23] = psspy.dist_branch_fault(154,3008,r"""1""",1, 230.0,[0.0,-0.2E+10]) # 3 phase line fault
            ierr[24] = psspy.change_channel_out_file(out) #save to same output file
            ierr[25] = psspy.run(0, 1.17,1,1,1) #run for 10 cycles (~0.17s) with trip
            ierr[26] = psspy.dist_clear_fault(1) # clears fault     
        elif disturbance_type == "bus_fault":
            ierr[23] = psspy.dist_3phase_bus_fault(154,0,1, 230.0,[0.0,-0.2E+10]) #Bus Fault
            ierr[24] = psspy.change_channel_out_file(out) #save to same output file
            ierr[25] = psspy.run(0, 1.17,1,1,1) #run for 10 cycles with trip
            ierr[26] = psspy.dist_clear_fault(1) #clears fault
        ierr[27] = psspy.change_channel_out_file(out) #save to same output file
        ierr[28] = psspy.run(0, runtime,1,1,1) #Run for 10 seconds (hopefully steady state)
        ierr[29] = psspy.delete_all_plot_channels() #Delete plot channels to get ready for next simulation
    print(ierr) #CHECK THIS OUTPUT: Makes sure this list is all zeros, others check what the error value means in the API document
    run_output = output.getvalue() #Collects the PSSe Output Bar info, this includes values of states
                                    # and important error messages
    
    #Check for errors written in output file
    current_error = 0
    if "Network not converged" in run_output:
        print('Network not converged') 
        current_error = 1
		#raise SystemExit #this will quit the program, if the program is called within a larger
        # program, like optimization, you will want to stop PSSe from running using this or 
        # have this rerun the program or skip this iteration's results
    elif "NaN" in run_output:
        print("NaN, network is no good")
        current_error = 1
		#raise SystemExit #this will quit the program, if the program is called within a larger
        # program, like optimization, you will want to stop PSSe from running using this or 
        # have this rerun the program or skip this iteration's results
    if current_error == 0 and "INITIAL CONDITIONS CHECK O.K." in run_output:
        print("No errors and initial conditions were good.")
 	
    #Gather the data and output to excel
    data = dyntools.CHNF(out) #getting data from channel.out file
    d,e,z=data.get_data() #gathering data from data in dictionary format
    #z contains the time series output data from each channel
    #e contains the label of each channel
    #d contains the header of the output file (case name, etc.)
    
    return d,e,z
#===================================================================================================================================================	



## Collect and Sort PSSe Output Data ====================================================================================================================================
def SortResults(d,e,z):
# This function takes the outputs from PSSe and assigns the data to pandas dataframes
# according to data type. This function will need to be updated if additional
# channels are added.

    #Initialize dataframes for each output type
    POWR = []
    POWR_index = 1
    FREQ = []
    FREQ_index = 1
    VOLT = []
    VOLT_index = 1 


    
    #Sort by channel type (machine electric power, bus frequency deviation, bus voltage angle, bus voltage)
    #Append dataframe with channel data, and increase index for inserting next channel's data
    for channel in range(1,53): #Check length of 'e' and run for all those channels
        channel_keys = re.split(' |\[|\]',e[channel]) #Parse name of channel
        if channel_keys[0] == 'POWR':
            if len(POWR) == 0:
                POWR = pd.DataFrame(z[channel],columns =[channel_keys[1]], index= z['time'])
            else:
                POWR.insert(POWR_index,channel_keys[1],z[channel],allow_duplicates = True) #location, column name, data, allow duplicates
                POWR_index = POWR_index+1
        elif channel_keys[0] == 'FREQ':
            if len(FREQ) == 0:
                FREQ = pd.DataFrame(z[channel],columns =[channel_keys[1]], index= z['time'])
            else:            
                FREQ.insert(FREQ_index,channel_keys[1],z[channel],allow_duplicates = True) #location, column name, data, allow duplicates
                FREQ_index = FREQ_index+1    
        elif channel_keys[0] == 'VOLT':
            if len(VOLT) == 0:
                VOLT = pd.DataFrame(z[channel],columns =[channel_keys[1]], index= z['time'])
            else:    
                VOLT.insert(VOLT_index,channel_keys[1],z[channel],allow_duplicates = True) #location, column name, data, allow duplicates
                VOLT_index = VOLT_index+1   

    

    #Sort DataFrames by Column    
    POWR = POWR.sort_index(axis = 1)  
    FREQ = FREQ.sort_index(axis = 1) 
    VOLT = VOLT.sort_index(axis = 1) 

    
    #Machine electric power for just GVEA POI is not just GVEA machines, pull those out here
    #POWR.iloc[:,0:3]

    return POWR,FREQ,VOLT
#===================================================================================================================================================


## Plot Output Results ===================================================================================================
def PlotResults(plot_file,POWR,FREQ,VOLT,excel_file,left_limit,right_limit):
    
    png_POWR = r"""C:\Users\lebredeson\Documents\PSSe_Demo\Outputs\POWR_%s""" %plot_file
    png_FREQ = r"""C:\Users\lebredeson\Documents\PSSe_Demo\Outputs\FREQ_%s""" %plot_file
    png_VOLT = r"""C:\Users\lebredeson\Documents\PSSe_Demo\Outputs\VOLT_%s""" %plot_file
    xlsx = r"""C:\Users\lebredeson\Documents\PSSe_Demo\Outputs\%s""" %excel_file

    plt.figure
    POWR.plot(linewidth=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('MW')
    plt.xlim(left=left_limit,right = right_limit) 
    plt.savefig(png_POWR, format='png', dpi=600)
    plt.figure
    FREQ.plot(linewidth=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('Devation from 60Hz')
    plt.xlim(left=left_limit,right = right_limit)
    plt.savefig(png_FREQ, format='png', dpi=600)
    plt.figure
    VOLT.plot(linewidth=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('V (pu)')
    plt.xlim(left=left_limit,right = right_limit)
    plt.savefig(png_VOLT, format='png', dpi=600)  
    
    #Save Data to excel file
    with pd.ExcelWriter(xlsx) as writer:  
        POWR.to_excel(writer, sheet_name='POWR')
        FREQ.to_excel(writer, sheet_name='FREQ')
        VOLT.to_excel(writer, sheet_name='VOLT')

#=============================================================================================================   

## Run Dynamic Simulation, Collect Data, and Plot ================================================================================================================================

if __name__ == "__main__":

    t0 = time.time()    
    
    
    #Run Test
    
    ## Demo system 'SAVNW' - found in PSSe examples
     
    #Specify model files
    # if using .snp file put 'None' for dyr_file, if using .dyr file put 'None for snp_file
    # .sav file (steady-state base case)
    sav_file_name = 'savnw'
    sav_file = '%s.sav' %sav_file_name
    # .dyr file (dynamic models)
    dyr_file = 'savnw.dyr'
    # .snp file (snapshot which include dynamic models and setup conditions (which can be mid-simulation))
    snp_file = 'None'
    
    #Disturbance Type
    # Type options:
    #"line_fault": Create a three phase to ground fault for transmission line
    #              from bus 154 to bus 3008 at time = 1s. The system is run for 
    #              10 cycles (0.17s) with the fault, then the fault is cleared 
    #              and the system is run until the specified time 'runtime'.
    #"bus_fault": Create a three phase to ground fault at bus 154 at time = 1s.
    #             The system is run for 10 cycles (0.17s) with the fault, then 
    #             the fault is cleared and the system is run until the specified
    #             time 'runtime'.
    #disturbance_type = "bus_fault"
    disturbance_type = "line_fault"
    

    #Pick Channels
    #channel_option = 'only_Area1' #Creat output channels for machines and buses only in Area1 of savnw system
    channel_option = 'All' # Create output channels for all machines and buses in savnw system
    
    #Choose length of total run and x-axis limits for plots
    runtime = 20 #length of simulation run
    left_limit = 0 #left axis-limit for plots
    right_limit = 20 #right axis-limit for plots
    
    #Case Name
    case_name = '%s_%s_%s_%s' %(sav_file_name,disturbance_type,channel_option,runtime)
    
    #Add case name to full file name
    plot_file = 'SP2020_%s.png' %case_name
    excel_file = 'SP2020_%s.xlsx' %case_name
    out_file = 'SP2020_%s.outx' %case_name
    
    # #Run Program for full Railbelt System
    d,e,z = Run_SIM(sav_file,snp_file,dyr_file,out_file,disturbance_type,channel_option,runtime)
    if channel_option == 'All':
        POWR,FREQ,VOLT = SortResults(d,e,z)
        PlotResults(plot_file,POWR,FREQ,VOLT,excel_file,left_limit,right_limit)
    # elif channel_option =="only_Area1":
    #     POWR_GVEA,FREQ_GVEA,ANGL_GVEA,VOLT_GVEA = SortResults_GVEAonly(d,e,z)
    #     PlotResults_GVEAonly(plot_file,POWR_GVEA,FREQ_GVEA,ANGL_GVEA,VOLT_GVEA,excel_file,left_limit,right_limit)
 
   
    #Calculate how long it took to run for reference
    t1 = time.time()
    total  = t1-t0
    print(total) #time for simulation
    
#======================================================================================================================================================	





