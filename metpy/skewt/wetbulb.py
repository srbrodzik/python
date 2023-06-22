import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from sympy import *

sys.path.append('/home/disk/meso-home/brodzik/python/metpy/skewt')

#Function to calculate wetbulb temperature given pressure, dewpoint, and
#temperature. Uses the psychrometric formula from the American
#Meteorological Society glossary. 
#Tetens's formula is used for vapor pressure calculations.

def calc_wetbulb_temp(T,Td,P):

    # Call with:
    # from wetbulb import *
    # wbt = calc_wetbulb_temp(T,Td,P)

    #%matplotlib

    # Test inputs
    #T = 9.6    #degC
    #Td = -15.1 #degC
    #P = 627.9  #hPa

    #T = 12.2
    #Td = 10.9
    #P = 840.6
    
    # Define the expression whose roots we want to find

    Tw=Symbol('Tw')   #Declare the Tw symbol, required for symbolic toolbox operations
    eAct = 6.11*10.**((7.5*Td)/(237.3+Td)); #Actual vapor pressure is calculated from dewpoint

    if T > 0:
        func = lambda Tw: eAct - ( 6.11*10**((7.5*Tw)/(237.3+Tw))-6.60*10**(-4)*(1+0.00115*Tw)*P*(T-Tw) )
    else:
        func = lambda Tw: eAct - ( 6.11*10**((7.5*Tw)/(237.3+Tw))-5.82*10**(-4)*(1+0.00115*Tw)*P*(T-Tw) )
    
    # Plot it to figure out initial guess

    #Tw = np.linspace(-100, 50, 201)

    #plt.plot(Tw, func(Tw))
    #plt.xlabel("Tw")
    #plt.ylabel("expression value")
    #plt.grid()
    #plt.show()

    # Use the numerical solver to find the roots

    ############### Need to see if this initial guess always works #####################
    Tw_initial_guess = 20.
    Tw_solution = fsolve(func, Tw_initial_guess)

    #print "The solution is Tw = %f" % Tw_solution
    #print "at which the value of the expression is %f" % func(Tw_solution)

    return Tw_solution[0]

def plot_wb(df, out_path, out_fname, figtitle, vlim):

    debug = False

    if debug:
        print('FOR WETBULB:')
        print('   out_path  =',out_path)
        print('   out_fname =',out_fname)
        print('   figtitle  =',figtitle)
        print('   vlim      =',vlim)

    # Get wet bulb temperature
    df['wbt'] = df.apply(lambda row: calc_wb(row), axis=1)

    # Create plot
    make_vertical_T_plot(df,out_path,out_fname,figtitle,vlim)
    
def calc_wb(row):
     return calc_wetbulb_temp(row['temperature'],row['dewpoint'],row['pressure'])
    
def make_vertical_T_plot(df,out_path,out_fname,figtitle,vlim):

    #Takes in the sounding data with the calculated wet bulb temperature, vertical boundries for the y-axis, and saves a figure
    #of the vertical temperature with wet bulb profile.
    
    #Gets vectors of the T,Tw,and H variables
    #T = df['temperature'].values * units.degC
    #Tw = df['wbt'].values * units.degC
    #H = df['height'].values * units.m
    T = df['temperature'].values
    Tw = df['wbt'].values
    H = df['height'].values
   
    #Creates the main plotting axis (left)
    fig, axL = plt.subplots()
    
    axL.set_xlabel("Degrees (C)$^\circ$")
    axL.set_ylabel('Height (km)')
    axL.set_ylim(0,int(vlim))
    
    #Plots T and Tw vs H in km
  
    axL.axvline(0,color = 'c', linestyle = 'dashed')
    
    axL.plot(Tw,H/1000.0,'b',linewidth = 2)
    axL.plot(T,H/1000.0,'r',  linewidth = 3)
    
    #Array of zeros for filling the region where T is above freezing
    x=np.zeros(len(T))
    axL.fill_betweenx(H/1000.0,x,T, where= T > x, facecolor = 'lightpink')

    axL.grid()
    
    #Creating the right axis labeled in feet
    axR = axL.twinx()
    axR.set_ylabel('Height (ft\')') 
    axR.set_ylim(0,int(vlim)) #Maintaining the vertical limit to keep both axes equal
    axR.set_yticks(np.arange(0,int(vlim)+1)) #Y ticks corresponding to the left axis (may be changed if the left axis tick marks are specified)
    axR.set_yticklabels(int(np.round(i*KM2FEET)) for i in np.arange(0,int(vlim)+1)) #Converts km into ft and rounds to the nearest foot
    
    #Creates the legend on the left axis plots (ALL)
    axL.legend(["0 C$^\circ$","Wet-Bulb Temperature","Temperature"])
    fig.tight_layout()
    
    #title = 'Sounding for '+site+' at '+date+'/'+time+'Z'
    title = figtitle
    fig.suptitle(title,y=0.99)
    #Saves the figure so no need to return anything
    fig.savefig(out_path+'/'+out_fname)

