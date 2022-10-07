#Function to calculate wetbulb temperature given pressure, dewpoint, and
#temperature. Uses the psychrometric formula from the American
#Meteorological Society glossary. 
#Tetens's formula is used for vapor pressure calculations.

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from sympy import *

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
