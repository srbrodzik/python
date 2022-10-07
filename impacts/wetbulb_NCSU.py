#Function to calculate wetbulb temperature given pressure, dewpoint, and
#temperature. Uses the psychrometric formula from the
#notes to MEA312: Atmospheric Thermodynamics at NC State. 
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

    # Define the expression whose roots we want to find

    epsilon = 0.622
    Lv = 2.5 * (10**6)     #J/kg
    Cp = 1005         #J/kg
    psychro = (Cp*P)/(epsilon*Lv); #Psychrometric constant
    Tw=Symbol('Tw')   #Declare the Tw symbol, required for symbolic toolbox operations
    #T, P, e should be inputs
    eAct = 6.11*10.**((7.5*Td)/(237.3+Td)); #Actual vapor pressure is calculated from dewpoint

    func = lambda Tw: eAct - ( 6.11*10**((7.5*Tw)/(237.3+Tw))-psychro*(T-Tw) )

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
