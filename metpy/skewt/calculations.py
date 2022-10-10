# Borrowed equations from Brody's old skewplot.py code and subroutines 

import numpy as np
from numpy import exp
import math

Epsilon=0.622         # Epsilon=R_s_da/R_s_v; The ratio of the gas constants
verbose = False

def calcs(df):

        h0 = df['height'][np.argmin(abs(df['temperature']-0))]
        h10 = df['height'][np.argmin(abs(df['temperature']+10))]
        h20 = df['height'][np.argmin(abs(df['temperature']+20))]
        h30 = df['height'][np.argmin(abs(df['temperature']+30))]
        h40 = df['height'][np.argmin(abs(df['temperature']+40))]
        if verbose:
                print('h0 =',h0,'h10 =',h10,'h20 =',h20,'h30 =',h30,'h40 =',h40)

        lcl_height = 216*(df['temperature'][0]-df['dewpoint'][0])
        if verbose:
                print('lcl_height =',lcl_height)
        
        wcd = h0 - lcl_height
        if verbose:
                print('wcd =',wcd)
        
        # SRB: removed this formula for mixing ratio and reverted to orig definition
        #mr = MixRatio(SatVap(df['dewpoint'][0]), df['pressure']*100)
        mr = MixRatio(SatVap(df['dewpoint']), df['pressure']*100)
        pdiff = -1.0*np.diff(df['pressure'])
        # precipitable water
        pw = np.sum(np.array([mr[_]*100.0*pdiff[_]/9.8 for _ in range(pdiff.shape[0]) ]))
        if verbose:
                print('pw =',pw)
        
        # crude estimate of wet bulb temp
        ##tw = df['temperature'][0]- (df['temperature'][0]-df['dewpoint'][0])/3.
        
        # now some shear calculations
        wspd6km = df['speed'][np.argmin(abs(df['height']-6000))]
        wdir6km = df['direction'][np.argmin(abs(df['height']-6000))]

        udiff = wspd6km*np.cos(np.radians(270-wdir6km)) - df['speed'][3]*np.cos(np.radians(270-df['direction'][3]))
        vdiff = wspd6km*np.sin(np.radians(270-wdir6km)) - df['speed'][3]*np.sin(np.radians(270-df['direction'][3]))
        # print udiff, vdiff
        shear6km = np.sqrt(udiff**2 + vdiff**2)
        if (math.isnan(shear6km)):
            shear6km = 0
        if verbose:
                print('shear6km =',shear6km)

        # 850mb-200mb Shear
        wspd850mb = df['speed'][np.argmin(abs(df['pressure']-850))]
        wdir850mb = df['direction'][np.argmin(abs(df['pressure']-850))]
        wspd200mb = df['speed'][np.argmin(abs(df['pressure']-200))]
        wdir200mb = df['direction'][np.argmin(abs(df['pressure']-200))]

        udiff = wspd200mb*np.cos(np.radians(270-wdir200mb)) - wspd850mb*np.cos(np.radians(270-wdir850mb))
        vdiff = wspd200mb*np.sin(np.radians(270-wdir200mb)) - wspd850mb*np.sin(np.radians(270-wdir850mb))
        # print udiff, vdiff
        shear850200mb = np.sqrt(udiff**2 + vdiff**2)
        if (math.isnan(shear850200mb)):
            shear850200mb = 0
        if verbose:
                print('shear850200mb =',shear850200mb)

        # SFC-700mb Shear
        wspd700mb = df['speed'][np.argmin(abs(df['pressure']-700))]
        wdir700mb = df['direction'][np.argmin(abs(df['pressure']-700))]

        udiff = wspd700mb*np.cos(np.radians(270-wdir700mb)) - df['speed'][3]*np.cos(np.radians(270-df['direction'][3]))
        vdiff = wspd700mb*np.sin(np.radians(270-wdir700mb)) - df['speed'][3]*np.sin(np.radians(270-df['direction'][3]))
        # print udiff, vdiff
        shear700mb = np.sqrt(udiff**2 + vdiff**2)
        if (math.isnan(shear700mb)):
            shear700mb = 0
        if verbose:
                print('shear700mb =',shear700mb)

        return h0,h10,h20,h30,h40,wcd,pw,lcl_height,shear6km,shear850200mb,shear700mb

def MixRatio(e,p):
        
    """
    Mixing ratio of water vapour
    INPUTS
    e (Pa) Water vapor pressure
    p (Pa) Ambient pressure

    RETURNS
    qv (kg kg^-1) Water vapor mixing ratio`
    """

    return Epsilon*e/(p-e)

def SatVap(tempc,phase="liquid"):
        
    """
    Calculate saturation vapour pressure over liquid water and/or ice.

    INPUTS: 
    tempc: (C)
    phase: ['liquid'],'ice'. If 'liquid', do simple dew point. If 'ice',
    return saturation vapour pressure as follows:

    Tc>=0: es = es_liquid
    Tc <0: es = es_ice


    RETURNS: e_sat  (Pa)

    SOURCE: http://cires.colorado.edu/~voemel/vp.html (#2:
    CIMO guide (WMO 2008), modified to return values in Pa)

    This formulation is chosen because of its appealing simplicity, 
    but it performs very well with respect to the reference forms
    at temperatures above -40 C. At some point I'll implement Goff-Gratch
    (from the same resource).
    """
        
    over_liquid=6.112*exp(17.67*tempc/(tempc+243.12))*100.
    over_ice=6.112*exp(22.46*tempc/(tempc+272.62))*100.

    if phase=="liquid":
        return over_liquid
    elif phase=="ice":
        return where(tempc<0,over_ice,over_liquid)
    else:
        raise NotImplementedError


