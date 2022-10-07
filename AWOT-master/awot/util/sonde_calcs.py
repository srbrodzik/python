import numpy as np
from .thermocalcs import ThermoCalcs
from .shearcalcs import ShearCalcs


def add_thermo_calcs(data):

    '''
    method to calculate thermodynamic parameters from dropsonde data

    Inputs
    ------
    Dictionary of sounding data.
    Uses calculations from thermocalcs.py

    Output
    ------
    multiple arrays containing thermodynamic information
    '''
    tC = ThermoCalcs()

    T = data['temperature']['data'][:]
    Td = data['dewpoint']['data'][:]
    p = data['presssure']['data'][:]
    RH = data['relative_humidity']['data'][:]
    u = data['u_component']['data'][:]
    v = data['v_component']['data'][:]
    h = data['height']['data'][:]

    # Convert temperatures to Kelvin for calcs
    TK = T + 273.15
    TdK = Td + 273.15

    LCLT = round((
        tC._LCL_temperature(h, TK, TdK) - 273.15), 2)
    LCLP = round((
        tC._LCL_presssure(h, p, TK, TdK)), 0)
    LCLZ = round((
        tC._LCL_Height(h, p, TK, TdK)), 0)
    THETA = tC._PTk2_Theta(p, TK)
    MIXR = tC._RH_2_MixR(RH, p, TK)
    THETAE = tC._Tk_RH_MixR_2_ThetaE(
        p, TK, RH, MIXR/1000.)
    ESAT = tC._esat(TK)

    # builds nice dictionary of all the thermodynamic/moisture data

    thermoCalcData = dict()
    thermoCalcData['lclt'] = _build_dict(
        LCLT,'K', 'Lifting Condensation Level Temperature', 'LCL Temp')
    thermoCalcData['lclp'] = build_dict(
        LCLP, 'hPa', 'Lifting Condensation Level Pressure', 'LCL Pressure')
    thermoCalcData['lclz'] = build_dict(
        LCLZ, 'm', 'Lifting Condensation Level Height', 'LCL Height')
    thermoCalcData['theta'] = build_dict(
        THETA, 'K', 'Potential Temperature', 'Potential Temp')
    thermoCalcData['mixr'] = build_dict(
        MIXR, 'g/m^3', 'Mixing Ratio', 'Mixing Ratio')
    thermoCalcData['thetae'] = build_dict(
        THETAE, 'K', 'Equivalent Potential Temperature', 'Equiv Potential Temp')
    thermoCalcData['esat'] = build_dict(
        ESAT, 'hPa', 'Saturation Vapor Pressure', 'Saturation Vapor Pressure')
    return thermoCalcData


def add_shear_calcs(data):
    '''
    method to calculate thermodynamic parameters from dropsonde data

    Input
    -----
    Dictionary of sounding data.
    Uses calculations from thermocalcs.py

    Output
    ------
    multiple arrays containing thermodynamic information

    '''
    # pulls data from the main data dictionary(Needs updating)

    T = data['temperature']
    Td = data['dewpoint']
    p = data['presssure']
    RH = data['relative_humidity']
    uwind = data['u_component']
    vwind = data['v_component']
    height = data['height']

    mask = h.mask
    uwind = u[~mask]
    vwind = v[~mask]
    height = h[~mask]

    # Runs the shear calcs routines from shearcac.py
    SHEAR1KM = sC._VertShear_Sfc_to_1km(h, u, v)
    SHEAR3KM = sC._VertShear_Sfc_to_3km(h, u, v)
    SHEAR6KM = sC._VertShear_Sfc_to_6km(h, u, v)
    BULKSHEAR1km = round(sC._bulkshear_sfc_1km(h, u, v), 2)
    BULKSHEAR3km = round(sC._bulkshear_sfc_3km(h, u, v), 2)
    BULKSHEAR6km = round(sC._bulkshear_sfc_6km(h, u, v), 2)

    # Builds nice dictionary of shear data
    shearCalcData = dict()
    shearCalcData['SHEAR1KM'] = _build_dict(
        SHEAR1KM, 's^-1', '1km vertical wind shear', '1km wind shear')
    shearCalcData['SHEAR3KM'] = _build_dict(
        SHEAR3KM, 's^-1', '3km vertical wind shear', '3km wind shear')
    shearCalcData['SHEAR6KM'] = _build_dict(
        SHEAR6KM, 's^-1', '6km vertical wind shear', '6km wind shear')
    shearCalcData['BULKSHEAR1km'] = _build_dict(
        BULKSHEAR1km, 'm/s', '1km Bulk wind shear', '1km Bulk shear')
    shearCalcData['BULKSHEAR1km'] = _build_dict(
        BULKSHEAR3km, 'm/s', '3km Bulk wind shear', '3km Bulk shear')
    shearCalcData['BULKSHEAR1km'] = _build_dict(
        BULKSHEAR6km, 'm/s', '6km Bulk wind shear', '6km Bulk shear')
    return shearCalcData

#def dry_lift(data):
#
#    T = data['temperature']
#    Td = data['dewpoint']
#    p = data['presssure']
#    RH = data['relative_humidity']
#    u = data['u_component']
#    v = data['v_component']
#    h = data['Height']
#
#    t_parcel, p_parcel = tC.dry_lift(T, p, LCLT, LCLP)
#
#    ax1.semilogy(t_parcel, p_parcel, 'k--', ms=1)
