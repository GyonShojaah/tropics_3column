#=============================================================================
# Module
#=============================================================================
DC_DTHETA = 5.70
OFFSET    = -55.73
DF_DC     = -1.162
import numpy as np
import mks as unit
from params import *
import sys

def cloud_forcing( t_layers, p_layers ):

    T_700hPa = t_layers[ np.where( p_layers <= 700.0*1e2 )[0][0] ]
    T_surf   = t_layers[ 0 ]
    theta_700hPa = potential_temperature( T_700hPa, 700.0*1e2 )
    theta_surf = potential_temperature( T_surf, p_layers[0] )
    dtheta  = theta_700hPa - theta_surf

#    print "700hPa", T_700hPa, theta_700hPa
#    print "surface", T_surf, theta_surf
    dc = DC_DTHETA * dtheta + OFFSET
    if ( dc < 0. ):
        dc = 0.
    elif ( dc > 100.0 ):
#        print "Cloud Cover Exceeding 100% ! ==> corrected to 100%"
        dc = 100.
#    print ""
#    print "dtheta", dtheta
#    print "dc", dc
#    print ""
    dF = DF_DC * dc - 6.21
    return dF


def cloud_forcing_flag( t_layers, p_layers ):

    T_700hPa = t_layers[ np.where( p_layers <= 700.0*1e2 )[0][0] ]
    T_surf   = t_layers[ 0 ]
    theta_700hPa = potential_temperature( T_700hPa, 700.0*1e2 )
    theta_surf = potential_temperature( T_surf, p_layers[0] )
    dtheta  = theta_700hPa - theta_surf

#    print "700hPa", T_700hPa, theta_700hPa
#    print "surface", T_surf, theta_surf
    dc = DC_DTHETA * dtheta + OFFSET
    flag_dc = 0
    if ( dc < 0. ):
        dc = 0.
    elif ( dc > 100.0 ):
#        print "Cloud Cover Exceeding 100% ! ==> corrected to 100%"
        flag_dc = -1
        dc = 100.
#    print ""
#    print "dtheta", dtheta
#    print "dc", dc
#    print ""
    dF = DF_DC * dc
    return dF, flag_dc

def potential_temperature( temp, pres ):

#    print ( unit.RR / c_p )
    theta = temp * ( P_surf / pres )**( unit.RR / ( c_p * mu_atm )  )
#    print "unit.RR / ( c_p * mu_atm )", unit.RR / ( c_p * mu_atm )
#    sys.exit()
    return theta



def cloud_forcing_w( a_w, z_layers, l_list ):
    

    z_0       = z_layers[int(l_list[0])]
    z_satur   = z_layers[int(l_list[1])]
    z_TI      = z_layers[int(l_list[2])]
    z_strato  = z_layers[int(l_list[3])]

    dR_strato = dR_w_TOA / a_w
    dR_TI     = ( ( dR_w_TOA - dR_w_satur ) / ( z_strato - z_satur ) *  ( z_TI - z_satur ) + dR_w_satur ) / a_w
    dR_0      = dR_w_satur / a_w

    return np.array( [ dR_0, dR_TI, dR_strato ])



def cloud_forcing_c( a_c, z_layers, l_list ):
    
    dR_all = dR_c / a_c
    return np.array( [ dR_all, dR_all, dR_all ])
