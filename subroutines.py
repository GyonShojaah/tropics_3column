import numpy as np
from params import *
from scipy.optimize import fsolve
import errors

#--------------------------------------------------------
def CCpressure( temp ):
    """
    Approximate Clausius-Clapeyron equation.
    Ref) A Short course in Cloud Physics eq(2.17)
    Latent heat is assumed to be constant.
    """
#    TTc = (TT-273.0)                                                                     
#    return 611.2*np.exp(17.67*TTc/(TTc+243.5))
    return 2.53e11 * np.exp( -5.42e3 / temp ) # Pa


#--------------------------------------------------------
def inv_CCpressure( p_H2O ):
    """
    Approximate Clausius-Clapeyron equation.
    Ref) A Short course in Cloud Physics eq(2.17)
    Latent heat is assumed to be constant.
    """
#    TTc = (TT-273.0)                                                                     
#    return 611.2*np.exp(17.67*TTc/(TTc+243.5))
#    return 2.53e11 * np.exp( -5.42e3 / temp ) # Pa
#    np.log( p_H2O / 2.53e11 ) = - 5.42e3 / temp
    return 5.42e3 / np.log( 2.53e11 / p_H2O )

#--------------------------------------------------------
def find_Tprof_dryadiabat( t_surf, z_layers ):
    """
    tempmerature profile assuming dry adiabat
    """
    Gamma = gg / c_p
    T_layers = t_surf - Gamma * z_layers
    return T_layers


#--------------------------------------------------------
def find_s( T, z ):
    """
    compute dry static energy
    """
    print "c_p * T, gg * z", c_p * T, gg * z
    s = c_p * T + gg * z
    return s


#--------------------------------------------------------
def find_saturation_level( t_surf, q_surf ):
    """
    compute saturation level
    """
    ratio = q_surf / MU_H2O * MU_atm
    def nulfunc( temp ):
        p_H2O = ratio * P_surf * ( temp / t_surf )**DELTA
        p_sat = CCpressure( temp )
        return p_H2O - p_sat
    t_sat = fsolve( nulfunc, t_surf )[0]
    Gamma = gg / c_p
    z_sat = ( t_surf - t_sat )/Gamma
    return z_sat


#--------------------------------------------------------
def find_saturation_point( temp, pres, qq ):
    """
    compute saturation level
    """
    ratio = qq / MU_H2O * MU_atm
    def nulfunc( t_sat ):
        p_H2O_1 = ratio * pres * ( t_sat / temp )**DELTA
        p_H2O_2 = CCpressure( t_sat )
        return p_H2O_1 - p_H2O_2
    t_sat = fsolve( nulfunc, temp )[0]
    p_sat = pres * ( t_sat / temp )**DELTA
    p_H2O = CCpressure( t_sat )

    return t_sat, p_sat



#--------------------------------------------------------
def find_RHprof_linear_with_P( rh0, t0, z_layers, alpha ):
    """
    compute profile of relative humidity, assuming that it changes linearly with pressure
    ****IT IS ONLY APPROXIMATION !!****
    in Miller's code, after Manabe and Wetherald, (JAS,1967)
         rh = rhs*(pm-0.02*psurf)/(psurf-0.02*psurf)
     =>  rh = rhs*(pm/psurf - 0.02)/(0.98)
     =>  rh = rhs*(12.5)*pm/psurf - rhs*0.25
     => rh/rh0 = (12.5)*pm/psurf - ( ALPHA - 1. )
    RH = RH0 - ALPHA * ( 1 - exp(-z/H) )
    SH  = ( RR * temp )/( MU_atm * gg )
    """
    ScaleHeight = ( RR * t0 ) / ( MU_atm * gg )
    rh_layers = rh0 * ( alpha * np.exp(-(z_layers-z_layers[0])/ScaleHeight) - ( alpha - 1. ) )
    rh_layers[np.where( rh_layers < 0. )] = np.zeros( len(np.where( rh_layers < 0. )) )
    return rh_layers


#--------------------------------------------------------
def find_Tprof_moistadiabat( t0, p0, z_layers ):
    """
    temperature profile assuming moist adiabat
    ****IT IS ONLY APPROXIMATION !!****
    """
    t_layers = np.zeros_like( z_layers ) + t0
    z_layers_diff = np.diff( z_layers )
    ScaleHeight = ( RR * t0 ) / ( MU_atm * gg )
    p_layers = p0 * np.exp( ( z_layers[0] - z_layers )  / ScaleHeight )
    epsilon = MU_H2O / MU_atm
    l_strato = len( z_layers )
    for ll in xrange( len(z_layers)-1 ): 
#        mixing_ratio = CCpressure( t_layers[ll] ) / p_layers[ll]
        mixing_ratio_mass = ( CCpressure( t_layers[ll] )*MU_H2O ) / ( p_layers[ll] * MU_atm )
        fac1 = 1.  + ( LL * mixing_ratio_mass )/( Rs_dry * t_layers[ll] ) 
        fac2 = c_p + ( LL**2 * mixing_ratio_mass * epsilon )/( Rs_dry * t_layers[ll]**2 )
        dtdz = - ( fac1 / fac2 ) * gg
        t_layers[ll+1] = t_layers[ll] + dtdz * z_layers_diff[ll]
        dpdz =  - ( MU_atm * p_layers[ll] / ( RR * t_layers[ll] ) ) * gg
        # correction
        if ( t_layers[ll+1] < T_strato ):
            t_layers[ll+1] = T_strato
            if ( l_strato > ll + 1 ):
                l_strato = ll + 1

    return t_layers, p_layers, l_strato


#--------------------------------------------------------
def connect_linear( x_grid, y0, y1 ):

    x0 = x_grid[0]
    x1 = x_grid[-1]
    kk = (y1-y0)/(x1-x0)
    y_grid = kk * ( x_grid - x0 ) + y0
    return y_grid


#--------------------------------------------------------
def find_Pprof( z_layers, t_layers, p0 ):
    """
    compute profile of pressure given the temperature profile

    dP/dz = - rho g = - ( MU_atm * P / ( R T ) ) g
    log P  = - ( MU_atm g / R ) * integrate[ 1. / T ] dz  ]
    """
    inv_t_layers = 1. / t_layers
    dlogP = np.zeros( len(z_layers) )
    for ll in xrange( len(z_layers) ):
        dlogP[ll] = np.trapz( inv_t_layers[:ll+1], x=z_layers[:ll+1] )
    p_layers = p0 * np.exp( -1.0 *( MU_atm * gg / RR ) * dlogP ) 
    return p_layers    


#--------------------------------------------------------
def find_theta_e( temp, pres, qq ):

    theta = temp * ( P_surf / pres )**( 1./DELTA )
    t_sat, p_sat = find_saturation_point( temp, pres, qq )
    theta_e = theta * np.exp( LL * qq / ( c_p * temp ) )
    return theta_e


#--------------------------------------------------------
def find_T_from_theta_e( theta_e, pres, qq, t0 ):

    def nulfunc( temp ):
        return find_theta_e( temp, pres, qq ) - theta_e
    TT = fsolve( nulfunc, t0 )[0]
    return TT


#--------------------------------------------------------
def find_Tprof_TI( z_layers, p_layers, p0, p1, T0, T1, q0, q1, beta ):

    p_sat_layers = beta * ( p_layers - p_layers[0] ) + p_layers[0]

    theta0_e = find_theta_e( T0, p0, q0 )
    theta1_e = find_theta_e( T1, p1, q1 )

    t_layers = np.zeros( len(z_layers) )
    q_layers = np.zeros( len(z_layers) )

    for zi in xrange( len(z_layers) ):

        def func( xx, p_sat ):
            theta_e = theta0_e * xx + theta1_e * ( 1. - xx )
            qq      = q0 * xx + q1 * ( 1. - xx )
            temp1   = find_T_from_theta_e( theta_e, p_layers[zi], qq, T0 )
            temp2   = inv_CCpressure( p_sat_layers[zi] * MU_atm * qq / MU_H2O )*( p_layers[zi] / p_sat_layers[zi] )**( 1./DELTA )
            return temp1 - temp2

        x_ans = fsolve( func, 1., args=p_sat_layers[zi] )[0]

        q_layers[zi] = q0 * x_ans + q1 * ( 1. - x_ans )
        theta_e = theta0_e * x_ans + theta1_e * ( 1. - x_ans )
        t_layers[zi] = find_T_from_theta_e( theta_e, p_layers[zi], q_layers[zi], T0 )

    return t_layers, q_layers



#--------------------------------------------------------
def find_saturation_q( pres, temp, p_sat ):
    """
    compute saturation level
    """
    t_sat = temp * ( p_sat / pres )**( 1./DELTA )
    p_H2O = CCpressure( t_sat )
    q =  ( p_H2O * MU_H2O ) / ( p_sat * MU_atm )
    return q


#--------------------------------------------------------
def rh2q( temp, pres, rh ):
    p_H2O = CCpressure( temp ) * rh
    q = ( p_H2O * MU_H2O ) / ( pres * MU_atm )
    return q


#--------------------------------------------------------
def q2rh( temp, pres, qq ):

    mixing_ratio = qq * MU_atm / MU_H2O
    p_H2O = pres * mixing_ratio
    rh = p_H2O / CCpressure( temp )
    return rh


#--------------------------------------------------------
def find_tropopause( T_layers ):
    """
    find tropopause and replace the temperature above it by T_strato
    """
    l_strato = np.where( T_layers < T_strato )[0][0]
    for ii in xrange( l_strato, len(T_layers) ):
        T_layers[ii] = T_strato
    return l_strato, T_layers


#--------------------------------------------------------
def find_M( R, s, qT, q_TI, F ):
    """
    find M
    """
    M = ( R[1] - R[2] - F )/( ( s[1] + LL*qT ) - ( s[0] + LL*q_TI ) )
    print "s[1], LL*qT, s[0], LL*q_TI", s[1], LL*qT, s[0], LL*q_TI
    return M


#--------------------------------------------------------
def find_Sw( M, R, s, s_Bplus_w, aa_w, aa_c ): 
    Mw, Mc = M
    Sw = ( aa_c/aa_w )*Mc*( s[0][0] - s[1][0] ) + ( R[0][0] - R[0][1] )
    return Sw


#--------------------------------------------------------
def find_Sc( M, R, s, s_Bplus_c ): 
    Mw, Mc = M
    Sc = R[1][0] - R[1][1]
    return Sc


#--------------------------------------------------------
def find_Ew( M, q_layers, l_TI, aa_w, aa_c ):
    Mw, Mc = M
    Ew = Mw * ( q_layers[0][0] - q_layers[0][l_TI] ) + ( aa_c / aa_w ) * Mc * ( q_layers[0][0] - q_layers[1][0] )
    return Ew


#--------------------------------------------------------
def find_Ec( M, q_layers, l_TI, F_mq ):
    Mw, Mc = M
    Ec = Mc * ( q_layers[1][0] - q_layers[1][l_TI] ) - F_mq
    return Ec


#--------------------------------------------------------
def check_surfbalance( R0, E, S, Fo ):
    check = R0 - LL * E - S + Fo
    return check



