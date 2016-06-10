### MKS unit ###
import numpy as np
import sys

from subroutines import *
from params import *
import rad_drv_socrates as rad_drv 
from cloud import *
import mks as unit


N_POOL  = 2

# Variables
MM        = np.zeros([N_POOL])
EE        = np.zeros([N_POOL])
SS        = np.zeros([N_POOL])

# Supportive Variables
ll       = np.zeros([N_POOL])
ss       = np.zeros([N_POOL,2])
qq       = np.zeros([N_POOL,2])
zz_satur = np.zeros([N_POOL])
l_satur  = np.zeros([N_POOL])
l_TI_btm = np.zeros([N_POOL])
l_TI_top  = np.zeros([N_POOL])



#=============================================================================
def myfunc( sst_in, qq_TI_c, aa_w, verbose, outfile ):


    #----------------------------------------------------
    # step 0: initial setup : prescribed energy flux
    #----------------------------------------------------
    aa_c      = 1. - aa_w
    FF_m_s    = FF_m_s_0/aa_c   # W m^-2
    LF_m_q    = LF_m_q_0/aa_c  # W m^-2
    FF_m_q    = LF_m_q/LL   # W m^-2

    #----------------------------------------------------
    # step 1: SST([N_POOL]) is passed
    #         temperature at the bottom of atmosphere is specified.
    #----------------------------------------------------
    SST = np.array([ sst_in[0], sst_in[0]-sst_in[1]**2. ])
    Tsurf = SST - dT_surf

    #----------------------------------------------------
    # step 2: compute statistic energy at the bottom
    #----------------------------------------------------
    qq.T[0] = rh2q( Tsurf, P_surf, RH_surf )
    ss.T[0] = find_s( Tsurf, 0. )

    #----------------------------------------------------
    # step 3: find saturation level
    #----------------------------------------------------
    for pool in xrange(2):
        zz_satur[pool] = find_saturation_level( Tsurf[pool], qq[pool][0] )

    #----------------------------------------------------
    # step 4: set up altitude grid
    #----------------------------------------------------
    zz_layers = np.r_[ np.linspace ( 0., zz_satur[1], 20 )[1:-1],
                       np.linspace( zz_satur[1], zz_satur[0], 4 )[:-1],
                       np.linspace( zz_satur[0], zz_TI_btm[1], 30 )[:-1], 
                       np.linspace( zz_TI_btm[1], zz_TI_top[1], 20 )[:-1], 
                       np.linspace( zz_TI_top[1], zz_TI_btm[0], 20 )[:-1], 
                       np.linspace( zz_TI_btm[0], zz_TI_top[0], 20 )[:-1], 
                       np.linspace( zz_TI_top[0], 60000, 140 ) ]
    for pool in xrange(2):
        l_satur[pool]  = np.argmin( abs( zz_layers - zz_satur[pool] ) )        
        l_TI_btm[pool] = np.argmin( abs( zz_layers - zz_TI_btm[pool] ) )
        l_TI_top[pool] = np.argmin( abs( zz_layers - zz_TI_top[pool] ) )

    #----------------------------------------------------
    # step 5: introduce parameters
    #----------------------------------------------------
    N_LAYER = len( zz_layers )
    #---------------
    if ( N_LAYER % 2 == 0 ):
        print 'ERROR! Set odd number of grids'
        sys.exit()
    #---------------
    TT_layers = np.zeros([N_POOL,N_LAYER])
    qq_layers = np.zeros([N_POOL,N_LAYER])
    RH_layers = np.zeros([N_POOL,N_LAYER])
    PP_layers = np.zeros([N_POOL,N_LAYER])

    #----------------------------------------------------
    # step 6: set the profile below saturation level
    #----------------------------------------------------
    # temperature, specific humidity
    for pool in xrange(2):
        TT_layers[pool][:l_satur[pool]+1] = find_Tprof_dryadiabat( Tsurf[pool], zz_layers[:l_satur[pool]+1] )
        qq_layers[pool] = np.tile( qq[pool][0], N_LAYER )
    # pressure
    pool = 0 # warm pool
    PP_layers[pool][:l_satur[pool]+1] = find_Pprof( zz_layers[:l_satur[pool]+1], TT_layers[pool][:l_satur[pool]+1], P_surf, Tsurf[pool] )

    #----------------------------------------------------
    # step 7: warm pool temperature profile above 
    #         saturation level (moist adiabat) and 
    #         find the tropopause
    #----------------------------------------------------
    pool = 0 # warm pool
    TT_layers[pool][l_satur[pool]:], PP_layers[pool][l_satur[pool]:], zz_strato, pp_strato = find_Tprof_moistadiabat( TT_layers[pool][l_satur[pool]], 
                                                                                                                     PP_layers[pool][l_satur[pool]], 
                                                                                                                     zz_layers[l_satur[pool]:] )
    l_strato = np.argmin( abs( zz_layers - zz_strato ) )


    #----------------------------------------------------
    # step 9: above TI, temperature profile of cold pool 
    #         is same as that of warm pool (without TI 
    #         correction, i.e., moist adiabat)
    #----------------------------------------------------
    TT_layers[1][l_TI_top[1]:] = TT_layers[0][l_TI_top[1]:]
    PP_layers[1] = PP_layers[0]


    #----------------------------------------------------
    # step 10: specific humidity and static energy 
    #          above tropopause (it is saturated there!)
    #----------------------------------------------------
    qq.T[1] = np.tile( rh2q( T_strato, pp_strato, 1.0 ), 2 )
    for pool in xrange(2):
        qq_layers[pool][l_strato:] = np.tile( qq[pool][1], N_LAYER - l_strato )
        ss[pool][1] = find_s( T_strato, zz_strato )


    #----------------------------------------------------
    # step 11: moisture profile above the TI
    #----------------------------------------------------
    pool = 0 # warm pool
    RH_layers[pool] = find_RHprof_linear_with_P( RH_surf, PP_layers[pool], ALPHA )
    qq_layers[pool][l_TI_top[pool]:] = rh2q( TT_layers[pool][l_TI_top[pool]:], PP_layers[pool][l_TI_top[pool]:], RH_layers[pool][l_TI_top[pool]:] )
    pool = 1 # cold pool
    RH_TI_c = q2rh( TT_layers[pool][l_TI_top[pool]], PP_layers[pool][l_TI_top[pool]], qq_TI_c )
    RH_layers[pool][l_TI_top[pool]:] = find_RHprof_linear_with_P( RH_TI_c, PP_layers[pool][l_TI_top[pool]:], ALPHA )
    qq_layers[pool][l_TI_top[pool]:] = rh2q( TT_layers[pool][l_TI_top[pool]:], PP_layers[pool][l_TI_top[pool]:], RH_layers[pool][l_TI_top[pool]:] )

    #----------------------------------------------------
    # step 12: T profiles between saturation level and 
    #          TI with mixing theory
    #----------------------------------------------------
    for pool in xrange(2):
        TT_layers[pool][l_satur[pool]:l_TI_btm[pool]+1], qq_layers[pool][l_satur[pool]:l_TI_btm[pool]+1] = find_Tprof_TI( zz_layers[l_satur[pool]:l_TI_btm[pool]+1], 
                                                                                                                          PP_layers[pool][l_satur[pool]:l_TI_btm[pool]+1], 
                                                                                                                          PP_layers[pool][l_satur[pool]], 
                                                                                                                          PP_layers[pool][l_TI_top[pool]], 
                                                                                                                          TT_layers[pool][l_satur[pool]], 
                                                                                                                          TT_layers[pool][l_TI_top[pool]], 
                                                                                                                          qq_layers[pool][l_satur[pool]], 
                                                                                                                          qq_layers[pool][l_TI_top[pool]], 
                                                                                                                          BETA )

    #----------------------------------------------------
    # step 13: profile within TI simply by linear interp.
    #----------------------------------------------------
    for pool in xrange(2):
        TT_layers[pool][l_TI_btm[pool]:l_TI_top[pool]+1] = connect_linear( zz_layers[l_TI_btm[pool]:l_TI_top[pool]+1], 
                                                                          TT_layers[pool][l_TI_btm[pool]], 
                                                                          TT_layers[pool][l_TI_top[pool]] )
        qq_layers[pool][l_TI_btm[pool]:l_TI_top[pool]+1] = connect_linear( zz_layers[l_TI_btm[pool]:l_TI_top[pool]+1], 
                                                                          qq_layers[pool][l_TI_btm[pool]], 
                                                                          qq_layers[pool][l_TI_top[pool]] )
                                                                          
    #----------------------------------------------------
    # step 14: radiation !!
    #----------------------------------------------------
    RR = rad_drv.get_R( SST, zz_satur, zz_strato, zz_layers, TT_layers, PP_layers, qq_layers, MU_atm, MU_H2O, SOL, FACTOR, ALBEDO )
    dR_add_w = cloud_forcing_w( aa_w, zz_layers, [0, zz_satur[0], zz_TI_top[0], zz_strato] )
    dR_add_c = cloud_forcing_c( aa_c )
    RR[0] += dR_add_w + dR_all
    RR[1] += dR_add_c + dR_all

       
    #----------------------------------------------------
    # step 15: obtain Mw and Mc
    #----------------------------------------------------
    qq_TI_w = qq_layers[0][l_TI_top[0]]
    MM[0] = find_M( RR[0], ss[0], qq[0][1], qq_TI_w, 0. )
    MM[1] = find_M( RR[1], ss[1], qq[1][1], qq_TI_c, FF_m_s )

    #----------------------------------------------------
    # step 16: obtain Sw and Sc
    #----------------------------------------------------
    SS[0] = find_Sw( MM, RR, ss, ss[0][0], aa_w, aa_c )
    SS[1] = find_Sc( RR ) # ss_Bplus is equal to ss (z_B)
    
    #----------------------------------------------------
    # step 17: Obtain Ew and Ec
    #----------------------------------------------------
    EE[0] = find_Ew( MM, qq_layers, qq_TI_w, aa_w, aa_c ) #??
    EE[1] = find_Ec( MM, qq_layers, qq_TI_c, FF_m_q ) #??
    
    #----------------------------------------------------
    # step 18: balanck check
    #----------------------------------------------------
    dE_w = check_surfbalance( RR[0][0], EE[0], SS[0], FF_o_w )
    dE_c = check_surfbalance( RR[1][0], EE[1], SS[1], FF_o_c )

    residual = np.array([ dE_w, dE_c ])


    #----------------------------------------------------------------
    if ( verbose == True ):

        with open( outfile, 'a') as f:

            f.write( str(aa_w) + "\t" )
            f.write( str(qq_TI_c*1e3) + "\t" )

            f.write( str(SST[0]-273.15) + "\t" )
            f.write( str(SST[1]-273.15) + "\t" )

            f.write( str(zz_strato*1e-3) + "\t" )

            f.write( str(RR[0][0]) + "\t" )
            f.write( str(RR[0][1]) + "\t" )
            f.write( str(RR[0][2]) + "\t" )
            f.write( str(RR[1][0]) + "\t" )
            f.write( str(RR[1][1]) + "\t" )
            f.write( str(RR[1][2]) + "\t" )
            f.write( str(MM[0])  + "\t" )
            f.write( str(MM[1])  + "\t" )
            f.write( str(SS[0])  + "\t" )
            f.write( str(SS[1])  + "\t" )
            f.write( str(LL*EE[0])  + "\t" )
            f.write( str(LL*EE[1])  + "\t" )
            f.write( str(residual[0]) + "\t" )
            f.write( str(residual[1]) + "\t" )
    #----------------------------------------------------------------

    return residual



