### MKS unit ###
import numpy as np
from subroutines import *
from params import *
#import rad_drv_libradtran as rad_drv 
import rad_drv_socrates as rad_drv 
from cloud import *
import mks as unit
import sys

N_LAYER = 299
N_POOL  = 2
zz_layers = np.r_[ 0.0, np.logspace(-1, np.log10(60), N_LAYER)[:-1] * 1e3 ]

# Variables
MM        = np.zeros([N_POOL])
EE        = np.zeros([N_POOL])
SS        = np.zeros([N_POOL])
RR        = np.zeros([N_POOL,3])
TT_layers = np.zeros([N_POOL,N_LAYER])
qq_layers = np.zeros([N_POOL,N_LAYER])

# Supportive Variables
RH_layers = np.zeros([N_POOL,N_LAYER])
PP_layers = np.zeros([N_POOL,N_LAYER])
ll       = np.zeros([N_POOL])
ss       = np.zeros([N_POOL,2])
zz_satur = np.zeros([N_POOL])
l_satur  = np.zeros([N_POOL])


#=============================================================================
# main
#=============================================================================
#if __name__ == "__main__":

def find_R( sst_in, qq_TI_c, aa_w, verbose=False ):

    #----------------------------------------------------
    # step 0a: initial setup : prescribed energy flux
    #----------------------------------------------------
    aa_c      = 1. - aa_w
    FF_m_s    = -8.6/aa_c   # W m^-2
    LF_m_q    = -10.2/aa_c  # W m^-2
    FF_m_q    = LF_m_q/LL   # W m^-2

    #----------------------------------------------------
    # step 0b: initial setup : location of trade inversion
    #----------------------------------------------------
    l_TI_btm = np.zeros(2)
    l_TI_tp  = np.zeros(2)
    l_TI_btm[0] = np.where( zz_layers > zz_TI_btm[0] )[0][0]
    l_TI_btm[1] = np.where( zz_layers > zz_TI_btm[1] )[0][0]
    l_TI_tp[0]  = np.where( zz_layers > zz_TI_tp[0] )[0][0]
    l_TI_tp[1]  = np.where( zz_layers > zz_TI_tp[1] )[0][0]

    #----------------------------------------------------
    # step 1: SST([N_POOL]) is passed
    #         temperature at the bottom of atmosphere is specified.
    #----------------------------------------------------
    SST = np.array([ sst_in[0], sst_in[0]-sst_in[1]**2. ])
    print ""
    print "sst_w, sst_c", SST[0]-273.15, SST[1]-273.15
    print ""

    Tsurf = SST - dT_surf

#    TT_layers.T[0] = SST - dT_surf

    #----------------------------------------------------
    # step 2: compute statistic energy at the bottom
    #----------------------------------------------------
    qq_layers.T[0] = rh2q( Tsurf, P_surf, RH_surf )
#    print "static energy at the bottom"
    ss.T[0]        = find_s( Tsurf, 0. )

    #----------------------------------------------------
    # step 3: find saturation level
    #----------------------------------------------------
    for pool in xrange(2):
        zz_satur[pool] = find_saturation_level( Tsurf[pool], 
                                                qq_layers[pool][0] )
        l_satur[pool]  = np.where( zz_layers > zz_satur[pool] )[0][0]
        
    #----------------------------------------------------
    # step 4: find T profiles below the saturation level
    #----------------------------------------------------
    for pool in xrange(2):
        TT_layers[pool][:l_satur[pool]+1] = find_Tprof_dryadiabat( Tsurf[pool], zz_layers[:l_satur[pool]+1] )
        qq_layers[pool] = np.tile( qq_layers[pool][0], N_LAYER )

    #----------------------------------------------------
    # step 5: find P profile in the warm pool below the saturation level
    #----------------------------------------------------
    pool = 0 # warm pool
    PP_layers[pool][:l_satur[pool]+1] = find_Pprof( zz_layers[:l_satur[pool]+1], 
                                                    TT_layers[pool][:l_satur[pool]+1], 
                                                    P_surf )
    
    #----------------------------------------------------
    # step 6: warm pool
    #           find T/RH profile above the saturation level (moist adiabat)
    #           and also find the tropopause  
    #----------------------------------------------------
    pool = 0 # warm pool
    # approximation HERE!!!!!!!!!
    TT_layers[pool][l_satur[pool]:], PP_layers[pool][l_satur[pool]:], l_strato = find_Tprof_moistadiabat( TT_layers[pool][l_satur[pool]], 
                                                                                                          PP_layers[pool][l_satur[pool]], 
                                                                                                          zz_layers[l_satur[pool]:] )
    l_strato = l_strato + l_satur[pool]
#    qq_layers[pool][l_satur[pool]:] = rh2q( TT_layers[pool][l_satur[pool]:], 
#                                            PP_layers[pool][l_satur[pool]:], 
#                                            np.ones( N_LAYER - l_satur[pool] ) ) 
#

    RH_layers[pool] = find_RHprof_linear_with_P( RH_surf, PP_layers[pool], ALPHA )
    qq_layers[pool][l_TI_tp[pool]:] = rh2q( TT_layers[pool][l_TI_tp[pool]:],
                                            PP_layers[pool][l_TI_tp[pool]:], 
                                            RH_layers[pool][l_TI_tp[pool]:] )


    #----------------------------------------------------
    # step 7: specific humidity above tropopause
    #----------------------------------------------------
    qT = rh2q( TT_layers[pool][l_strato], PP_layers[pool][l_strato], 1.0 )
    for pool in xrange(2):
        qq_layers[pool][l_strato:] = np.tile( qT, N_LAYER - l_strato )


    #----------------------------------------------------
    # step 8: above the trade inversion, the temperature profile of cold pool is same as that of warm pool
    #----------------------------------------------------
    TT_layers[1][l_TI_tp[1]:] = TT_layers[0][l_TI_tp[1]:]
    # just for now
    PP_layers[1] = PP_layers[0]

    #----------------------------------------------------
    # step 10: find moisture profile above the TI
    #----------------------------------------------------
    pool = 1
    RH_TI2 = q2rh( TT_layers[pool][l_TI_tp[pool]], 
                   PP_layers[pool][l_TI_tp[pool]], 
                   qq_TI_c )
    RH_layers[pool][l_TI_tp[pool]:] = find_RHprof_linear_with_P( RH_TI2, PP_layers[pool][l_TI_tp[pool]:], ALPHA )
            
    qq_layers[pool][l_TI_tp[pool]:] = rh2q( TT_layers[pool][l_TI_tp[pool]:],
                                            PP_layers[pool][l_TI_tp[pool]:], 
                                            RH_layers[pool][l_TI_tp[pool]:] )


    #----------------------------------------------------
    # step 10: find T profiles between saturation level and inversion with mixing theory
    #----------------------------------------------------
    for pool in xrange(2):

#        print "z_B, z_TI," , zz_layers[l_satur[pool]], zz_layers[l_TI_tp[pool]]
        TT_layers[pool][l_satur[pool]:l_TI_btm[pool]+1], qq_layers[pool][l_satur[pool]:l_TI_btm[pool]+1] = find_Tprof_TI( zz_layers[l_satur[pool]:l_TI_btm[pool]+1], 
                                                                                                                          PP_layers[pool][l_satur[pool]:l_TI_btm[pool]+1], 
                                                                                                                          PP_layers[pool][l_satur[pool]], 
                                                                                                                          PP_layers[pool][l_TI_tp[pool]], 
                                                                                                                          TT_layers[pool][l_satur[pool]], 
                                                                                                                          TT_layers[pool][l_TI_tp[pool]], 
                                                                                                                          qq_layers[pool][l_satur[pool]], 
                                                                                                                          qq_layers[pool][l_TI_tp[pool]], 
                                                                                                                          BETA )
    #----------------------------------------------------
    # step 11: determine temperature in TI simply by linear interpolation
    #----------------------------------------------------
    for pool in xrange(2):

        TT_layers[pool][l_TI_btm[pool]:l_TI_tp[pool]+1] = connect_linear( zz_layers[l_TI_btm[pool]:l_TI_tp[pool]+1], 
                                                                          TT_layers[pool][l_TI_btm[pool]], 
                                                                          TT_layers[pool][l_TI_tp[pool]] )
        qq_layers[pool][l_TI_btm[pool]:l_TI_tp[pool]+1] = connect_linear( zz_layers[l_TI_btm[pool]:l_TI_tp[pool]+1], 
                                                                          qq_layers[pool][l_TI_btm[pool]], 
                                                                          qq_layers[pool][l_TI_tp[pool]] )

    




#    if ( verbose == True ):
#        for ii in xrange( N_LAYER ):
#            print zz_layers[ii]*1e-3, TT_layers[0][ii], TT_layers[1][ii], qq_layers[0][ii]*1e3, qq_layers[1][ii]*1e3

    RR = rad_drv.get_R( SST, zz_satur, zz_layers[l_strato], zz_layers, TT_layers, PP_layers, qq_layers, MU_atm, MU_H2O, SOL, FACTOR, ALBEDO )
#    RR = rad_drv.get_R( SST, l_satur, l_strato, zz_layers, TT_layers, PP_layers, qq_layers, MU_atm, MU_H2O, FACTOR_W, FACTOR_C )

    dR_add_w = cloud_forcing_w( aa_w, zz_layers, [0, l_satur[0], l_TI_tp[0], l_strato] )
    dR_add_c = cloud_forcing_c( aa_c, zz_layers, [0, l_satur[1], l_TI_tp[1], l_strato]  )
   
    print "dR_add_w + dR_all", dR_add_w + dR_all
    print "dR_add_c + dR_all", dR_add_c + dR_all

    RR[0] += dR_add_w + dR_all
    RR[1] += dR_add_c + dR_all

    
    return np.r_[ RR[0], RR[1] ]
