# Module
#=============================================================================
from copy import deepcopy
import numpy as np
import datetime
import mks as unit
import util_interp
import errors
#import rad.gasext as gasext
#import rad.emission as emission
#import rad.scatter as scatter
import sys
import os
import commands

from scipy import constants

from params import *

#------------------------------------------------
# options
#------------------------------------------------
# Include molecular absorption ?
MOLABS_ON   = True
# Include molecular absorption ?
CNTNM_ON   = False
# Include Rayleigh scattering ?
RAYLEIGH_ON = True
# Include cloud scattering ?
CLD_ON      = False


ATMFILE_W  = "pass_libradtran/atmfile_w.dat"
ATMFILE_C  = "pass_libradtran/atmfile_c.dat"
INPUTFILE  = "pass_libradtran/tmp.inp"
OUTPUTFILE = "pass_libradtran/tmp.out"


#=============================================================================
# functions
#=============================================================================


#=============================================================================
def make_dic_n_layers( atmfile, zkm_layers ):
    """
    read atmospheric data for O3, O2, CO2, NO2 from afglt.dat
    """
    molname0 = ["alt", "pres", "temp", "air", "O3", "O2", "H2O", "CO2", "NO2"]
    dict_n_layers = {}
    data0 = np.loadtxt( atmfile ).T
    for ii in np.array([4,5,7,8]):
#        func_interp = util_interp.interp_1d( data0[0], data0[ii]/data0[3]*1e6, logx=False )
        func_interp = util_interp.interp_1d( data0[0], data0[ii], logx=False )
        dict_n_layers[molname0[ii]] = func_interp( zkm_layers )
    return data0, dict_n_layers


#=============================================================================
def write_atmfile( zkm_layers, t_layers, p_layers, q_layers, mu_atm, mu_H2O, dict_n_layers, atmprof0, tsurf, atmfile ):
    """
    Set up profile of atmosphere
    """

    n0_layers = p_layers / ( RR * t_layers ) * unit.NA * 1e-6    # m^-3 => cm^-3
    nH2O_layers = q2n( t_layers, p_layers, q_layers, mu_atm, mu_H2O ) # approximating mixingration as specific humidity

#    dict_n_layers = {}
#    for key in dict_ppm_layers :
#        dict_n_layers[key] = n0_layers * dict_ppm_layers[key] * 1e-6

    with open( atmfile, 'w') as f :

        f.write("#     z(km)      p(mb)        T(K)    air(cm-3)    o3(cm-3)     o2(cm-3)    h2o(cm-3)    co2(cm-3)     no2(cm-3)\n")
#        for zi in np.where( atmprof0[0] > zkm_layers[-1] )[0][10:] :
#            for ii in xrange( len( atmprof0 ) ):
#                f.write( "\t" + str(atmprof0[ii][zi]) )
#            f.write( "\n" )
        ScaleHeight_km = ( unit.RR * t_layers[-1] ) / ( mu_atm * gg ) * 1e-3

#        print "np.where( atmprof0[0] > zkm_layers[-1] )[0]", np.where( atmprof0[0] > zkm_layers[-1] )[0][10:]

        for zi in np.where( atmprof0[0] > zkm_layers[-1] )[0][10:] :
            f.write( "\t" + str( atmprof0[0][zi] ) )
            f.write( "\t" + str( p_layers[-1]*1e-2*np.exp( - ( atmprof0[0][zi] - zkm_layers[-1] ) / ScaleHeight_km ) ) )        # pressure in mbar !!
            f.write( "\t" + str( t_layers[-1] ) )             # temperature
            f.write( "\t" + str( n0_layers[-1]*np.exp( - ( atmprof0[0][zi] - zkm_layers[-1] ) / ScaleHeight_km ) ) )             # number density
            f.write( "\t" + str( atmprof0[4][zi] ) )  # O3
            f.write( "\t" + str( atmprof0[5][zi] ) )  # O2
            f.write( "\t" + str( nH2O_layers[-1]*np.exp( - ( atmprof0[0][zi] - zkm_layers[-1] ) / ScaleHeight_km ) ) )          # H2O
            f.write( "\t" + str( atmprof0[7][zi] ) ) # CO2
            f.write( "\t" + str( atmprof0[8][zi] ) ) # NO2
            f.write( "\n" )
        dzi = 1
        if ( len(zkm_layers) > 50 ): 
            dzi = len(zkm_layers)/25
        for zi in np.arange(1,len(zkm_layers))[::-1*dzi] :
            f.write( "\t" + str( zkm_layers[zi] ) )        # altitude in km !!
            f.write( "\t" + str( p_layers[zi]*1e-2 ) )        # pressure in mbar !!
            f.write( "\t" + str( t_layers[zi] ) )             # temperature
            f.write( "\t" + str( n0_layers[zi] ) )            # number density
            f.write( "\t" + str( dict_n_layers['O3'][zi] ) )  # O3
            f.write( "\t" + str( dict_n_layers['O2'][zi] ) )  # O2
            f.write( "\t" + str( nH2O_layers[zi] ) )          # H2O
            f.write( "\t" + str( dict_n_layers['CO2'][zi] ) ) # CO2
            f.write( "\t" + str( dict_n_layers['NO2'][zi] ) ) # NO2
            f.write( "\n" )

        # zero altitude
        f.write( "\t" + "0.0" )        # altitude in km !!
        f.write( "\t" + str( P_surf*1e-2 ) )        # pressure in mbar !!
        f.write( "\t" + str( tsurf ) )             # temperature
        f.write( "\t" + str( P_surf / ( tsurf * RR ) * unit.NA * 1e-6 ) )            # number density
        f.write( "\t" + str( dict_n_layers['O3'][0] ) )  # O3
        f.write( "\t" + str( dict_n_layers['O2'][0] ) )  # O2
        f.write( "\t" + str( nH2O_layers[0] ) )          # H2O
        f.write( "\t" + str( dict_n_layers['CO2'][0] ) ) # CO2
        f.write( "\t" + str( dict_n_layers['NO2'][0] ) ) # NO2
        f.write( "\n" )


#=============================================================================
def q2n( temp, pres, qq, mu_atm, mu_H2O ):
    mixing_ratio = qq * mu_atm / mu_H2O
    number_density = pres / ( unit.RR * temp ) * unit.NA * 1e-6   # m^-3 => cm^-3
    return number_density * mixing_ratio



#=============================================================================
def write_inputfile_mir( list_z, tsurf, atmfile ):

    with open( INPUTFILE, 'w') as f :
        f.write( "data_files_path /Users/yuka/libRadtran-1.7/data/" )
        f.write( "\natmosphere_file " + atmfile                     )
        f.write( "\nsolar_file /Users/yuka/libRadtran-1.7/examples/UVSPEC_LOWTRAN_THERMAL.TRANS" )
        f.write( "\nsurface_temperature " + str(tsurf) )
        f.write( "\nsource thermal" )
        f.write( "\nrte_solver disort2" )
        f.write( "\ntransmittance_wl_file /Users/yuka/libRadtran-1.7/examples/UVSPEC_LOWTRAN_THERMAL.TRANS" )
        f.write( "\ncorrelated_k LOWTRAN" )
#        f.write( "\nzout 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 15.0 20.0 25.0 30.0" )
        f.write( "\nzout" )
        for ii in xrange( len(list_z) ):
            f.write( " " + str(list_z[ii]) )
        f.write( "\noutput_user zout enet" )
        f.write( "\noutput integrate"                           )
        f.write( "\nquiet"                                      )

 

#=============================================================================
def write_inputfile_vis( list_z, atmfile ):

    with open( INPUTFILE, 'w') as f :
        f.write( "data_files_path /Users/yuka/libRadtran-1.7/data/" )
        f.write( "\natmosphere_file " + atmfile                     )
        f.write( "\nalbedo 0.2"                                   )
        f.write( "\nsza 0.0"                                     )
        f.write( "\nrte_solver disort2"                           )
        f.write( "\ncorrelated_k KATO2"                           )
#        f.write( "\nzout 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 15.0 20.0 25.0 30.0" )
        f.write( "\nzout" )
        for ii in xrange( len(list_z) ):
            f.write( " " + str(list_z[ii]) )
        f.write( "\noutput_user zout enet" )
        f.write( "\noutput sum"                                   )
        f.write( "\nquiet"                                        )






#=============================================================================
def read_outputfile( ):

    data = np.loadtxt( OUTPUTFILE ).T[1]
    return data


#=============================================================================
def calc_rad( type, list_l0, z_layers, t_layers, p_layers, q_layers, mu_atm, mu_H2O, tsurf, pool ):


    if ( pool == 0 ):
        atmfile = ATMFILE_W
    else:
        atmfile = ATMFILE_C

    #------------------------------------------------
    # set up atmospheric data profile
    #------------------------------------------------
    zkm_layers = z_layers * 1e-3

    #------------------------------------------------
    # set up atmospheric data profile
    #------------------------------------------------
    atmprof0, dict_n_layers = make_dic_n_layers( "/Users/yuka/libRadtran-1.7/data/atmmod/afglt.dat", zkm_layers )


    write_atmfile( zkm_layers, t_layers, p_layers, q_layers, mu_atm, mu_H2O, dict_n_layers, atmprof0, tsurf, atmfile )


    #------------------------------------------------
    # set input file to be passed to libRadtran
    #------------------------------------------------
    list_z = []
    for ii in xrange( len(list_l0) ):
        list_z.append( zkm_layers[list_l0[ii]] ) # km
    if ( type == "mir" ):
        write_inputfile_mir( list_z, tsurf, atmfile )
    elif ( type == "vis" ):
        write_inputfile_vis( list_z, atmfile )
    else :
        errors.exit_msg("Invalid radiation type.")

    #------------------------------------------------
    # execute libRadtran
    #------------------------------------------------
    check = os.system( 'uvspec < ' + INPUTFILE + " > " + OUTPUTFILE )
    if check != 0 :
        errors.exit_msg("Something is wrong with uvspec.")        

    #------------------------------------------------
    # read the output
    #------------------------------------------------
    array_rad = read_outputfile( )
    return array_rad



#=============================================================================
def get_R( sst_w, sst_c, qq_TI_c, aa_w, factor ):
#
#    Rw_0 = griddata( pnts, vals[0], np.array( [ sst_w, sst_c ] ), method="cubic" )[0]
#    Rw_B = griddata( pnts, vals[1], np.array( [ sst_w, sst_c ] ), method="cubic" )[0]
#    Rw_T = griddata( pnts, vals[2], np.array( [ sst_w, sst_c ] ), method="cubic" )[0]
#    Rc_0 = griddata( pnts, vals[3], np.array( [ sst_w, sst_c ] ), method="cubic" )[0]
#    Rc_B = griddata( pnts, vals[4], np.array( [ sst_w, sst_c ] ), method="cubic" )[0]
#    Rc_T = griddata( pnts, vals[5], np.array( [ sst_w, sst_c ] ), method="cubic" )[0]
#    return np.array( [[ Rw_0, Rw_B, Rw_T ], [Rc_0, Rc_B, Rc_T]] )

    #----------------------------------------------------
    # step 2: location of the trade inversion
    #----------------------------------------------------
    l_TI_btm = np.zeros(2)
    l_TI_tp  = np.zeros(2)

    l_TI_btm[0] = np.where( zz_layers > zz_TI_btm[0] )[0][0]
    l_TI_btm[1] = np.where( zz_layers > zz_TI_btm[1] )[0][0]
    l_TI_tp[0]  = np.where( zz_layers > zz_TI_tp[0] )[0][0]
    l_TI_tp[1]  = np.where( zz_layers > zz_TI_tp[1] )[0][0]

    #----------------------------------------------------
    # step 3: find q_B
    #----------------------------------------------------
    qq_layers.T[0] = rh2q( TT_layers.T[0], P_surf, RH_surf )
    ss.T[0]        = find_s( TT_layers.T[0], 0. )

    #----------------------------------------------------
    # step 4: find saturation level
    #----------------------------------------------------
    for pool in xrange(2):
        zz_satur[pool] = find_saturation_level( TT_layers[pool][0], 
                                                qq_layers[pool][0] )
        l_satur[pool]  = np.where( zz_layers > zz_satur[pool] )[0][0]
        
    #----------------------------------------------------
    # step 5: find T profile below the saturation level
    #----------------------------------------------------
    for pool in xrange(2):
        TT_layers[pool][:l_satur[pool]+1] = find_Tprof_dryadiabat( TT_layers[pool][0], zz_layers[:l_satur[pool]+1] )
        qq_layers[pool] = np.tile( qq_layers[pool][0], N_LAYER )

    #----------------------------------------------------
    # step 6: find P profile below the saturation level
    #----------------------------------------------------
    pool = 0 # warm pool
    PP_layers[pool][:l_satur[pool]+1] = find_Pprof( zz_layers[:l_satur[pool]+1], 
                                                    TT_layers[pool][:l_satur[pool]+1], 
                                                    P_surf )
    
    #----------------------------------------------------
    # step 7a: find T profile for warm pool above the saturation level (moist adiabat)
    #----------------------------------------------------
    pool = 0 # warm pool
    # approximation HERE!!!!!!!!!
    TT_layers[pool][l_satur[pool]:], PP_layers[pool][l_satur[pool]:], l_strato = find_Tprof_moistadiabat( TT_layers[pool][l_satur[pool]], 
                                                                                                          PP_layers[pool][l_satur[pool]], 
                                                                                                          zz_layers[l_satur[pool]:] )

    # correction
    l_strato = l_strato + l_satur[pool]

#    TT_layers[pool][l_satur[pool]:], PP_layers[pool][l_satur[pool]:] = find_Tprof_moistadiabat( TT_layers[pool][l_satur[pool]], 
#                                                                                                PP_layers[pool][l_satur[pool]], 
#                                                                                                zz_layers[l_satur[pool]:] ) 
    qq_layers[pool][l_satur[pool]:] = rh2q( TT_layers[pool][l_satur[pool]:], 
                                            PP_layers[pool][l_satur[pool]:], 
                                            np.ones( N_LAYER - l_satur[pool] ) ) 

#    #----------------------------------------------------
#    # step 8: find tropopause
#    #----------------------------------------------------
#    pool = 0
#    if ( TT_layers[0][-1] > T_strato ) :
#        if ( verbose == True ):
#            with open( OUTPUTFILE, 'a') as f:
#                f.write( "# invalid ! TT_layers[0][-1] > T_strato " )
#        return np.array([ 0.0, 0.0 ])
#
#    l_strato, TT_layers[pool] = find_tropopause( TT_layers[pool] )


    #----------------------------------------------------
    # step 9: specific humidity above tropopause
    #----------------------------------------------------
    qT = rh2q( TT_layers[pool][l_strato], PP_layers[pool][l_strato], 1.0 )
    # specific humidity above tropopause
    for pool in xrange(2):
        qq_layers[pool][l_strato:] = np.tile( qT, N_LAYER - l_strato )


    #----------------------------------------------------
    # step 7b: above the trade inversion, the temperature profile of cold pool is same as that of warm pool
    #----------------------------------------------------
    TT_layers[1][l_TI_tp[1]:] = TT_layers[0][l_TI_tp[1]:]


    # just for now
    PP_layers[1] = PP_layers[0]


    #----------------------------------------------------
    # step 11: find T profile between saturation level and inversion with mixing theory
    #----------------------------------------------------
    for pool in xrange(2):

        if ( pool == 0 ) : # warm pool
            qq_layers[pool][l_TI_tp[pool]] =  rh2q( TT_layers[pool][l_TI_tp[pool]], 
                                                    PP_layers[pool][l_TI_tp[pool]], 
                                                    RH_surf )

        elif ( pool == 1 ) : # cold pool
            qq_layers[pool][l_TI_tp[pool]] = qq_TI_c

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
    # step 12: determine temperature in TI simply by linear interpolation
    #----------------------------------------------------
    for pool in xrange(2):

        TT_layers[pool][l_TI_btm[pool]:l_TI_tp[pool]+1] = connect_linear( zz_layers[l_TI_btm[pool]:l_TI_tp[pool]+1], 
                                                                          TT_layers[pool][l_TI_btm[pool]], 
                                                                          TT_layers[pool][l_TI_tp[pool]] )
        qq_layers[pool][l_TI_btm[pool]:l_TI_tp[pool]+1] = connect_linear( zz_layers[l_TI_btm[pool]:l_TI_tp[pool]+1], 
                                                                          qq_layers[pool][l_TI_btm[pool]], 
                                                                          qq_layers[pool][l_TI_tp[pool]] )

    


    #----------------------------------------------------
    # step 10: find moisture profile above the TI
    #----------------------------------------------------
    for pool in xrange(2):

        if ( pool == 0 ) : # warm pool
            RH_TI2 = RH_surf
        elif ( pool == 1 ) : # cold pool
            RH_TI2 = q2rh( TT_layers[pool][l_TI_tp[pool]], 
                           PP_layers[pool][l_TI_tp[pool]], 
                           qq_layers[pool][l_TI_tp[pool]] )
            
        RH_layers[pool][l_TI_tp[pool]:] = find_RHprof_linear_with_P( RH_TI2, 
                                                                     TT_layers[pool][l_TI_tp[pool]], 
                                                                     zz_layers[l_TI_tp[pool]:], 
                                                                     ALPHA )

        qq_layers[pool][l_TI_tp[pool]:] = rh2q( TT_layers[pool][l_TI_tp[pool]:],
                                                PP_layers[pool][l_TI_tp[pool]:], 
                                                RH_layers[pool][l_TI_tp[pool]:] )



    if ( verbose == True ):
        for ii in xrange( N_LAYER ):
            print zz_layers[ii]*1e-3, TT_layers[0][ii], TT_layers[1][ii], qq_layers[0][ii]*1e3, qq_layers[1][ii]*1e3

    #----------------------------------------------------
    # step 13: compute radiation
    #----------------------------------------------------
    print "in the warm pool"
    rad_vis_w = factor*calc_rad( "vis", [0, l_satur[0], l_strato], zz_layers, TT_layers[0], PP_layers[0], qq_layers[0], mu_atm, mu_H2O, sst_w, 0  )
    print "rad_vis", rad_vis_w
    rad_mir_w = calc_rad( "mir", [0, l_satur[0], l_strato], zz_layers, TT_layers[0], PP_layers[0], qq_layers[0], mu_atm, mu_H2O, sst_w, 0  )
    print "rad_mir", rad_mir_w
    print "NET: ", rad_mir_w + rad_vis_w
    
    print ""
    
    print "in the cold pool"
    rad_vis_c = factor*calc_rad( "vis", [0, l_satur[1], l_strato], zz_layers, TT_layers[1], PP_layers[1], qq_layers[1], mu_atm, mu_H2O, sst_c, 1 )
    print "rad_vis", rad_vis_c
    rad_mir_c = calc_rad( "mir", [0, l_satur[1], l_strato], zz_layers, TT_layers[1], PP_layers[1], qq_layers[1], mu_atm, mu_H2O, sst_c, 1 )
    print "rad_mir", rad_mir_c
    print "NET: ", rad_mir_c + rad_vis_c


    R_w = rad_mir_w + rad_vis_w
    R_c = rad_mir_c + rad_vis_c

    print ""
    print "np.array( [R_w, R_c] )", np.array( [R_w, R_c] )
    print ""
    return np.array( [R_w, R_c] )


