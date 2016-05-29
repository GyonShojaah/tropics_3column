# Module
#=============================================================================
from copy import deepcopy
import numpy as np
import datetime
import mks as unit
import util_interp
import errors
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


ATMFILE_W  = "pass_libradtran/atmfile_q5_w.dat"
ATMFILE_C  = "pass_libradtran/atmfile_q5_c.dat"
INPUTFILE_MIR  = "pass_libradtran/mir.inp"
INPUTFILE_VIS  = "pass_libradtran/vis.inp"
OUTPUTFILE_MIR = "pass_libradtran/mir.out"
OUTPUTFILE_VIS = "pass_libradtran/vis.out"


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
        if ( len(zkm_layers) > 30 ): 
            dzi = len(zkm_layers)/50
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

    with open( INPUTFILE_MIR, 'w') as f :
        f.write( "data_files_path /Users/yuka/libRadtran-1.7/data/" )
        f.write( "\natmosphere_file " + atmfile                     )
        f.write( "\nsolar_file /Users/yuka/libRadtran-1.7/examples/UVSPEC_LOWTRAN_THERMAL.TRANS" )
        f.write( "\nsurface_temperature " + str(tsurf) )
        f.write( "\nsource thermal" )
        f.write( "\nrte_solver disort2" )
        f.write( "\ntransmittance_wl_file /Users/yuka/libRadtran-1.7/examples/UVSPEC_LOWTRAN_THERMAL.TRANS" )
        f.write( "\ncorrelated_k LOWTRAN" )
###        f.write( "\nzout 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 15.0 20.0 25.0 30.0" )
        f.write( "\nzout" )
        for ii in xrange( len(list_z) ):
            f.write( " " + str(list_z[ii]) )
        f.write( "\noutput_user zout enet" )
        f.write( "\noutput integrate"                           )
        f.write( "\nquiet"                                      )

 

#=============================================================================
def write_inputfile_vis( list_z, atmfile, albd ):

    with open( INPUTFILE_VIS, 'w') as f :
        f.write( "data_files_path /Users/yuka/libRadtran-1.7/data/" )
        f.write( "\natmosphere_file " + atmfile                     )
        f.write( "\nalbedo " + str(albd)                            )
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
def read_outputfile( outputfile ):

    zout = np.loadtxt( outputfile ).T[0]
    data = np.loadtxt( outputfile ).T[1]
    return zout, data


#=============================================================================
def calc_rad( type, list_l0, z_layers, t_layers, p_layers, q_layers, mu_atm, mu_H2O, tsurf, pool, factor, albedo ):


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
        input  = INPUTFILE_MIR
        output = OUTPUTFILE_MIR
    elif ( type == "vis" ):
        write_inputfile_vis( list_z, atmfile, albedo )
        input  = INPUTFILE_VIS
        output = OUTPUTFILE_VIS
    else :
        errors.exit_msg("Invalid radiation type.")


    check = os.system( 'uvspec < ' + input + " > " + output )
    if check != 0 :
        errors.exit_msg("Something is wrong with uvspec.")


    #------------------------------------------------
    # read the output
    #------------------------------------------------
    array_zout, array_rad = read_outputfile( output )


    print "array_zout", array_zout
    print "factor", factor
    print "array_rad",  array_rad
    return array_zout, factor*array_rad



#=============================================================================
def get_R( sst, l_satur, l_strato, z_layers, T_layers, P_layers, q_layers, mu_atm, mu_H2O, sol, factor, albedo ):


    #----------------------------------------------------
    # step 13: compute radiation
    #----------------------------------------------------
    print "in the warm pool"
    zout, rad_vis_w = calc_rad( "vis", [0, l_satur[0], l_strato], z_layers, T_layers[0], P_layers[0], q_layers[0], mu_atm, mu_H2O, sst[0], 0, factor[0], albedo  )
    print "rad_vis", rad_vis_w
    zout, rad_mir_w = calc_rad( "mir", [0, l_satur[0], l_strato], z_layers, T_layers[0], P_layers[0], q_layers[0], mu_atm, mu_H2O, sst[0], 0, 1.0, albedo  )
    print "rad_mir", rad_mir_w
    print "NET: ", rad_mir_w + rad_vis_w
    
    print ""
    
    print "in the cold pool"
    zout, rad_vis_c = calc_rad( "vis", [0, l_satur[1], l_strato], z_layers, T_layers[1], P_layers[1], q_layers[1], mu_atm, mu_H2O, sst[1], 1, factor[1], albedo )
    print "rad_vis", rad_vis_c
    zout, rad_mir_c = calc_rad( "mir", [0, l_satur[1], l_strato], z_layers, T_layers[1], P_layers[1], q_layers[1], mu_atm, mu_H2O, sst[1], 1, 1.0, albedo )
    print "rad_mir", rad_mir_c
    print "NET: ", rad_mir_c + rad_vis_c


    R_w = rad_mir_w + rad_vis_w
    R_c = rad_mir_c + rad_vis_c

#    print "-----------------------------------------------------------------------------------"
#    print "z[km] \t warm, vis \t warm, mir \t cold, vis \t cold, mir \n"
#    for zi in xrange( len(zout) ):
#        print zout[zi], rad_vis_w[zi], rad_mir_w[zi], rad_vis_c[zi], rad_mir_c[zi]
#    print "-----------------------------------------------------------------------------------"

    return np.array( [R_w, R_c] )


