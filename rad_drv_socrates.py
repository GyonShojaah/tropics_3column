# Module
#=============================================================================

import netCDF4
import numpy as np
import util_interp
import errors
import sys
import subprocess
from pass_socrates.read_write_socrates_netcdf import *

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

SOCRATES_DIR = ""
BASE_NAME = [ "atm_w" , "atm_c" ]
ALBEDO_LW = 0.05
ALBEDO_SW = 0.05

DIM_LON = 1
DIM_LAT = 1

molmassratio = { 
    'o2'  : 0.231400,
    'co2' : 0.539380e-03, 
    'n2o' : 0.0, 
    'ch4' : 0.0 
    }


plev0_layers =  np.array( [ .502500E+01,  .120040E+02,  .168350E+02,  .236250E+02,
                     .331340E+02,  .465250E+02,  .652160E+02,  .915110E+02,
                     .128300E+03,  .180090E+03,  .252470E+03,  .354200E+03,
                     .496860E+03,  .696480E+03,  .978940E+03,  .136970E+04,
                     .192370E+04,  .269880E+04,  .379010E+04,  .531150E+04,
                     .738890E+04,  .976910E+04,  .121390E+05,  .145290E+05,
                      .169040E+05,  .192760E+05,  .216550E+05,  .240380E+05,
                     .264050E+05,  .287870E+05,  .311610E+05,  .335390E+05,
                     .359160E+05,  .382950E+05,  .406720E+05,  .430520E+05,
                     .454270E+05,  .478090E+05,  .501820E+05,  .525600E+05,
                     .549440E+05,  .573150E+05,  .596920E+05,  .620740E+05,
                     .644470E+05,  .668250E+05,  .692030E+05,  .715820E+05,
                     .739580E+05,  .763350E+05,  .787130E+05,  .810960E+05,
                     .834680E+05,  .858460E+05,  .882230E+05,  .906030E+05,
                     .929780E+05,  .953560E+05,  .977330E+05,  .100110E+06 ] )
o30_layers = np.array( [ .611567E-06,  .114955E-05,  .158121E-05,  .195664E-05,
                          .244117E-05,  .285659E-05,  .361983E-05,  .457777E-05,
                          .584955E-05,  .771771E-05,  .102454E-04,  .129101E-04,
                          .145514E-04,  .144036E-04,  .133629E-04,  .116207E-04,
                          .993219E-05,  .808276E-05,  .567424E-05,  .385413E-05,
                          .218847E-05,  .118164E-05,  .899844E-06,  .743435E-06,
                          .589429E-06,  .433302E-06,  .354990E-06,  .302726E-06,
                          .251060E-06,  .213030E-06,  .192483E-06,  .178017E-06,
                          .159612E-06,  .147516E-06,  .135683E-06,  .126216E-06,
                          .116192E-06,  .107514E-06,  .103504E-06,  .972589E-07,
                          .923623E-07,  .893382E-07,  .848012E-07,  .812518E-07,
                          .786883E-07,  .753030E-07,  .721148E-07,  .693525E-07,
                          .666582E-07,  .641278E-07,  .616621E-07,  .613340E-07,
                          .599852E-07,  .582105E-07,  .567639E-07,  .556471E-07,
                          .545617E-07,  .531814E-07,  .518557E-07,  .505964E-07 ] )
                          
l2_satur = [0,0]


#=============================================================================
# functions
#=============================================================================

#=============================================================================



#=============================================================================
def read_outputfile( basename, list_l0 ):
 
#    print 'Reading outputfile ' + basename + '...' 
#    print 'at ', list_l0
    ncfile_r = netCDF4.Dataset( SOCRATES_DIR+basename+'.nflx', 'r', format='NETCDF3_CLASSIC' )
    nflx_r = ncfile_r.variables['nflx']
    data = nflx_r[list_l0,0,0]
#    print 'data', data
    return data


#=============================================================================
def calc_rad_mir( pool, list_l0 ):

    # $RAD_SCRIPT/Cl_run_cdf -B 7460_28 -s $RAD_DIR/data/spectra/hadgem1/sp_lw_hadgem1_3 -R 1 9 -I -G 5 0 -t 12 +R -v 11 -g 2 -c -C 2 -K 1 -d 1 -i 1    
    basis = np.array([0.])
    write_latlon( SOCRATES_DIR, ALBEDO_LW, BASE_NAME[pool], 'surf', 'alb', 'None', 'Albedo weights', basis=basis)

#    print ''
    cmd = 'Cl_run_cdf -s $RAD_DATA/spectra/ga7/sp_lw_ga7 -R 1 9 -I -g 2 -c -t 12 -v 13 -C 5 -B '+SOCRATES_DIR+BASE_NAME[pool]
#    print cmd
    subprocess.call( cmd, shell=True )
#    if check != 0 :
#        errors.exit_msg("Something is wrong with Cl_run_cdf.")

    cmd = 'mv '+BASE_NAME[pool]+'.nflx '+BASE_NAME[pool]+'_mir.nflx '
    subprocess.call( cmd, shell=True )
    cmd = 'mv '+BASE_NAME[pool]+'.uflx '+BASE_NAME[pool]+'_mir.uflx '
    subprocess.call( cmd, shell=True )
    cmd = 'mv '+BASE_NAME[pool]+'.dflx '+BASE_NAME[pool]+'_mir.dflx '
    subprocess.call( cmd, shell=True )
    array_rad = read_outputfile( BASE_NAME[pool]+'_mir', list_l0 )
#    print 'end THERMAL'

    return array_rad


#=============================================================================
def calc_rad_sol( pool, list_l0 ):
    
    basis = np.array([0.])
    write_latlon( SOCRATES_DIR, ALBEDO_SW, BASE_NAME[pool], 'surf', 'alb', 'None', 'Albedo weights', basis=basis)
    cmd = 'Cl_run_cdf -s $RAD_DATA/spectra/ga7/sp_sw_ga7  -R 1 6 -S -g 2 -c -t 16 -r -v 13 -C 5 -B '+SOCRATES_DIR+BASE_NAME[pool]
#    print cmd
    check = subprocess.call( cmd, shell=True )
#    check = os.system( 'Cl_run_cdf -s $RAD_DATA/spectra/ga7/sp_sw_ga7  -R 1 6 -S -g 2 -c -t 16 -r -v 13 -C 5 -B '+BASE_NAME[pool] )
#    if check != 0 :
#        errors.exit_msg("Something is wrong with Cl_run_cdf.")
    cmd = 'mv '+BASE_NAME[pool]+'.nflx '+BASE_NAME[pool]+'_sol.nflx '
    subprocess.call( cmd, shell=True )
    cmd = 'mv '+BASE_NAME[pool]+'.uflx '+BASE_NAME[pool]+'_sol.uflx '
    subprocess.call( cmd, shell=True )
    cmd = 'mv '+BASE_NAME[pool]+'.dflx '+BASE_NAME[pool]+'_sol.dflx '
    subprocess.call( cmd, shell=True )
    cmd = 'mv '+BASE_NAME[pool]+'.sflx '+BASE_NAME[pool]+'_sol.sflx '
    subprocess.call( cmd, shell=True )
    cmd = 'mv '+BASE_NAME[pool]+'.vflx '+BASE_NAME[pool]+'_sol.vflx '
    subprocess.call( cmd, shell=True )
    array_rad = read_outputfile( BASE_NAME[pool]+'_sol', list_l0 )
#    print 'end SOLAR'
#    print ''
    return array_rad



#=============================================================================
def get_R( sst, z_satur, z_strato, z_layers_in, T_layers_in, P_layers_in, q_layers_in, mu_atm, mu_H2O, sol, factor, albedo ):


    # from top to bottom

    z_layers = z_layers_in[::-1]
    T_layers = np.zeros_like(T_layers_in)
    P_layers = np.zeros_like(P_layers_in)
    q_layers = np.zeros_like(q_layers_in)

    for pool in xrange(2) :
        T_layers[pool] = T_layers_in[pool][::-1]
        P_layers[pool] = P_layers_in[pool][::-1]
        q_layers[pool] = q_layers_in[pool][::-1]


    l_strato = np.argmin( abs( z_layers[1::2] - z_strato ) )
    l_satur  = [0, 0]
    for pool in xrange(2):
        l_satur[pool]  = np.argmin( abs( z_layers[1::2] - z_satur[pool] ) )

#    print ""
#    print "Writing atmospheric files..."
    for pool in xrange(2) :

        # write_latlon(path, x, prefix, suffix, name, unit, title,
        #     p_lev=None, basis=None):

        # solar irradiance
        write_latlon( SOCRATES_DIR, sol*factor[pool], BASE_NAME[pool], 'stoa', 
                      'stoa', 'WM-2', 'Solar Irradiance' )
        
        # solar zenith angle
        write_latlon( SOCRATES_DIR, 0., BASE_NAME[pool], 'szen', 
                      'szen', 'degree', 'Solar zenith angle' ) # ??


        # surface
#        write_latlon( SOCRATES_DIR, albedo, BASE_NAME[pool], 'surf', 
#                      'surf', 'no unit', 'surf', basis=basis ) # ??
#
#        write_latlon( SOCRATES_DIR, ALBEDO_SW, BASE_NAME[pool], 'surfsw', 'alb', 'None',
#                     'Albedo weights', basis=basis)

        # temperature
        plev_layers = P_layers[pool][1::2]
        write_latlon( SOCRATES_DIR, T_layers[pool][1::2], BASE_NAME[pool], 't', 
                      't', 'K', 'Temperature', p_lev=plev_layers ) # ??

        plev_boundary = P_layers[pool][0::2]
        write_latlon( SOCRATES_DIR, T_layers[pool][0::2], BASE_NAME[pool], 'tl', 
                      'tl', 'K', 'Temperature on levels', p_lev=plev_boundary ) # ??

        plev_surf = np.array( [ P_layers[pool][-1] ] )
        write_latlon( SOCRATES_DIR, sst[pool], BASE_NAME[pool], 'tstar', 
                      'tstar', 'K', 'Surface temperature', p_lev=plev_surf ) # ??

        # specific humidity
        write_latlon( SOCRATES_DIR, q_layers[pool][1::2], BASE_NAME[pool], 'q', 
                      'q', 'g/g', 'Specific humidity', p_lev=plev_layers ) # ??


        # gases
        unit_layers = np.ones( len( plev_layers ) )
        for key in molmassratio :
            gas_layers = np.ones( len( plev_layers ) )*molmassratio[key]
            write_latlon( SOCRATES_DIR, gas_layers, BASE_NAME[pool], key, 
                          key, '', 'gas mixing ratio', p_lev=plev_layers ) 

        # O3
        func_interp = util_interp.interp_1d( plev0_layers, o30_layers, bounds_error=False, fill_value=0. )
        O3_layers   = func_interp( plev_layers )
        write_latlon( SOCRATES_DIR, O3_layers, BASE_NAME[pool], 'o3', 
                      'o3', '', 'o3 mixing ratio', p_lev=plev_layers ) 


    l_end = len( P_layers[pool][0::2] ) - 1
    l_list = np.array([l_end, l_satur[0], l_strato])

#    print ''
#    print 'l_list, z_list', l_list, z_layers[1::2][l_list]
#    print ""
#    print "Calculating radiation fields in the warm pool..."
    rad_vis_w = calc_rad_sol( 0, l_list  )
#    print "rad_vis", rad_vis_w
    rad_mir_w = calc_rad_mir( 0, l_list )
#    print "rad_mir", rad_mir_w
#    print "NET: ", rad_mir_w + rad_vis_w
    
#    print ""
#    print "Calculating radiation fields in the cold pool..."
    rad_vis_c = calc_rad_sol( 1, l_list )
#    print "rad_vis", rad_vis_c
    rad_mir_c = calc_rad_mir( 1, l_list )
#    print "rad_mir", rad_mir_c
#    print "NET: ", rad_mir_c + rad_vis_c


    R_w = rad_mir_w + rad_vis_w
    R_c = rad_mir_c + rad_vis_c

#    print "-----------------------------------------------------------------------------------"
#    print "z[km] \t warm, vis \t warm, mir \t cold, vis \t cold, mir \n"
#    for zi in xrange( len(zout) ):
#        print zout[zi], rad_vis_w[zi], rad_mir_w[zi], rad_vis_c[zi], rad_mir_c[zi]
#    print "-----------------------------------------------------------------------------------"

    return np.array( [R_w, R_c] )


