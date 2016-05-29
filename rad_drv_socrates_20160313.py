# Module
#=============================================================================

import netCDF4
import numpy as np
import util_interp
import errors
import sys
import os


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

BASE_NAME = [ "pass_socrates/atm_w" , "pass_socrates/atm_c" ]

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
def add_lonlat( file_w ):
    """
    --------------------------------------------------------------------
    called from write_irradiance, write_TPprof, write_gasprof
    --------------------------------------------------------------------
    create demensions for longitude and latitude
    and set these variables
    --------------------------------------------------------------------
    """
    file_w.createDimension('lon', DIM_LON )
    file_w.createDimension('lat', DIM_LAT )
    lon_w = file_w.createVariable('lon','f8',('lon'))
    lat_w = file_w.createVariable('lat','f8',('lat'))
    lon_w[:] = np.zeros( DIM_LON ) # dummy
    lon_w.units = 'degree'
    lon_w.title = 'Longitude'
    lat_w[:] = np.zeros( DIM_LAT ) # dummy
    lat_w.units = 'degree'
    lat_w.title = 'Latitude'
    return file_w


#=============================================================================
def add_plev( file_w, plev_layers ):
    """
    --------------------------------------------------------------------
    called from write_irradiance, write_TPprof, write_gasprof
    --------------------------------------------------------------------
    create demensions for longitude and latitude
    and set these variables
    --------------------------------------------------------------------
    """
    dim_plev  = len( plev_layers )
    file_w.createDimension('plev', dim_plev )
    plev_w    = file_w.createVariable('plev','f8',('plev'))
    plev_w[:] = plev_layers[:]
    return file_w


#=============================================================================
def add_variable2D( file_w, paramname, paramvalue, stype ):
    """
    --------------------------------------------------------------------
    called from write_irradiance, write_TPprof, write_gasprof
    --------------------------------------------------------------------
    create demensions for longitude and latitude
    and set these variables
    --------------------------------------------------------------------
    """
    param_w      = file_w.createVariable( paramname, stype, ('lon','lat') )
    param_w[:,:] = np.ones( [ DIM_LAT, DIM_LON ] ) * paramvalue
    return file_w


#=============================================================================
def add_variable3D( file_w, paramname, paramvalue_layers, stype ):
    """
    --------------------------------------------------------------------
    called from write_irradiance, write_TPprof, write_gasprof
    --------------------------------------------------------------------
    create demensions for longitude and latitude
    and set these variables
    --------------------------------------------------------------------
    """
    param_w      = file_w.createVariable( paramname, stype, ('plev','lon','lat') )
    print 'len(paramvalue_layers)', len(paramvalue_layers)
    param_w[:,0,0] = paramvalue_layers
    return file_w


#=============================================================================
def write_szen( basename ):
    """
    --------------------------------------------------------------------
    write irradiance to .nc file
    --------------------------------------------------------------------
    """
    ncfile_w    = netCDF4.Dataset( basename+'.szen', mode='w', format='NETCDF3_CLASSIC') 
    ncfile_w    = add_lonlat( ncfile_w )
    ncfile_w    = add_variable2D( ncfile_w, 'szen', 0., 'f8' )
    ncfile_w.close()


#=============================================================================
def write_irradiance( basename, sol ):
    """
    --------------------------------------------------------------------
    write irradiance to .nc file
    --------------------------------------------------------------------
    """
    ncfile_w    = netCDF4.Dataset( basename+'.stoa', mode='w', format='NETCDF3_CLASSIC') 
    ncfile_w    = add_lonlat( ncfile_w )
    ncfile_w    = add_variable2D( ncfile_w, 'stoa', sol, 'f8' )
    ncfile_w.close()


#=============================================================================
def write_TPprof( basename, p_layers, t_layers, tsurf ):

    # surface temperature
    ncfile_w    = netCDF4.Dataset( basename+'.tstar',mode='w', format='NETCDF3_CLASSIC') 
    ncfile_w    = add_lonlat( ncfile_w )
    ncfile_w    = add_variable2D( ncfile_w, 'tstar', tsurf, 'f8' )
    ncfile_w.close()

    # vertical temperature profile
    plev_layers = p_layers[1::2]
    print 'len(plev_layers)', len(plev_layers)
    ncfile_w    = netCDF4.Dataset( basename+'.t',mode='w', format='NETCDF3_CLASSIC') 
    ncfile_w    = add_lonlat( ncfile_w )
    ncfile_w    = add_plev( ncfile_w, plev_layers )
    ncfile_w    = add_variable3D( ncfile_w, 't', t_layers[1::2], 'f8' )
    ncfile_w.close()

    # vertical temperature profile at the layer bounderies
    plev_layers = p_layers[0::2]
    ncfile_w    = netCDF4.Dataset( basename+'.tl',mode='w', format='NETCDF3_CLASSIC') 
    ncfile_w    = add_lonlat( ncfile_w )
    ncfile_w    = add_plev( ncfile_w, plev_layers )
    ncfile_w    = add_variable3D( ncfile_w, 'tl', t_layers[0::2], 'f8' )
    ncfile_w.close()


#=============================================================================
def write_H2Oprof( basename, p_layers, q_layers ): 

    # vertical temperature profile
    plev_layers = p_layers[1::2]
    ncfile_w    = netCDF4.Dataset( basename+'.q',mode='w', format='NETCDF3_CLASSIC') 
    ncfile_w    = add_lonlat( ncfile_w )
    ncfile_w    = add_plev( ncfile_w, plev_layers )
    ncfile_w    = add_variable3D( ncfile_w, 't', q_layers[1::2], 'f8' )
    ncfile_w.close()


#=============================================================================
def write_gasprof( basename, p_layers ): 

    plev_layers = p_layers[1::2]
    unit_layers = np.ones( len( plev_layers ) )
    for key in molmassratio :
        ncfile_w    = netCDF4.Dataset( basename+'.'+key,mode='w', format='NETCDF3_CLASSIC') 
        ncfile_w    = add_lonlat( ncfile_w )
        ncfile_w    = add_plev( ncfile_w, plev_layers )
        ncfile_w    = add_variable3D( ncfile_w, key, unit_layers*molmassratio[key] , 'f8' )
        ncfile_w.close()


#=============================================================================
def write_O3prof( basename, p_layers ): 

    plev_layers = p_layers[1::2]
    func_interp = util_interp.interp_1d( plev0_layers, o30_layers, bounds_error=False, fill_value=0. )
    O3_layers   = func_interp( plev_layers )
    ncfile_w    = netCDF4.Dataset( basename+'.o3',mode='w', format='NETCDF3_CLASSIC') 
    ncfile_w    = add_lonlat( ncfile_w )
    ncfile_w    = add_plev( ncfile_w, plev_layers )
    ncfile_w    = add_variable3D( ncfile_w, 'o3', O3_layers, 'f8' )
    ncfile_w.close()


#=============================================================================
def write_surf( basename, albd ):
    """
    --------------------------------------------------------------------
    write irradiance to .nc file
    --------------------------------------------------------------------
    """
    ncfile_w    = netCDF4.Dataset( basename+'.surf', mode='w', format='NETCDF3_CLASSIC') 
    ncfile_w    = add_lonlat( ncfile_w )
    ncfile_w.createDimension('basis', 1 )
    basis_w = ncfile_w.createVariable('basis','f8',('basis'))
    basis_w[:] = np.zeros( 1 ) # dummy
    albd_w      = ncfile_w.createVariable( 'alb', 'f8', ('basis','lon','lat') )
    albd_w[0,0,0] = albd
    ncfile_w.close()

    ncfile_w    = netCDF4.Dataset( basename+'.surf_lw', mode='w', format='NETCDF3_CLASSIC') 
    ncfile_w    = add_lonlat( ncfile_w )
    ncfile_w.createDimension('basis', 1 )
    basis_w = ncfile_w.createVariable('basis','f8',('basis'))
    basis_w[:] = np.zeros( 1 ) # dummy
    albd_w      = ncfile_w.createVariable( 'alb', 'f8', ('basis','lon','lat') )
    albd_w[0,0,0] = albd
    ncfile_w.close()

    ncfile_w    = netCDF4.Dataset( basename+'.surf_sw', mode='w', format='NETCDF3_CLASSIC') 
    ncfile_w    = add_lonlat( ncfile_w )
    ncfile_w.createDimension('basis', 1 )
    basis_w = ncfile_w.createVariable('basis','f8',('basis'))
    basis_w[:] = np.zeros( 1 ) # dummy
    albd_w      = ncfile_w.createVariable( 'alb', 'f8', ('basis','lon','lat') )
    albd_w[0,0,0] = albd
    ncfile_w.close()




#=============================================================================
def read_outputfile( basename, list_l0 ):

    ncfile_r = netCDF4.Dataset( basename+'.nflx', 'r', format='NETCDF3_CLASSIC' )
    nflx_r = ncfile_r.variables['nflx']
    data = nflx_r[ list_l0 ] 
    return data


#=============================================================================
def calc_rad_mir( pool, list_l0 ):
    
    print 'Cl_run_cdf -s $RAD_DATA/spectra/dev/sp_lw_9_jm3 -R 1 9 -I -g 2 -c -t 12 -v 13 -C 5 -B ' + BASE_NAME[pool] 
    check = os.system( 'Cl_run_cdf -s $RAD_DATA/spectra/dev/sp_lw_9_jm3 -R 1 9 -I -g 2 -c -t 12 -v 13 -C 5 -B '+BASE_NAME[pool] )
    if check != 0 :
        errors.exit_msg("Something is wrong with Cl_run_cdf.")

    array_rad = read_outputfile( BASE_NAME[pool], list_l0 )
    return array_rad


#=============================================================================
def calc_rad_sol( pool, list_l0 ):
    
    print 'Cl_run_cdf -s $RAD_DATA/spectra/dev/sp_sw_6_jm2  -R 1 6 -S -g 2 -c -t 16 -v 13 -C 5 -B '+BASE_NAME[pool]
    check = os.system( 'Cl_run_cdf -s $RAD_DATA/spectra/dev/sp_sw_6_jm2  -R 1 6 -S -g 2 -c -t 16 -v 13 -C 5 -B '+BASE_NAME[pool] )
    if check != 0 :
        errors.exit_msg("Something is wrong with Cl_run_cdf.")

    array_rad = read_outputfile( BASE_NAME[pool], list_l0 )
    return array_rad



#=============================================================================
def get_R( sst, l_satur_in, l_strato_in, z_layers_in, T_layers_in, P_layers_in, q_layers_in, mu_atm, mu_H2O, sol, factor, albedo ):

    z_layers = z_layers_in[::-1]
    T_layers = T_layers_in[::-1]
    P_layers = P_layers_in[::-1]
    q_layers = q_layers_in[::-1]
    print "z_layers", z_layers

    l_satur  = [0,0]
    for pool in xrange(2):
        l_satur[pool]  = len(z_layers) - l_satur_in[pool]
    l_strato  = len(z_layers) - l_strato_in 


    for pool in xrange(2):
        l2_satur[pool]  = np.argmin( abs( z_layers[1::2] - z_layers[l_satur[pool]] ) )
    l2_strato = np.argmin( abs( z_layers[1::2] - z_layers[l_strato] ) )

    print ""
    print "Writing atmospheric files..."
    for pool in xrange(2) :
        write_irradiance( BASE_NAME[pool], sol * factor[pool] )
        write_szen( BASE_NAME[pool] )
        write_TPprof(  BASE_NAME[pool], P_layers[pool], T_layers[pool], sst[pool] )
        write_H2Oprof( BASE_NAME[pool], P_layers[pool], q_layers[pool] )
        write_gasprof( BASE_NAME[pool], P_layers[pool] )
        write_O3prof(  BASE_NAME[pool], P_layers[pool] )
        write_surf(  BASE_NAME[pool], albedo )

    print ""
    print "Calculating radiation fields in the warm pool..."
    rad_vis_w = calc_rad_sol( 0, [0, l2_satur[0], l2_strato]  )
    print "rad_vis", rad_vis_w
    rad_mir_w = calc_rad_mir( 0, [0, l2_satur[0], l2_strato] )
    print "rad_mir", rad_mir_w
    print "NET: ", rad_mir_w + rad_vis_w
    
    print ""
    print "Calculating radiation fields in the cold pool..."
    rad_vis_c = calc_rad_sol( 1, [0, l_satur[1], l_strato] )
    print "rad_vis", rad_vis_c
    rad_mir_c = calc_rad_mir( 1, [0, l_satur[1], l_strato] )
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


