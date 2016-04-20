### MKS unit ###
#import find_R
from scipy.optimize import fsolve
from scipy.optimize import minimize
import numpy as np
from params import *
#import func_to_minimize_ver4
import func_to_minimize as func_to_minimize

#from subroutines import *
#from params import *
#from radiation import *
#from cloud import *
#import mks as unit
#import sys

OUTPUTFILE = "survey_F0.85"
#OUTPUTFILE = "survey_test0"
#qq_TI_c = 4e-3
#aa_w    = 0.4

SST0 = np.array([ 308, 4.5]) 
#SST0 = np.array([ 293., 290.]) 

#SST0   = np.array([ 290.15, 289.]) 
SST_NUM = 4
SST_array = np.linspace(280, 310, SST_NUM)
SST_mesh_w, SST_mesh_c = np.meshgrid( SST_array, SST_array )
points = np.c_[ SST_mesh_w.flatten(), SST_mesh_c.flatten() ]

#=============================================================================
#def make_gridR( q_TI, a_w ):
#
#    values = np.zeros( [6, len(points)] )
#    #-----------------------------------------------------------
#    for ii in xrange( len(points) ):
#
#        sst_w = points[ii][0]
#        sst_c = points[ii][1]
#        
#        if ( sst_c <= sst_w ):
#            values.T[ii] = find_R.find_R( np.array([sst_w, sst_c]), q_TI, a_w, verbose=False )
#
#    #-----------------------------------------------------------
#    return points, values




#=============================================================================
# main
#=============================================================================
if __name__ == "__main__":

    
    aw_array = [0.33, 0.4, 0.5, 0.67]
    qq_array = [2e-3, 4e-3, 6e-3, 8e-3]

#    aw_array = [0.4]
#    qq_array = [4e-3]

    with open( OUTPUTFILE, 'w') as f:
        f.write( "# aa_w \t qq_T_c \t SST_w \t SST_c \t RR_w_0 \t RR_w_1 \t RR_w_2 \t RR_c_0 \t RR_c_1 \t RR_c_2 \t MM_w \t MM_c \t SS_w \t SS_c \t EE_w \t EE_c \t residuals \n" )

    for aa_w in aw_array :

        for qq_TI_c in qq_array :      

#            print ""
#            print "CREATING GRIDDATA OF RADIATION..."
#            print ""
#            points, values = make_gridR( qq_TI_c, aa_w )

            print ""
            print "SOLVING THE EQUATION for...( aa_w, qq_TI_c )=", aa_w, qq_TI_c
            print ""

            print ""
            print "#----------------------------------------"
            print " # for qq_TI =", qq_TI_c, ", aa_w =", aa_w
            print "#----------------------------------------"

            params = ( qq_TI_c, aa_w, False, OUTPUTFILE )
#            output = fsolve( func_to_minimize_ver2.myfunc, SST0, args=params )

#            bnds = ((273.15, 273.15), (273.15, 373.15))
#            cons = ({'type': 'ineq', 'fun': lambda x:  x[0] - x[1] })

            output = fsolve( func_to_minimize.myfunc, SST0, args=params )
            residual = func_to_minimize.myfunc( output, qq_TI_c, aa_w, True, OUTPUTFILE )
#            output = fsolve( func_to_minimize_ver2.myfunc, SST0, args=params )
#            residual = func_to_minimize_ver2.myfunc( output, qq_TI_c, aa_w, True, OUTPUTFILE )

            print output
            print "residual", residual

            with open( OUTPUTFILE, 'a') as f:
#                f.write( str(output['success']) + "\n" )
                f.write( "\n" )
#                f.write( str(aa_w) + "\t" )
#                f.write( str(qq_TI_c) + "\t" )
#                f.write( str(output['x'][0]) + "\t" )
#                f.write( str(output['x'][1]) + "\n" )




