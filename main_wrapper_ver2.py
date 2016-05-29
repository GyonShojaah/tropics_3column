### MKS unit ###
import find_R
from scipy.optimize import fsolve
from scipy.optimize import minimize
import numpy as np
from params import *
import func_to_minimize_ver2 as func_to_minimize

OUTPUTFILE = "survey_F1.0-0.9"
#OUTPUTFILE = "survey_test0"

#SST0 = np.array( [ 34. + 273.15, np.sqrt( 34. - 30. ) ] )
SST0 = np.array([ 300.0, 1.0])
#SST0 = np.array([ 293., 290.]) 

SST_array_w    = np.linspace(290.15, 320.15, 5)
SST_array_diff = np.linspace(0,3,5)
SST_mesh_w, SST_mesh_diff = np.meshgrid( SST_array_w, SST_array_diff )
points = np.c_[ SST_mesh_w.flatten(), SST_mesh_diff.flatten() ]


#=============================================================================
# main
#=============================================================================
def make_gridR( q_TI, a_w ):

    values = np.zeros( [6, len(points)] )
    #-----------------------------------------------------------
    for ii in xrange( len(points) ):

        sst_in_w    = points[ii][0]
        sst_in_diff = points[ii][1]
        
        values.T[ii] = find_R.find_R( np.array([sst_in_w, sst_in_diff]), q_TI, a_w, verbose=False )

    #-----------------------------------------------------------
    return points, values



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


            print ""
            print "#---------------------------------------------------------------"
            print "CREATING GRIDDATA OF RADIATION..."
            print "#---------------------------------------------------------------"
            print ""
            points, values = make_gridR( qq_TI_c, aa_w )

            print ""
            print "#---------------------------------------------------------------"
            print "# SOLVING THE EQUATION for...( aa_w, qq_TI_c )=", aa_w, qq_TI_c
            print "#---------------------------------------------------------------"
            print ""

            params = ( qq_TI_c, aa_w, points, values, False, OUTPUTFILE )

            output = fsolve( func_to_minimize.myfunc, SST0, args=params )
            print output

            residual = func_to_minimize.myfunc( output, qq_TI_c, aa_w, points, values, True, OUTPUTFILE )
            print "residual", residual

            with open( OUTPUTFILE, 'a') as f:
                f.write( "\n" )

        with open( OUTPUTFILE, 'a') as f:
            f.write( "\n\n" )


