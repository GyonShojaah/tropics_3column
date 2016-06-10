
### MKS unit ###
from scipy.optimize import fsolve
from scipy.optimize import minimize
import numpy as np
from params import *
import func_to_minimize as func_to_minimize
import sys

#from subroutines import *
#from params import *
#from radiation import *
#from cloud import *
#import mks as unit
#import sys

#qq_TI_c = 4e-3
#aa_w    = 0.4

#SST0 = np.array([ 302.15, 298.15]) 
#SST0 = np.array([ 302.15, 2.0]) 
#SST0 = np.array([ 293., 290.]) 


OUTPUTFILE = "checkgrid_aw0.33_qTIc2"
#OUTPUTFILE = "outputfile"

SST = np.array( [ 39. + 273.15, np.sqrt( 39. - 37. ) ] )
#SST = np.array( [ 34. + 273.15, np.sqrt( 34. - 30. ) ] )
#SST = np.array( [ 302., np.sqrt( 1.8 ) ] )


#=============================================================================
# main
#=============================================================================
if __name__ == "__main__":

    
#    aw_array = [0.33, 0.4, 0.5, 0.67]
#    qq_array = [2e-3, 4e-3, 6e-3, 8e-3

    aa_w = 0.33
    qq_TI_c = 2e-3
    
    sst_w   = np.linspace( 25., 45., 20 ) + 273.15
    sst_dif = np.linspace( np.sqrt(0.1), np.sqrt(10.), 20 )

#    sst_w = np.array( [ 30. + 273.15 ] )
#    sst_dif = np.array( [ 1. ] )
#    print sst_dif
#    sys.exit()

    for ii in xrange(len(sst_w)):
        for jj in xrange(len(sst_dif)):

            SST_in = np.array( [ sst_w[ii], sst_dif[jj] ] )

            SST = np.array([ SST_in[0], SST_in[0]-SST_in[1]**2. ])
#            print "SST", SST

            params = ( qq_TI_c, aa_w, True, OUTPUTFILE )
            residual = func_to_minimize.myfunc( SST_in, *params )

            print ""

            with open( OUTPUTFILE, 'a') as f:
                f.write( "\n" )


        print ""

        with open( OUTPUTFILE, 'a') as f:
            f.write( "\n" )

    


