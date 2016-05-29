
### MKS unit ###
from scipy.optimize import fsolve
from scipy.optimize import minimize
import numpy as np
from params import *
import func_to_minimize_Ronly as func_to_minimize

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


OUTPUTFILE = "check_test"
#OUTPUTFILE = "outputfile"

SST = np.array( [ 39. + 273.15, np.sqrt( 39. - 37. ) ] )
#SST = np.array( [ 34. + 273.15, np.sqrt( 34. - 30. ) ] )
#SST = np.array( [ 302., np.sqrt( 1.8 ) ] )


#=============================================================================
# main
#=============================================================================
if __name__ == "__main__":

    
#    aw_array = [0.33, 0.4, 0.5, 0.67]
#    qq_array = [2e-3, 4e-3, 6e-3, 8e-3]

    aa_w = 0.5
    qq_TI_c = 6e-3
    
    sst_w   = np.linspace( 20., 50., 20 ) + 273.15
    sst_dif = np.logspace( np.log10(0.1), np.log10(15), 20 )

    for ii in xrange(len(sst_w)):
        for jj in xrange(len(sst_dif)):

            SST_in = np.array( [ sst_w[ii], sst_dif[jj] ] )
            func_to_minimize.myfunc( SST_in, qq_TI_c, aa_w, True, OUTPUTFILE )


    


