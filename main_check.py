### MKS unit ###
from scipy.optimize import fsolve
from scipy.optimize import minimize
import numpy as np
from params import *
import func_to_minimize

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


OUTPUTFILE = "check_F1.0_0.67-5e-3"
#OUTPUTFILE = "outputfile"

#SST = np.array( [ 33.5 + 273.15, np.sqrt( 33.5 - 28. ) ] )
#SST = np.array( [ 34. + 273.15, np.sqrt( 34. - 30. ) ] )
SST = np.array( [ 302., np.sqrt( 1.8 ) ] )


#=============================================================================
# main
#=============================================================================
if __name__ == "__main__":

    
#    aw_array = [0.33, 0.4, 0.5, 0.67]
#    qq_array = [2e-3, 4e-3, 6e-3, 8e-3]

    aa_w = 0.67
    qq_TI_c = 5e-3
    
    with open( OUTPUTFILE, 'w') as f:
        f.write( "# aa_w \t qq_T_c \t SST_w \t SST_c \t RR_w_0 \t RR_w_1 \t RR_w_2 \t RR_c_0 \t RR_c_1 \t RR_c_2 \t MM_w \t MM_c \t SS_w \t SS_c \t EE_w \t EE_c \t z_tropopause \t residuals \n" )

    func_to_minimize.myfunc( SST, qq_TI_c, aa_w, True, OUTPUTFILE )

    with open( OUTPUTFILE, 'a') as f:
        f.write( "\n" )

    


