### MKS unit ###
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.optimize import root
import numpy as np
from params import *
import func_to_minimize as func_to_minimize


#=============================================================================
# main
#=============================================================================
if __name__ == "__main__":
    
    aw_array = [0.33, 0.4, 0.5, 0.67]
#    qq_array = [2e-3, 4e-3, 6e-3, 8e-3]
    qq_array = [1e-3, 4e-3, 6e-3, 9e-3]

#    aw_array = [0.4]
#    qq_array = [4e-3]


    with open( OUTPUTFILE, 'w') as f:
        f.write( "# aa_w \t qq_T_c \t SST_w \t SST_c \t z_tropopause \t RR_w_0 \t RR_w_1 \t RR_w_2 \t RR_c_0 \t RR_c_1 \t RR_c_2 \t MM_w \t MM_c \t SS_w \t SS_c \t EE_w \t EE_c \t residuals \n" )

    for aa_w in aw_array :

        for qq_TI_c in qq_array :      

            SST_ini = SST0

            print ""
            print "#---------------------------------------------------------------"
            print "# SOLVING THE EQUATION for...( aa_w, qq_TI_c )=", aa_w, qq_TI_c
            print "#---------------------------------------------------------------"
            print ""

            params = ( qq_TI_c, aa_w, False, OUTPUTFILE )

#            output = fsolve( func_to_minimize.myfunc, SST0, args=params )
            output = root( func_to_minimize.myfunc, SST_ini, args=params )['x']
            print output

            SST_ini = output

            residual = func_to_minimize.myfunc( output, qq_TI_c, aa_w, True, OUTPUTFILE )
            print "residual", residual

            with open( OUTPUTFILE, 'a') as f:
                f.write( "\n" )

        with open( OUTPUTFILE, 'a') as f:
            f.write( "\n\n" )


