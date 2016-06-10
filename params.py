from scipy import constants
import numpy as np

OUTPUTFILE = "surveyR_F0.95-0.92"
FACTOR = np.array( [ (0.95/np.pi), (0.92/np.pi) ] )

SST0 = np.array( [ 45. + 273.15, np.sqrt( 1. ) ] )
SOL      = 1360.6718     # W/m^2
ALBEDO   = 0.05


#---------------------------------
# Physical Constants
#---------------------------------
RH_surf = 0.75
RR      = constants.R    
Rs_dry  = 287.0          # J kg^-1 K^-1
#LL      = 2.26476e6      # J kg^-1 
LL      = 2.50e6      # J kg^-1 
c_p     = 1.006e3        # J kg^-1 K^-1

#---------------------------------
# Planetary Parameters
#---------------------------------
P_surf   = 0.98e5      # Pa
gg       = constants.g  # 
MU_H2O   = 18.0e-3      # kg
MU_atm   = 28.8e-3      # kg
T_strato = 195          # K ( Miller's paper )
GAMMA    = 1.4
DELTA    = ( GAMMA / ( GAMMA-1. ) ) # P/P0 = (T/T0)**DELTA for dry adiabat


#---------------------------------
# parameters in hands
#---------------------------------
dT_surf  = 1.           # K ( Miller's paper )
ALPHA    = 1.02         
BETA     = 1.2          # index for mixing
zz_TI_top = np.array([2000.0, 1500.0])       # m
delta_zz_TI = 200. # m
zz_TI_btm = zz_TI_top - delta_zz_TI


#---------------------------------
# Radiative forcing
#---------------------------------
dR_all   = -5.
dR_c     = -15.
dR_w_TOA = -7.
dR_w_satur = -27.

#---------------------------------
# Divergence
#---------------------------------
FF_o_w    = -23.6  # W m^-2
FF_o_c    = -23.6  # W m^-2
FF_m_s_0  = -8.6   # W m^-2
LF_m_q_0  = -10.2  # W m^-2


#---------------------------------
# computational parameters (arbitrary)
#---------------------------------
PRECISION = 1e-8 # arbitrary
dT_dloop  = 1.0  # K
