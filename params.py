from scipy import constants
import numpy as np

FACTOR1 = (1.0/np.pi)
#FACTOR1 = (0.8/np.pi)
#FACTOR1 = (1.0/np.pi)


#---------------------------------
# Physical Constants
#---------------------------------
RH_surf = 0.75
RR      = constants.R    
Rs_dry  = 287.0          # J kg^-1 K^-1
LL      = 2.26476e6      # J kg^-1 
c_p     = 1.006e3        # J kg^-1 K^-1

#---------------------------------
# Planetary Parameters
#---------------------------------
P_surf   = 1.013e5      # Pa
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
#BETA     = 1.05          # index for mixing
#BETA     = 1.6          # index for mixing
zz_TI_tp = np.array([2000.0, 1500.0])       # m
delta_zz_TI = 200. # m
zz_TI_btm = zz_TI_tp - delta_zz_TI


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
FF_o_w    = -23.6       # W m^-2
FF_o_c    = -23.6       # W m^-2
FF_m_s_0  = -8.6   # W m^-2
LF_m_q_0  = -10.2  # W m^-2



#---------------------------------
# computational parameters (arbitrary)
#---------------------------------
PRECISION = 1e-8 # arbitrary
dT_dloop  = 1.0  # K
