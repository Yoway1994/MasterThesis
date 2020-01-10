import numpy as np
from numpy import pi, exp
import matplotlib.pyplot as plt

# Basic phyiscis parameters 

h = 6.6636* (10**-34) # Planck's constant (J*s)
c = 3* (10**8) # Speed of light (m/s)

# Cr4+:YAG crystal parameters
sigma_a = 22* (10**-23) # pump absorption cross-section at 1064nm (m^2) 
sigma_e = 6* (10**-23) # signal emission cross-section (m^2)
f_p = 0.17
f_L = 0.68 # the ratio between ESA and SE
sigma_esa = sigma_e* f_L # signal exicited state absorption
sigma_p_esa = sigma_a* f_p # pump ESA

QE = 0.094
tau_f = 4.5* (10**-6) # fluorescence lifetime(s)
tau_r = tau_f/ QE

FoM = 100 # crystal figure of merit
#alpha_p0 = 112 # small-signal differential pump absorption coefficient
alpha_p0 = 75
alpha_L0 = alpha_p0/ FoM # material losses at lambda_l
Nt = alpha_p0/ sigma_a # CDF doping concentration (m^-3)

# crystal fiber parameters
n_core = 1.82
n_inclad = 1.64

L = 0.06 # (m)
nz = 11 # number of fiber sections
z_span = [0, L] # span for ivp
z = np.linspace(0, L, nz) # series for ivp, matlab 0:dz:L
dz = L/ (nz -1)

lambda_p = 1064* (10**-9) # pumping wavelength (m)
lambda_l = 1443* (10**-9) # signal wavelength (m)

r_core = 10* (10**-6) # radius of core (m)
r_p = 10* (10**-6) # radius of pump mode field (m)
r_s = 10* (10**-6) # radius of signal mode field (m)
A_core = (r_core**2)* pi 
A_p = (r_p**2)* pi
A_s = (r_s**2)* pi

gammap = A_core/ A_p # pump confinement factor
gammas = A_core/ A_s # signal confinement factor

alpha_p = 2 # fiber attenuation constant for pump (m^-1)
alpha_s = 4.5 # fiber attenuation constant for signal (m^-1)

# multi-wavelength version settings
       
lambda_start = 1100* (10**-9) # ASE spectrum starting wavelength (m)
lambda_end = 1700* (10**-9) # ASE spectrum ending wavelength (m)
dlambda = 1* (10**-9) # wavelength interval (m)
lambda_array = np.linspace(lambda_start, lambda_end, 601) # wavelength array


# CDF Lineshape

A_f = 1.929* (10**-12)
A_g = 0.351* (10**18)
K = 0.127249
lambda0 = 1.4* (10**-6)

lineshape_f = A_f/ (lambda_array**2)* exp( -((lambda_array - lambda0)/ (K* lambda_array))**2)
lineshape_f = lineshape_f/ (np.trapz(lineshape_f, lambda_array))       # Normalization: integration equals to unity
lineshape_g = A_g* (lambda_array**3)* exp( -((lambda_array - lambda0)/ (K* lambda_array))**2)
sigma_s_e = sigma_e* lineshape_g # emission cross-section fx of wavelengths

# Laser cavity parameters

R1_p = 0.08 # pump reflectance at CDF input end 
R2_p = 0.08 # pump reflectance at CDF output end 
R1_s = 0.08 # reflectance at CDF input end 
R2_s = 0.08 # signal reflectance at CDF output end

cpe_s = 0.42 # couple efficiency for signal
cpe_p = 0.8 # couple efficiency for pump
cpe_pf = 1

pump_f = 1.4 # (W)
pump_b = 0 # (W)

pf = pump_f* (1 - R1_p)* exp( -gammap* sigma_a* Nt* z) # spatial-resolved forward pump
pb = pump_b* (1 - R2_p)* exp( -gammap* sigma_a* Nt* z[::-1]) # spatial-resolved backward pump
sf = np.zeros((len(z), len(lambda_array))) # (W)
sb = np.zeros((len(z), len(lambda_array))) # (W)

r_sp = (1 - (n_inclad/ n_core))/ 2 # Captured ratio of spontaneous emission

# experiment data
pump_exp = [0.025, 0.060, 0.120, 0.204, 0.305, 0.430, 0.560, 0.710, 0.865, 1.010, 1.160, 1.310]
ASE_exp = [21.6, 51.3, 97.7, 153.3, 221.4, 288, 358, 411, 458, 493, 518, 520]
RsP_exp = np.array([0.0182, 0.0578, 0.167, 0.495, 1.64, 5.1, 12.7, 27, 45.5, 67.6, 87, 105]) * 10**-3