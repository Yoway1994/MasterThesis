# modules import
from numpy import pi, exp
import numpy as np

# Basic phyiscis parameters 
h = 6.6636* (10**-34) # Planck's constant (J*s)
c = 3* (10**8) # Speed of light (m/s)

# crystal fiber parameters #
# length
L = 0.2 # (m)
nz = 1001 # number of fiber sections
z_span = [0, L] # span for ivp
z = np.linspace(0, L, nz) # series for ivp, matlab 0:dz:L
dz = L/ (nz - 1)

# core
n_core = 1.81
n_inclad = 1.64

r_core = 8.5* (10**-6) # radius of core (m)
r_p = 8.5* (10**-6) # radius of pump mode field (m)
r_s = 8.5* (10**-6) # radius of signal mode field (m)
A_core = (r_core**2)* pi # area of core
A_p = (r_p**2)* pi # area of pump mode field
A_s = (r_s**2)* pi # area of signal mode field 
gammap = A_core/ A_p # pump confinement factor
gammas = A_core/ A_s # signal confinement factor

# propogaion
alpha_p = 1 # fiber attenuation constant for pump (m^-1)
alpha_s = 1 # fiber attenuation constant for signal (m^-1)

# Cr4+:YAG crystal parameters
# wavelength
lambda_p = 1064* (10**-9) # pumping wavelength (m)
lambda_l = 1447* (10**-9) # signal lasing wavelength (m)
lambda_c = 1420* (10**-9) # center wavelength of emission cross-section (m)
sigma_a = 22* (10**-23) # pump absorption cross-section at 1064nm (m^2) 
sigma_e = 5.934* (10**-23) # signal emission cross-section (m^2)
f_p = 0.04 # the ESA ratio of pump wavelength
f_L = 0.4 # the ratio between ESA and SE
sigma_p_esa = sigma_a* f_p
sigma_esa = sigma_e* f_L # signal exicited state absorption

# lifetime
tau_f = 4.5* (10**-6) # fluorescence lifetime(s)
# tau_r = (8* pi* (n_core**2)* c* np.trapz(lambda2,sigma_s_e2./(lambd_l**4)))**(-1)
QE = 0.13
tau_r = tau_f/ QE

# concentration
#alpha_p0 = 150 # small-signal differential pump absorption coefficient
#Nt = alpha_p0/ sigma_a # CDF doping concentration (m^-3)
Nt = 10* 10**23
#FoM = 100 # crystal figure of merit
#alpha_L0 = alpha_p0/ FoM # material losses at lambda_l

# lineshape function
K_e = 0.134
lineshape_e = exp(-((lambda_l - lambda_c)/ (K_e* lambda_l))**2) # normalized lineshape
lineshape_esa = lineshape_e #
sigma_s_e = sigma_e* lineshape_e # emission cross-section (m^2)
sigma_s_esa = sigma_esa* lineshape_esa #
dlambda = 0.1* (10**-9) # the slot width at center wavelength

r_sp = (1 - (n_inclad/ n_core))/2 # 
Ssp = (4* pi* n_core*(n_core - n_inclad)* h* (c**2)* sigma_s_e/ (lambda_l**5))* r_sp
# Laser cavity parameters #
# input power
pump_f = 20 # maxium input forward pump power(W)
pump_b = 0 # (W)

# cavity
R1_p = 0.085 # pump reflectance at CDF input end 
R2_p = 0.085 # pump reflectance at CDF output end 
R1_s = 0.085 # signal reflectance at CDF input end 
R2_s = 0.085 # signal reflectance at CDF output end
Ti1_p = 1 # single-pass transmission between IC and input fiber end (for pump)
Ti2_p = 1 # single-pass transmission between OC and output fiber end (for pump)
Ti1_s = 1 # single-pass transmission between IC and input fiber end (for signal)
Ti2_s = 1 # single-pass transmission between OC and output fiber end (for signal)

T_filter = 1 # transmission of filter

# couple efficiency
cpe_s = 1 # couple efficiency for signal
cpe_p = 1 # couple efficiency for pump

lambda_start = 1100* (10**-9) # ASE spectrum starting wavelength (m)
lambda_end = 1700* (10**-9) # ASE spectrum ending wavelength (m)
dlambda = 1* (10**-9) # wavelength interval (m)
lambda_array = np.linspace(lambda_start, lambda_end, 601) # wavelength array
