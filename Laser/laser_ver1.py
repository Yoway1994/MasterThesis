# inital power setting 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from numpy import pi, exp, log10
from laser_p import *
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
###################################################################################################################
#laser cavity function
###################################################################################################################
# totalstimulation emission rates (s^-1)
def W21(Ss):
	w21 = np.dot(Ss, (sigma_s_e* lambda_l))/ (A_s* h* c)
	# Ss* (gammas* sigma_s_e* lambda_l)/ (A_s* h* c)
	return w21

# pumping rate (s^-1)
def R03(Pp):
    r03 = Pp* (gammap* sigma_a* lambda_p)/ (A_p* h* c)
    return r03

# N2 distribution
def calcN2(Pp, Ss, SASE):
    cn2 = Nt* (R03(Pp)/ (W21(Ss + SASE) + (1/ tau_f) + R03(Pp))) 
    return cn2
	
def f(z0):
	k = interp1d(z, calcN2((Ppf+Ppb), (Ssf+Ssb), (SASEf+SASEf)), fill_value="extrapolate")
	return k(z0)
	
###################################################################################################################	
# pump absorption coefficient (m^-1)
def a_p(z):
    ap = gammap* (sigma_a* (Nt - f(z)) + sigma_p_esa* f(z))
    return ap

# pump evolution
def dPpdz(z, Pp):
    dPpdz = -Pp* (a_p(z) + alpha_p) 
    return dPpdz
	
def dPpdz_back(z, Pp):
    dPpdz = Pp* (a_p(z) + alpha_p) 
    return dPpdz	
###################################################################################################################
# gain coefficient (m^-1)
def g(z):
    g = gammas* f(z)* (sigma_s_e - sigma_s_esa)
    return g 

# signal evolution
def dSsdz(z, Ss):
	dSsdz = Ss* (g(z) - alpha_s)
	return dSsdz
	
def dSsdz_back(z, Ss):
    dSsdz = Ss* (-g(z) + alpha_s)
    return dSsdz
###################################################################################################################
# SE Power per mode per polarization per direction per unit length (SMF or FMF)
def Psp(z):
    Psp = f(z)* h* (c**2)* sigma_s_e/ (lambda_l**3)* dlambda
    return Psp

# ASE evolution
def dSASEdz(z, SASE):
	dSASEdz = SASE* (g(z) - alpha_s)+ 2* Psp(z)
	print(SASE[-1])
	# dSASEdz = SASE* (g(z) - alpha_s)
	return dSASEdz

def dSASEdz_back(z, SASE):
	dSASEdz = -SASE*(g(z) - alpha_s)- 2* Psp(z)
	# dSASEdz = SASE*(-g(z) + alpha_s)
	return dSASEdz
###################################################################################################################
def Ssf_out(Ssf):
	s_out = Ssf[-1]* (1-R2_s)* Ti2_s
	return s_out 

########################################################################################

Nth = alpha_s/gammas/(sigma_s_e-sigma_s_esa) # N2 threshold for g>alpha_s

# initial value

Ssf_i = np.zeros(len(z)) # Forward signal
Ssb_i = np.zeros(len(z)) # Backward signal
SASEf_i = np.zeros(len(z)) # Forward ASE (W)
SASEb_i = np.zeros(len(z)) # Backward ASE (W)
#Ppf_i = np.zeros(len(z))
#Ppb_i = np.zeros(len(z))

p_interval = 21 # number of input pump points
Pump_f = np.linspace(0, pump_f, p_interval) # the series of iuput pump power from 0 ~ 2W
Pump_b = np.linspace(0, pump_b, p_interval) # pump power series
power = []

########################################################################################
for ip in range(p_interval):
	# the initial pump power spatial distribution along crystal fiber 
	Ppf_i = Pump_f[ip]* (1 - R1_p)* Ti1_p* exp(-gammap* Nt* sigma_a* z) 
	Ppb_i = Pump_b[ip]* (1 - R2_p)* Ti2_p* exp(-gammap* Nt* sigma_a* (L - z))
	global Ssf, Ssb, SASEf, SASEb, Ppf, Ppb
	Ssf = np.copy(Ssf_i)
	Ssb = np.copy(Ssb_i)
	SASEf = np.copy(SASEf_i)
	SASEb = np.copy(SASEb_i)
	Ppf = np.copy(Ppf_i)
	Ppb = np.copy(Ppb_i)
	# propagation inside cavity
	for i in range(100):
		print(i)
		# forward Runge-Kutta
		sol_dPpdz = solve_ivp(dPpdz, z_span, Ppf_i + Ppb[0]* R1_p* (Ti1_p**2), t_eval = z)
		sol_dSsdz = solve_ivp(dSsdz, z_span, [(Ssb[0] + SASEb[0])* R1_s* (Ti1_s**2)], t_eval = z)
		sol_dSASEdz = solve_ivp(dSASEdz, z_span, [SASEb[0]* R1_s* (Ti1_s**2)], t_eval = z)
		Ppf = sol_dPpdz.y[0]
		Ssf = sol_dSsdz.y[0]
		SASEf = sol_dSASEdz.y[0]

		# backward Runge-Kutta
		sol_dPpdz_back = solve_ivp(dPpdz, z_span[::-1], Ppb_i + Ppf[-1]* R2_p* (Ti2_p**2), t_eval = z[::-1])
		sol_dSsdz_back = solve_ivp(dSsdz, z_span[::-1], [(Ssf[-1] + SASEf[-1])* R2_s* (Ti2_s**2)], t_eval = z[::-1])
		sol_dSASEdz_back = solve_ivp(dSASEdz_back, z_span[::-1], [SASEf[-1]* R1_s* (Ti1_s**2)], t_eval = z[::-1])
		Ppb = sol_dPpdz_back.y[0]
		Ssb = sol_dSsdz_back.y[0]
		SASEb = sol_dSASEdz_back.y[0]
		#
		
	# data storge
	dict = {}
	dict['pump'] = Pump_f[ip]
	dict['Ppf'] = Ppf
	dict['Ppb'] = Ppb
	dict['Ssf'] = Ssf
	dict['Ssb'] = Ssb
	dict['SASEf'] = SASEf
	dict['SASEb'] = SASEb
	dict['N2'] = calcN2((Ppf+Ppb), (Ssf+Ssb), (SASEf+SASEb))
	# plt.plot(z, calcN2((Ppf+Ppb), (Ssf+Ssb), (SASEf+SASEb)))
	power.append(dict)

#######################################################################################
Gain = []
sf = []
N2_list = []

for i in range(len(power)):
	#
	gain = power[i].get('N2')* sigma_e* (1 - f_L) # differential gain, g(lambda,z)
	Gain.append(10*log10(exp(np.trapz(gammas*gain, z))))
	#
	sf.append(power[i].get('Ssf')[-1])
	#
	N2_list.append(power[i].get('N2'))
	
# Roundtrip gain and loss
Loss_rt = -10*log10(R1_s* (Ti1_s**2)* R2_s* (Ti2_s**2)* exp(-alpha_s*2*L)) # round-trip loss, dB	
Loss_rt = np.ones(len(Pump_f))
plt.figure('Gain & Loss')
plt.plot(Pump_f, Gain, Pump_f, Loss_rt)
plt.xlabel('Input pump power (mW)')
plt.ylabel('Gain & round-trip Loss(dB)')
# 
plt.figure('')
plt.plot(Pump_f, sf)
#
plt.figure('N2')
for i in N2_list:
	plt.plot(z, i)
#######################################################################################
plt.show()
	