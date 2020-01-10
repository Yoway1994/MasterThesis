# modules import
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, exp, log10
from laser_p import *
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

# laser rate equation

# pumping rate
def R03(Ppf, Ppb):
	r03 = (sigma_a* lambda_p* (Ppf + Ppb))/ (A_p* h* c)
	return r03
	
# total stimulation emission rates	
def W21(Ssf, Ssb, SASEf, SASEb):
	w21 = (dlambda* np.dot((Ssf + Ssb + SASEf + SASEb), (sigma_s_e* lambda_l)))/ (A_s* h* c)
	return w21

# N2 density	
def N2_func(Ppf, Ppb, Ssf, Ssb, SASEf, SASEb):
	n2 = Nt* (R03(Ppf, Ppb)/ (W21(Ssf, Ssb, SASEf, SASEb) + (1/tau_f) + R03(Ppf, Ppb)))
	return n2
	
def f(z0):
	k = interp1d(z, N2_func(Ppf_, Ppb_, Ssf_, Ssb_, SASEf_, SASEb_), fill_value="extrapolate")
	return k(z0)

# power evalation

def dPpfdz(z, Ppf):
	dPpfdz = (-Ppf)* (gammap* (sigma_a* (Nt - f(z)) + (sigma_p_esa* f(z))) + alpha_p)
	return dPpfdz
	
def dPpbdz(z, Ppb): 
	dPpbdz = Ppb* (gammap* (sigma_a* (Nt - f(z)) + (sigma_p_esa* f(z))) + alpha_p)
	return dPpbdz

def dSsfdz(z, Ssf):
	dSsfdz = Ssf* (gammas* (sigma_s_e* (1 - f_L)* f(z)) - alpha_s)
	return dSsfdz
	
def dSsbdz(z, Ssb):
	dSsbdz = (-Ssb)* (gammas* (sigma_s_e* (1 - f_L)* f(z)) - alpha_s)
	return dSsbdz

def dSASEfdz(z, SASEf):
	dSsfdz = SASEf* (gammas* (sigma_s_e* (1 - f_L)* f(z)) - alpha_s) + (r_sp* h* c/ lambda_array* A_core* f(z)/ tau_r* lineshape_f)
	return dSsfdz
	return dSsfdz
	
def dSASEbdz(z, SASEb):
	dSsbdz = (-SASEb)* (gammas* (sigma_s_e* (1 - f_L)* f(z)) - alpha_s) - (r_sp* h* c/ lambda_array* A_core* f(z)/ tau_r* lineshape_f)
	return dSsfdz
	return dSsbdz

def power_evalation(Ppf, Ppb, Ssf, Ssb, SASEf, SASEb, input_p):	
	# forward 
	sol_dPpfdz = solve_ivp(dPpfdz, z_span, input_p* (1 - R1_p) + R1_p* Ppb[0], t_eval = z)	
	sol_dSsfdz = solve_ivp(dSsfdz, z_span, [(Ssb[-1] + SASEb[-1])], t_eval = z)
	sol_dSASEfdz = solve_ivp(dSASEfdz, z_span, [SASEb[-1]], t_eval = z)
	Ppf = sol_dPpfdz.y[0]
	Ssf = sol_dSsfdz.y[0]
	SASEf = sol_dSASEfdz.y[0]
	# backward
	sol_dPpbdz = solve_ivp(dPpbdz, z_span[::-1], pb* (1 - R2_p) + R2_p* Ppf[-1], t_eval = z[::-1])
	sol_dSsbdz = solve_ivp(dSsbdz, z_span[::-1], [R2_s*(Ssf[-1] + SASEf[-1])], t_eval = z[::-1])
	sol_dSASEbdz = solve_ivp(dSsbdz, z_span[::-1], [SASEf[-1]], t_eval = z[::-1])
	Ppb = sol_dPpbdz.y[0]
	Ssb = sol_dSsbdz.y[0]
	SASEb = sol_dSASEbdz.y[0]
	Ppb = Ppb[::-1]
	Ssb = Ssb[::-1]
	SASEb = SASEb[::-1]
	
	return Ppf, Ppb, Ssf, Ssb, SASEf, SASEb
#
# power calculating function
# output ASE power
def Psf_out(Ssf):
	out = (1 - R1_s)* Ssf[-1]
	return out

	
	
	
# iteration

power = [] # list for saving result
tol = 1* (10**-10) # tolerance
N_Run = 20 # maximum iteration times

sf = np.zeros(len(z))
sb = np.zeros(len(z))
sasef = np.zeros(len(z))
saseb = np.zeros(len(z))


for p_f in np.linspace(0, pump_f, 21):
	global Ppf_, Ppb_, Ssf_, Ssb_, SASEf_, SASEb_
	pf = p_f* (1 - R1_p)* exp( -gammap* sigma_a* Nt* z)
	pb = pump_b* (1 - R2_p)* exp( -gammap* sigma_a* Nt* z) # pump_b fixed
	Ppf_ = np.copy(pf)
	Ppb_ = np.copy(pb)
	Ssf_ = np.copy(sf)
	Ssb_ = np.copy(sb)
	SASEf_ = np.copy(sasef)
	SASEb_ = np.copy(saseb)
	
	# tolerance control testing
	for i in range(N_Run):
		Ppf_, Ppb_, Ssf_, Ssb_, SASEf_, SASEb_ = power_evalation(Ppf_, Ppb_, Ssf_, Ssb_, SASEf_, SASEb_, pf)
		print(i)
		
	# store data into power list
	dict = {}
	dict['pump'] = p_f
	dict['Ppf'] = Ppf_
	dict['Ppb'] = Ppb_
	dict['Ssf'] = Ssf_
	dict['Ssb'] = Ssb_
	dict['SASEf'] = SASEf_
	dict['SASEb'] = SASEb_
	dict['N2'] = N2_func(Ppf_, Ppb_, Ssf_, Ssb_, SASEf_, SASEb_)
	power.append(dict)
	
pump_x = []
Ssf_y = []
Gain_y = []	
	
for i in power:
	pump_x.append(i.get('pump'))
	Ssf_y.append(i.get('Ssf')[-1])
	
	gain = i.get('N2')* sigma_e* (1 - f_L) # differential gain, g(lambda,z)
	Gain_y.append(10*log10(exp(np.trapz(gammas*gain, z))))

plt.figure('Ssf')
plt.plot(pump_x, np.array(Ssf_y)*(1 - R2_s))
plt.plot()
plt.figure('Gain')
plt.plot(pump_x, Gain_y)
plt.plot()
plt.show()