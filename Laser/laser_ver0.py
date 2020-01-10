# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from numpy import pi, exp
from laser_p import *

def Psp(N2):
	psp = N2* h* (c**2)* sigma_s_e* dlambda/ (lambda_l**3)
	return psp
	
def R03(Pp):
	r03 = Pp* gammap* sigma_a* lambda_p/(A_p* h* c)
	return r03
	
def W21(Ss):
	w21 = Ss* gammas* sigma_s_e* lambda_l/(A_s* h* c)
	return w21

def g(N2):
	gg = gammas* N2* (sigma_s_e - sigma_s_esa)
	return gg

def a_p(N2):
	ap = gammap* (sigma_a* (Nt- N2) + sigma_p_esa* N2)
	return ap
	
def calcN2(Pp, Ss, SASE):
	calcN2 = Nt* R03(Pp)/(W21(Ss + SASE) + 1/tau_f + R03(Pp))
	return calcN2
	
def dPpdz(z0, Pp, N2):
	dppdz = -Pp*(a_p(N2) + alpha_p)
	return dppdz
	
def dSsdz(z0, Ss, N2):
	dssdz = Ss*(g(N2) - alpha_s)
	return dssdz
	
def dSASEdz(z0, SASE, N2):
	dsasedz = SASE*(g(N2) - alpha_s) + 2*Psp(N2)
	return dsasedz

def Ssf_out(ssf):
	out = (1 - R1_s)* ssf[-1]
	return out

N_Run = 20000
tol = 10**-5
Pump_f = np.linspace(0, pump_f, 11)

sf = []

for i, i_p in enumerate(Pump_f):
	Ssf = np.zeros(len(z))
	Ssb = np.zeros(len(z))
	SASEf = np.zeros(len(z))
	SASEb = np.zeros(len(z))
	# pump initialization
	Ppf = i_p* (1 - R1_p)* Ti1_p* exp(-gammap* Nt* sigma_a* z)
	Ppb = pump_b* (1 - R2_p)* Ti2_p* exp(-gammap* Nt* sigma_a* (L - z))
    # iteration
	for it in range(N_Run):
        # forward integration
        # Boundary condition
		Ppf[0] = i_p* (1-R1_p)* Ti1_p + Ppb[0]* R1_p* (Ti1_p**2)
		Ssf[0] = Ssb[0]* R1_s* Ti1_s**2 + SASEb[0]* R1_s* (Ti1_s**2)
		SASEf[0] = 0
		# Update N2
		N2 = calcN2((Ppf + Ppb), (Ssf + Ssb), (SASEf + SASEb))
		# Integration
		for iz in range(nz-1):
			Ppf[iz+1] = Ppf[iz] + dz* dPpdz(z, Ppf[iz], N2[iz])
			Ssf[iz+1] = Ssf[iz] + dz* dSsdz(z, Ssf[iz], N2[iz])
			SASEf[iz+1] = SASEf[iz] + dz* dSASEdz(z, SASEf[iz], N2[iz])
		# backward integration
        # Boundary condition            
		Ppb[-1] = pump_b* (1-R2_p)* Ti2_p + Ppf[-1]* R2_p* (Ti2_p**2)
		Ssb[-1] = Ssf[-1]* R2_s* Ti2_s**2 + SASEf[-1]* R2_s* (Ti2_s**2)
		SASEb[-1] = 0
		# Update N2
		N2 = calcN2((Ppf + Ppb), (Ssf + Ssb), (SASEf + SASEb))
		# Integration
		for iz in range(nz):
			Ppb[iz-1] = Ppb[iz] + dz* dPpdz(z, Ppb[iz], N2[iz])
			Ssb[iz-1] = Ssb[iz] + dz* dSsdz(z, Ssb[iz], N2[iz])
			SASEb[iz-1] = SASEb[iz] + dz* dSASEdz(z, SASEb[iz], N2[iz])
		# breaking condition
		print(it)
	
	sf.append(Ssf[-1])


plt.plot(Pump_f, sf)
plt.show()