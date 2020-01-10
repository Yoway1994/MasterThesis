# coding: utf-8

import os
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

def RK45(z, y0, n2, dfdx, d):
	x0 = z[iz]
	k1 = dfdx(z, y0, n2)
	x1 = x0 + d/2
	y1 = y0 + k1*d/2
	
	k2 = dfdx(x1, y1, n2)
	y2 = y0 + k2*d/2
	x2 = x0 + d/2
	
	k3 = dfdx(x2, y2, n2)
	y3 = y0 + k3* d
	x3 = x0 + d
    
	k4 = dfdx(x3, y3, n2)
    
	k = (k1 + 2*k2 + 2*k3 + k4)/6
	y = y0 + k*d
	x = x0 + d
	return y

N_Run = 2000
tol = 10**-2
Pump_f = np.linspace(0, pump_f, 11)[np.newaxis, :]
Pump_b = np.linspace(0, pump_b, 11)[np.newaxis, :]
z = np.linspace(0, L, nz)[:, np.newaxis]
Ssf = np.zeros(len(z))
Ssb = np.zeros(len(z))
SASEf = np.zeros(len(z))
SASEb = np.zeros(len(z))
# pump initialization
Ppf = Pump_f* (1 - R1_p)* Ti1_p* exp(-gammap* Nt* sigma_a* z)
Ppb = Pump_b* (1 - R2_p)* Ti2_p* exp(-gammap* Nt* sigma_a* (L - z))

for it in range(N_Run):
	Ppf[0, :] = Ppf[0]* (1-R1_p)* Ti1_p + Ppb[0, :]* R1_p* (Ti1_p**2)
	Ssf[0] = Ssb[0]* R1_s* Ti1_s**2 + SASEb[0]* R1_s* (Ti1_s**2)
	SASEf[0] = 0
	# Update N2
	N2 = calcN2((Ppf + Ppb), (Ssf + Ssb), (SASEf + SASEb))
	# Integration
	for iz in range(nz-1):
		Ppf[iz+1] = RK45(z, Ppf[iz], N2[iz], dPpdz, dz)
		Ssf[iz+1] = RK45(z, Ssf[iz], N2[iz], dSsdz, dz)
		SASEf[iz+1] = RK45(z, SASEf[iz], N2[iz], dSASEdz, dz)
		# backward integration
		# Boundary condition            
	Ppb[-1] = pump_b* (1-R2_p)* Ti2_p + Ppf[-1]* R2_p* (Ti2_p**2)
	Ssb[-1] = Ssf[-1]* R2_s* Ti2_s**2 + SASEf[-1]* R2_s* (Ti2_s**2)
	SASEb[-1] = 0
		# Update N2
	N2 = calcN2((Ppf + Ppb), (Ssf + Ssb), (SASEf + SASEb))
		# Integration
	for iz in range(nz):
		Ppb[iz-1] = RK45(z, Ppb[iz], N2[iz], dPpdz, dz)
		Ssb[iz-1] = RK45(z, Ssb[iz], N2[iz], dSsdz, dz)
		SASEb[iz-1] = RK45(z, SASEb[iz], N2[iz], dSASEdz, dz)

