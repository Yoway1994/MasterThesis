#!/usr/bin/env python
# coding: utf-8

# In[11]:


import numpy as np
import matplotlib.pyplot as plt
hemi = np.array([0.1, 4.3, 11, 18, 23.5, 31.5, 39, 46, 53, 59, 64])
coli = np.array([4.4, 9.5, 14.8, 19.3, 24.5, 29.1, 33.6, 37.9])
hp = [20, 44, 79, 112, 145, 175, 196]
tl = [39, 35, 31.2, 25.7, 21, 15.4, 9.9, 5.4, 2.1]

pump_hemi = np.linspace(50, 550, 11)
pump_coli = np.linspace(150, 500, 8)
pump_hp = [435, 800, 1330, 1830, 2300, 2660, 2900]
pump_tl = [1261, 1142, 1010, 874, 739, 610, 481, 370, 274]

Pth_hp = 189
Pth_hemi = 51
Pth_coli = 101
Pth_tl = 210

x0 = np.linspace(0, Pth_tl, 10)
y0 = np.zeros(len(x0))

x = np.linspace(Pth_tl, 1300, 100)
y = 0.0381*(x-Pth_tl)

plt.figure()
ax = plt.gca()
sp = ['right', 'left', 'top', 'bottom']
for i in sp:
    ax.spines[i].set_color('black')
    ax.spines[i].set_linewidth(2)
ax.tick_params(axis = "y", direction = "in")
ax.tick_params(axis = "x", direction = "in")
plt.plot(pump_tl, tl, 'o', color = 'black', marker = 's')
plt.plot(x0, y0, x, y, color = 'blue', linestyle = '--')
plt.legend(loc = 'best', labels = ['experiment','simulation'])
plt.xlabel('Pump power (mW)', fontsize = 14)
plt.ylabel('Output power (mW)', fontsize = 14)
plt.show()

