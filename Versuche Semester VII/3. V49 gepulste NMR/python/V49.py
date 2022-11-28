import numpy as np
import matplotlib.pyplot as plt
import uncertainties as unc
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
import scipy.constants as const
import sympy
import os
from tabulate import tabulate
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)


if os.path.exists("build") == False:
    os.mkdir("build")

if os.path.exists("build/plots") == False:
    os.mkdir("build/plots")



########MESSWERTE#######

#Justage

f = 21.73738 #Hz
phase = 1040

#shim?

deg_180 = 5.28 #mmicro seconds
deg_90 = 0.5 * deg_180

Temp_1 = 20.8
Temp_2 = 20.8

#T1-MEssung
t_1= np.array([1, 1.5, 2, 3, 5, 7.5, 10, 15, 20, 27.5, 35, 45, 60, 75, 150, 300, 500, 750, 1000, 1400, 1800, 2300, 3000, 4500, 6000, 7500, 9999]) #in ms
U_11 = np.array([3.075, 3.05, 3.05, 3.05, 3.1, 3.1, 3.075, 3.075, 3.05, 3.04, 3, 2.975, 2.925, 2.875, 2.72, 2.4, 2, 1.44]) * 0.5 # in V
U_12 = np.array([0.6, 0.375, 0.11, -0.26, -0.57, -1.085, -1.275, -1.425, -1.55])
U_1 = np.concatenate((U_11, U_12), axis=0)


#T2-Messung
t_2, U_2 = np.genfromtxt("python/data/scope_1.csv",skip_header= 2, skip_footer= 2220, unpack = True)

#Diffusionsmessung

t_diff = np.array([0.1, 0.5, 1, 2.5, 5, 7.5, 10, 12.5, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]) #in ms
U_diff = np.array([1.5, 1.487, 1.487, 1.468, 1.375, 1.15, 0.837, 0.475, 0.212, 0.160, 0.110, 0.073, 0.047, 0.041, 0.033, 0.029, 0.02, 0.027]) # in V

tau5_diff, Real_diff, Imag_diff = np.genfromtxt("python/data/scope_4.csv",skip_header= 5, unpack = True)


#####FUNKTIONEN#######

def rel_abw(theo,a):
    c = (theo - a)/theo
    print(f"Relative Abweichung in Prozent: {noms(c) * 100 :.4f} \pm {stds(c) * 100 :.5f}\n")





#####RECHNUNGEN#######

###### T1-Auswertung #####
print("\n\t##### T1- Auswertung ######")
def fit_T1(t, a, b, c):
    return a*(1-2*np.exp(-t/b))+c

params1, cov1 = curve_fit(fit_T1, t_1, U_1, p0 = (-1.6, 2839, -0.638 ))
cov1 = np.sqrt(np.diag(abs(cov1)))

uparams1 = unp.uarray(params1, cov1)

print(f"U_0 = {noms(uparams1[0]):.4f} \pm {stds(uparams1[0]):.4f} Volt\nT1 = {noms(uparams1[1]):.4f} \pm {stds(uparams1[1]):.4f} ms\nU_1 = {noms(uparams1[2]):.4f} \pm {stds(uparams1[2]):.4f} Volt\n")


x1 = np.linspace(t_1[0], t_1[-1], 1000)

plt.figure()
plt.plot(t_1, U_1,"x",label="Messwerte")
plt.plot(x1, fit_T1(x1, *params1))
#plt.xscale('log')
plt.xlabel(r"$\tau$ $/$ $ms$")
plt.ylabel(r"$U$ $/$ $V$")
plt.xscale("log")
#plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
#plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/T1.pdf")



#####RECHNUNGEN#######

###### T2-Auswertung #####
print("\n\t##### T2- Auswertung ######")
def fit_T1(t, a, b, c):
    return a*(1-2*np.exp(-t/b))+c

params1, cov1 = curve_fit(fit_T1, t_1, U_1, p0 = (-1.6, 2839, -0.638 ))
cov1 = np.sqrt(np.diag(abs(cov1)))

uparams1 = unp.uarray(params1, cov1)

print(f"U_0 = {noms(uparams1[0]):.4f} \pm {stds(uparams1[0]):.4f} Volt\nT1 = {noms(uparams1[1]):.4f} \pm {stds(uparams1[1]):.4f} ms\nU_1 = {noms(uparams1[2]):.4f} \pm {stds(uparams1[2]):.4f} Volt\n")



########Grafiken########



plt.figure()
plt.plot(t_2, U_2,"x",label="Messwerte")
#plt.xscale('log')
plt.xlabel(r"$\tau$ $/$ $ms$")
plt.ylabel(r"$U$ $/$ $V$")
#plt.xscale("log")
#plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
#plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/plot2.pdf")