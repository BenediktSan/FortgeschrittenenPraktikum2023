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
from tabulate import tabulate
from scipy.signal import find_peaks



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
U_diff = np.array([1.5, 1.487, 1.487, 1.468, 1.375, 1.15, 0.837, 0.475, 0.212, 0.160, 0.110, 0.073, 0.047, 0.041, 0.033, 0.029, 0.024, 0.02, 0.027]) # in V

tau5_diff, Real_diff, Imag_diff = np.genfromtxt("python/data/echo_gradient.csv",skip_header= 5, unpack = True)


#####FUNKTIONEN#######

def rel_abw(theo,a):
    c = (theo - a)/theo
    print(f"Relative Abweichung in Prozent: {noms(c) * 100 :.4f} \pm {stds(c) * 100 :.5f}\n")


def printer(a,b,c,d,name):
    print(f"### {name} ###")
    table1 ={'Messreihe 1': a, 'Messreihe 2': b,  'Messreihe 3': c, 'Messreihe 4': d }
    print("\n", tabulate(table1, tablefmt = "latex_raw"))  

def printer(a,b,name):
    print(f"### {name} ###")
    table1 ={'Messreihe 1': a, 'Messreihe 2': b }
    print("\n", tabulate(table1, tablefmt = "latex_raw"))  


#####RECHNUNGEN#######








###### T1-Auswertung #####
print("\n\t##### T1- Auswertung ######")

def fit_exp(t, a, b, c):
    return a*(1-2*np.exp(-t/b))+c

params1, cov1 = curve_fit(fit_exp, t_1, U_1, p0 = (-1.6, 2839, -0.638 ))
cov1 = np.sqrt(np.diag(abs(cov1)))

uparams1 = unp.uarray(params1, cov1)

print(f"U_0 = {noms(uparams1[0]):.4f} \pm {stds(uparams1[0]):.4f} Volt\nT1 = {noms(uparams1[1]):.4f} \pm {stds(uparams1[1]):.4f} ms\nU_1 = {noms(uparams1[2]):.4f} \pm {stds(uparams1[2]):.4f} Volt\n")


x1 = np.linspace(t_1[0], t_1[-1], 1000)

plt.figure()
plt.plot(t_1, U_1,"x",label="Messwerte")
plt.plot(x1, fit_exp(x1, *params1), label = "Fit-Funktion")
#plt.xscale('log')
plt.xlabel(r"$\tau$ $/$ $ms$")
plt.ylabel(r"$U$ $/$ $V$")
plt.xscale("log")
#plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
#plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/T1.pdf")









###### T2-Auswertung #####
print("\n\t##### T2- Auswertung ######")


U_2_peaks, peak_heights = find_peaks(U_2, height = 0.25)
t_2_new = t_2[U_2_peaks]
U_2_new = U_2[U_2_peaks]

params2, cov2 = curve_fit(fit_exp,2 * t_2_new, U_2_new, p0 = (-0.44, 1.97, 0.4 ))
cov2 = np.sqrt(np.diag(abs(cov2)))
uparams2 = unp.uarray(params2, cov2)

print(f"U_0 = {noms(uparams2[0]):.4f} \pm {stds(uparams2[0]):.4f} Volt\nT2 = {noms(uparams2[1]):.4f} \pm {stds(uparams2[1]):.4f} s\nU_1 = {noms(uparams2[2]):.4f} \pm {stds(uparams2[2]):.4f} Volt\n")

x2 = np.linspace(t_2_new[0], t_2_new[-1])

plt.figure()
#plt.plot(t_2, U_2,".",label="Messwerte")
plt.plot(t_2_new, U_2_new, "x", label = "Peaks")
plt.plot(x2, fit_exp(2 * x2, *params2), label = "Fit-Funktion")
#plt.xscale('log')
plt.xlabel(r"$\tau$ $/$ $s$")
plt.ylabel(r"$U$ $/$ $V$")
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/T2.pdf")






#######DIFFUSIONSKONSTANTE#########
print("\n\t##### Dffusionskonstante ######")


G = 0.1203  #aus der fourier.py
gyro = 2.67*10**8 #für Protonen T/s
global T2 
T2 = uparams2[1] * 10**(3) # in ms

def fit_diff(t, a, b, c):
    #T2 = 3.9523
    T_2 = 3.9523 * 10**(3)
    return a * np.exp(-(2*t)/T_2) * np.exp(-t**3/b)+c

params3, cov3 = curve_fit(fit_diff, t_diff, U_diff, p0 = (1.46, 1677, 0.026))
cov3 = np.sqrt(np.diag(abs(cov3)))
uparams3 = unp.uarray(params3, cov3)

D = 3 / (2 * uparams3[1] * 10**(-9) * gyro**2 * G**2) * 10**9

print(f"U_0 = {noms(uparams3[0]):.4f} \pm {stds(uparams3[0]):.4f} Volt\nT3 = {noms(uparams3[1]):.4f} \pm {stds(uparams3[1]):.4f} (ms)^3\nU_1 = {noms(uparams3[2]):.4f} \pm {stds(uparams3[2]):.4f} Volt\n")
print(f"D = {noms(D):.4f} \pm {stds(D):.4f}")
x3 = np.linspace(t_diff[0], t_diff[-1])

plt.figure()
#plt.plot(t_2, U_2,".",label="Messwerte")
plt.plot(t_diff, U_diff, "x", label = "Messwerte")
plt.plot(x3, fit_diff(x3, *params3), label = "Fit-Funktion")
#plt.xscale('log')
plt.xlabel(r"$\tau$ $/$ $ms$")
plt.ylabel(r"$U$ $/$ $V$")
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/diff1.pdf")


















########Grafiken########




######MESSDATEN##########


#printer(t_2[U_2_peaks], U_2[U_2_peaks] ,"peaks")
#printer(t_diff, U_diff ,"Diffusion")



