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

#___________________________________________________________GLOBALE KONSTANTEN__________________________________________________________________
c = const.c 
h = const.h
mu_B = const.physical_constants["Bohr magneton"][0]

L = 120e-3
nr = 1.4567 
nb = 1.4635
d = 4e-3
lr = 643.8e-9
lb = 480e-9

print(f'c = {c} \nh = {h} \nmu_B = {mu_B}\n')

#_______________________________________________________________HYSTERESE________________________________________________________________________

I, B = np.genfromtxt('python/data/hysterese.txt', unpack=True)

para1, cov1 = np.polyfit(I, B, 1, cov=True)
errors1 = np.sqrt(np.diag(cov1))
print(f'LINEAR \n a = {para1[0]} \pm {errors1[0]}, \n b = {para1[1]} \pm {errors1[1]} \n ')

para2, cov2 = np.polyfit(I, B, 2, cov=True)
errors2 = np.sqrt(np.diag(cov2))
print(f'QUADRATISCH \n a = {para2[0]} \pm {errors2[0]}, \n b = {para2[1]} \pm {errors2[1]} \n c = {para2[2]} \pm {errors2[2]}\n ')

para3, cov3 = np.polyfit(I, B, 3, cov=True)
errors3 = np.sqrt(np.diag(cov3))
print(f'KUBISCH\n a = {para3[0]} \pm {errors3[0]}, \n b = {para3[1]} \pm {errors3[1]} \n c = {para3[2]} \pm {errors3[2]}\n d = {para3[3]} \pm {errors3[3]}\n ')

x_plot = np.linspace(0, 8, 10000)

plt.figure()
plt.plot(I, B, 'k.', label='Messpunkte')
plt.plot(x_plot, para3[3] + x_plot*para3[2] + x_plot**2 *para3[1] + x_plot**3 *para3[0],'r--', label='Kubischer Fit')
plt.xlabel(r'$I \mathbin{/} \si{\ampere}$')
plt.ylabel(r'$B \mathbin{/} \si{\milli\tesla}$')
plt.legend()
plt.tight_layout()
plt.savefig('build/plots/hyst.pdf')

#___________________________________AUSWERTUNG_________________________________________________________________________________________________
def dispersion(lamb, d, n):
    ff = lamb**2/(2*d)
    fs = np.sqrt(n**2 - 1)
    return ff/fs

def Aufl(lamb, L, n):
    return L/lamb*(n**2 - 1)

def dell(ds, Ds, disp):
    fr = ds/Ds
    return 0.5*fr*disp

def gg(dlam, lamb, B):
    co = h*c/mu_B
    n = B*lamb**2
    return dlam*co/n
    
def mag(I, a0, a1, a2, a3):
    poly = a3*I**3 + a2*I**2 + a1*I +a0
    return poly*10**(-3)

def relab(ex,th):
    return 1 - ex/th

#________________________________ROT_______________________________________________________________________________
Ds , ds = np.genfromtxt('python/data/rot.txt', unpack=True)

#cut off
ds = ds[:-1]
ds_r = unp.uarray(ds, 10*np.ones(len(ds)))

Dsm_n = np.zeros(len(ds))

for i in range(len(ds)):
    Dsm_n[i] = 0.5*(Ds[i]+Ds[i+1])

Dsm = unp.uarray(Dsm_n, 10*np.ones(len(ds)))

print(f'DELTA s = {Dsm} \n \ndelta s = {ds_r}\n') 

#Umrechnen in del lamb 
disp_r = dispersion(lr, d, nr)
A = Aufl(lr, L, nr)

dsm = ufloat(np.mean(ds), 10)
Dsmm = ufloat(np.mean(Dsm_n), 10)
B_r = mag(8, para3[3], para3[2], para3[1], para3[0])

dell_r = dell(dsm, Dsmm, disp_r)

g_r = gg(dell_r, lr, B_r)

print(f' dispersionsgebiet {disp_r}\n Auflösung {A}\n\n delta s = {dsm}\n DELTA s = {Dsmm}\n B = {B_r}\n\n delta l = {dell_r}\n\n g = {g_r}\n abw {relab(g_r,1)}' )


#_________________________________________________BLAU SIGMA___________________________________________________________________________________

Dbs , dbs =  np.genfromtxt('python/data/blausigma.txt', unpack=True)

#cut off
dbs = dbs[:-1]
ds_b = unp.uarray(dbs, 10*np.ones(len(dbs)))

Dbs_n = np.zeros(len(dbs))

for i in range(len(dbs)):
    Dbs_n[i] = 0.5*(Dbs[i]+Dbs[i+1])

Dbsm = unp.uarray(Dbs_n, 10*np.ones(len(dbs)))

print(f'DELTA s = {Dbsm} \n \ndelta s = {ds_b}\n') 

#Umrechnen in del lamb 
disp_b = dispersion(lb, d, nb)
A = Aufl(lb, L, nb)

dbsm = ufloat(np.mean(dbs), 10)
Dbsmm = ufloat(np.mean(Dbs_n), 10)
B_b = mag(3.4, para3[3], para3[2], para3[1], para3[0])

dell_b = dell(dbsm, Dbsmm, disp_b)

g_bs = gg(dell_b, lb, B_b)

print(f' dispersionsgebiet {disp_b}\n Auflösung {A}\n\n delta s = {dbsm}\n DELTA s = {Dbsmm}\n B = {B_b}\n\n delta l = {dell_b}\n\n g = {g_bs}\n abw {relab(g_bs,1.75)}' )

#____________________________________________________________BLAU PI_____________________________________________________________________________

Dbp , dbp = np.genfromtxt('python/data/blaupi.txt', unpack=True)

#cut off
dbp = dbp[:-1]
ds_bp = unp.uarray(dbp, 10*np.ones(len(dbp)))

Dbp_n = np.zeros(len(dbp))

for i in range(len(dbp)):
    Dbp_n[i] = 0.5*(Dbp[i]+Dbp[i+1])

Dbpm = unp.uarray(Dbp_n, 10*np.ones(len(dbp)))

print(f'DELTA s = {Dbpm} \n \ndelta s = {ds_bp}\n') 

#Umrechnen in del lamb 

dbpm = ufloat(np.mean(dbp), 10)
Dbpmm = ufloat(np.mean(Dbp_n), 10)
B_bp = mag(8, para3[3], para3[2], para3[1], para3[0])

dell_bp = dell(dbpm, Dbpmm, disp_b)

g_bp = gg(dell_bp, lb, B_bp)

print(f' delta s = {dbpm}\n DELTA s = {Dbpmm}\n B = {B_bp}\n\n delta l = {dell_bp}\n\n g = {g_bp}\n abw {relab(g_bp,0.5)}' )