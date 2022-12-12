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



########KOnstanten#######
n = 3.3543 #brechungsindex GaAs
B_maximum = 0.422 #tesla
L_rein = 5.11e-3 #metre
L_12 = 1.36e-3 #metre
L_28 = 1.296e-3 #metre
N_12 = 1.2e22 #\per\metre\cubic
N_28 = 2.8e22 #\per\metre\cubic

#####Funktionen#######

def fortable(eins, zwei, drei, vier):
    for i in range(len(eins)):
        print(f'{eins[i]:.3f} & {zwei[i]:.3f} & {drei[i]:.3f} & {vier[i]:.3f} \ \ \n ')
    return

def line(x, a, b):
    return a*x + b

def effmass(steig, N):
    konst_z = const.e**3 *10**12
    konst_n = 8*np.pi**2 * const.epsilon_0*const.c**3*n
    return unp.sqrt(konst_z/konst_n * N*B_maximum/steig)

########Grafiken########



print("______________________________________________________________B-FELD_________________________________________________________________________\n ")

dz, B = np.genfromtxt("python/data/bfeld.txt", unpack=True) # dz in mm, B in mT

B_max = np.max(abs(B)) #mT
print(f'maximales B-Feld: {B_max} mT\n')

plt.figure()
plt.plot(dz, abs(B), 'k.', label='Messwerte')
plt.hlines(B_max, -10, 10, color='red', linestyle='--', label='maximale Feldstärke')
# plt.xlabel(r'$z $/$ \si{\milli\metre}$')
# plt.ylabel(r'$B $/$ \si{\milli\tesla}$')
plt.legend()
plt.tight_layout()
plt.savefig('build/plots/bfeld.pdf')

print("____________________________________________________________REINES GaAs______________________________________________________________________\n")
def firstplot(name, lab, col):
    lamb, t1, wm1, t2, wm2 = np.genfromtxt('python/data/'+name+'.txt', unpack=True)

    theta_1 = (t1 + wm1/60)*(2*np.pi/360)   # in radiant
    theta_2 = (t2 + wm2/60)*(2*np.pi/360)   # in radiant

    theta = 0.5*(theta_1 - theta_2)

    plt.figure()
    plt.plot(lamb**2, theta, '.', color=col, label='Messwerte '+lab+'')
    plt.xlabel(r'$\lambda^2 $/$ \si{\micro\metre\squared}$')
    plt.ylabel(r'$\theta$')
    plt.xticks(lamb**2, ['$1.06^2$', '$1.29^2$', '$1.45^2$', '$1.72^2$', '$1.96^2$', '$2.156^2$', '$2.34^2$', '$2.51^2$', '$2.65^2$'])
    # plt.xlabel('lamb / microm^2')
    # plt.ylabel('theta')
    plt.legend()
    plt.tight_layout()
    plt.savefig('build/plots/'+name+'.pdf')

    print('___________________________'+name+'_____________'+lab+'_______________________')
    fortable(lamb, theta_1, theta_2, theta)
    return theta

theta_rein = firstplot('hochrein', 'reines GaAs', 'black')

theta_12 = firstplot('probe1', 'Probe 1', 'red')      #r'$N=\SI{1.2e18}{1\per\centi\metre\tothe{3}$',
theta_28 = firstplot('probe2', 'Probe 2', 'orange')   #r'$N=\SI{2.8e18}{1\per\centi\metre\tothe{3}$',


lamlam = np.array([1.06, 1.29, 1.45, 1.72, 1.96, 2.156, 2.34, 2.51, 2.65])


def secondplot(name, theta_probe, L_probe, N_probe, col):
    diff_theta = theta_probe/L_probe - theta_rein/L_rein

    params, cov = curve_fit(line, lamlam**2, diff_theta)
    errors = np.sqrt(np.diag(cov))

    a = ufloat(params[0], errors[0])
    b = ufloat(params[1], errors[1])

    l_plot = np.linspace(1, 2.65**2, 1000)

    print(f' _________________________{name}_______________________\n a = {a} \n b = {b} \n \n ')

    plt.figure()
    plt.plot(lamlam**2, diff_theta, '.', color=col, label='Messwerte')
    plt.plot(l_plot, line(l_plot, a.n, b.n), 'k--', label='Ausgleichsrechnung')
    plt.ylabel(r'$\theta_\text{frei} $/$ \si{\per\metre}$')
    plt.xlabel(r'$\lambda^2 $/$ \si{\micro\metre}$')
    plt.xticks(lamlam**2,  ['$1.06^2$', '$1.29^2$', '$1.45^2$', '$1.72^2$', '$1.96^2$', '$2.156^2$', '$2.34^2$', '$2.51^2$', '$2.65^2$'])
    plt.legend()
    plt.tight_layout()
    plt.savefig('build/plots/'+name+'_diff.pdf')

    print(N_probe)
    m_eff = effmass(abs(a), N_probe)

    print(f'effektive Masse m* = {m_eff} kg')
    return

secondplot('probe1', theta_12, L_12, N_12, 'red')
secondplot('probe2', theta_28, L_28, N_28, 'orange')