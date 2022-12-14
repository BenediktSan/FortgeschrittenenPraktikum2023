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

def fortable(eins, zwei, drei, vier, fünf):
    for i in range(len(eins)):
        print(f'{eins[i]:.3f} & {zwei[i]:.3f} & {drei[i]:.3f} & {vier[i]:.3f} & {fünf[i]:.3f} \ \ \n ')
    return

def line(x, a, b):
    return a*x + b

def effmass(steig, N):
    konst_z = const.e**3 
    konst_n = 8*np.pi**2 * const.epsilon_0*const.c**3*n
    return unp.sqrt(konst_z/konst_n * N*B_maximum/steig)

def relab(exp, theo):
    return 1 - exp/theo

########Grafiken########



print("______________________________________________________________B-FELD_________________________________________________________________________\n ")

dz, B = np.genfromtxt("python/data/bfeld.txt", unpack=True) # dz in mm, B in mT

B_max = np.max(abs(B)) #mT
print(f'maximales B-Feld: {B_max} mT\n')

# plt.figure()
# plt.plot(dz, abs(B), 'k.', label='Messwerte')
# plt.hlines(B_max, -10, 10, color='red', linestyle='--', label='maximale Feldstärke')
# plt.xlabel(r'$z $/$ \si{\milli\metre}$')
# plt.ylabel(r'$B $/$ \si{\milli\tesla}$')
# plt.legend()
# plt.tight_layout()
# plt.savefig('build/plots/bfeld.pdf')

fortable(dz, B, np.zeros(len(B)), np.zeros(len(B)), np.zeros(len(B)))

print("____________________________________________________________REINES GaAs______________________________________________________________________\n")
def firstplot(name, lab, col, L_probe):
    lamb, t1, wm1, t2, wm2 = np.genfromtxt('python/data/'+name+'.txt', unpack=True)

    theta_1 = (t1 + wm1/60)*(2*np.pi/360)   # in radiant
    theta_2 = (t2 + wm2/60)*(2*np.pi/360)   # in radiant

    theta = 0.5*(theta_1 - theta_2)

    # plt.figure()
    # plt.plot(lamb**2, theta, '.', color=col, label='Messwerte '+lab+'')
    # plt.xlabel(r'$\lambda^2 $/$ \si{\micro\metre\squared}$')
    # plt.ylabel(r'$\theta$')
    # plt.xticks(lamb**2, ['$1.06^2$', '$1.29^2$', '$1.45^2$', '$1.72^2$', '$1.96^2$', '$2.156^2$', '$2.34^2$', '$2.51^2$', '$2.56^2$'])
    # # plt.xlabel('lamb / microm^2')
    # # plt.ylabel('theta')
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig('build/plots/'+name+'.pdf')

    print('___________________________'+name+'_____________'+lab+'_______________________')
    fortable(lamb, theta_1, theta_2, theta, theta/L_probe)
    return theta

theta_rein = firstplot('hochrein', 'reines GaAs', 'black', L_rein)

theta_12 = firstplot('probe1', 'Probe 1', 'red', L_12)      #r'$N=\SI{1.2e18}{1\per\centi\metre\tothe{3}$',
theta_28 = firstplot('probe2', 'Probe 2', 'orange', L_28)   #r'$N=\SI{2.8e18}{1\per\centi\metre\tothe{3}$',


lamlam = np.array([1.06, 1.29, 1.45, 1.72, 1.96, 2.156, 2.34, 2.51, 2.56])

# plt.figure()
# plt.plot(lamlam**2, theta_rein/L_rein, '.', color='black', label=r'$\che{GaAs}$')
# plt.plot(lamlam**2, theta_12/L_12, '.', color='red', label=r'$\che{InGaAs}, N = \SI{1.2e18}{\per\centi\metre\tothe{3}}$ ')
# plt.plot(lamlam**2, theta_28/L_28, '.', color='orange', label=r'$\che{InGaAs}, N = \SI{2.8e18}{\per\centi\metre\tothe{3}$')
# plt.xlabel(r'$\lambda^2 $/$ \si{\micro\metre\squared}$')
# plt.ylabel(r'$\theta_\text{frei} $/$ \si{1\per\metre}$')
# plt.xticks(lamlam**2, ['$1.06^2$', '$1.29^2$', '$1.45^2$', '$1.72^2$', '$1.96^2$', '$2.156^2$', '$2.34^2$', '$2.51^2$', '$2.56^2$'])
# # plt.xlabel('lamb / microm^2')
# # plt.ylabel('theta')
# plt.legend()
# plt.tight_layout()
# plt.savefig('build/plots/firstplot.pdf')


def secondplot(name, theta_probe, L_probe, N_probe, col):
    diff_theta = theta_probe/L_probe - theta_rein/L_rein

    params, cov = curve_fit(line, lamlam**2, diff_theta)
    errors = np.sqrt(np.diag(cov))

    a = ufloat(params[0], errors[0])
    b = ufloat(params[1], errors[1])

    l_plot = np.linspace(1, 2.65**2, 1000)

    print(f' _________________________{name}_______________________\n a = {a} \n b = {b} \n \n ')

    # plt.figure()
    # plt.plot(lamlam**2, diff_theta, '.', color=col, label='Messwerte')
    # plt.plot(l_plot, line(l_plot, a.n, b.n), 'k--', label='Ausgleichsrechnung')
    # plt.ylabel(r'$\theta_\text{frei} $/$ \si{\per\metre}$')
    # plt.xlabel(r'$\lambda^2 $/$ \si{\micro\metre}$')
    # plt.xticks(lamlam**2,  ['$1.06^2$', '$1.29^2$', '$1.45^2$', '$1.72^2$', '$1.96^2$', '$2.156^2$', '$2.34^2$', '$2.51^2$', '$2.56^2$'])
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig('build/plots/'+name+'_diff.pdf')

    print(N_probe)
    m_eff = effmass(abs(a)*10**12, N_probe)

    print(f'effektive Masse m* = {m_eff} kg = {m_eff/const.m_e} m_e \n abweichung {relab(m_eff/const.m_e, 0.067)}')
    return

secondplot('probe1', theta_12, L_12, N_12, 'red')
secondplot('probe2', theta_28, L_28, N_28, 'orange')