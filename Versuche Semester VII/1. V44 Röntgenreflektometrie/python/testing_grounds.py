import numpy as np
import matplotlib.pyplot as plt
import uncertainties as unc
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
import scipy.constants as const
import sympy
import os
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from scipy.signal import find_peaks

if os.path.exists("build") == False:
    os.mkdir("build")

if os.path.exists("build/plots") == False:
    os.mkdir("build/plots")



########MESSWERTE#######

###Justage###
detectorscan, detectorscan_I = np.genfromtxt("python/data/detektorscan.UXD",skip_header= 57, unpack = True)
z1, z1_I = np.genfromtxt("python/data/z1.UXD",skip_header= 57, unpack = True)
x1_2theta, x1_I = np.genfromtxt("python/data/x1.UXD",skip_header= 57, unpack = True)
just_rock1_2theta, just_rock1_I = np.genfromtxt("python/data/just_rock1.UXD",skip_header= 57, unpack = True)
z2_2theta, z2_I = np.genfromtxt("python/data/z2.UXD",skip_header= 57, unpack = True)
just_rock2_2theta, just_rock2_I = np.genfromtxt("python/data/just_rock2.UXD",skip_header= 57, unpack = True)
z3_2theta, z3_I = np.genfromtxt("python/data/z3.UXD",skip_header= 57, unpack = True)



###Messwerte###



diffus_2theta, diffus_I = np.genfromtxt("python/data/diffus.UXD",skip_header= 56, unpack = True)
refl_2theta, refl_I = np.genfromtxt("python/data/reflektivitat.UXD",skip_header= 57, unpack = True)
just_rock3_2theta, just_rock3_I = np.genfromtxt("python/data/just_rock_03.UXD",skip_header= 57, unpack = True)
test_2theta, test_I = np.genfromtxt("python/data/test.UXD",skip_header= 57, unpack = True)

###Theoriewerte###

dicke = 20 # in mm
lam = 1.54e-10      #K-alpha Linie
n = 1 - 7.6e-6 + 1.54e-8j*141/(4*np.pi)
k = 2*np.pi / lam
print("\t\t k = ", k)

#Silizium

sil_rho = 20 * 10**(-6) # per metre**2
sil_del = 7.6 * 10**(-6)
sil_mu = 8600 #per metre
sil_alpha = 0.174 # degree

#Polyester

poly_rho = 9.5 * 10**(-6) # per metre**2
poly_del = 3.5 * 10**(-6)
poly_mu = 400 #per metre
poly_alpha = 0.153 # degree

#####RECHNUNGEN#######














###Funktionen###


def geometrie_corr(I, alpha, alpha_g, beam_width,D):
    I_corr = np.zeros(np.size(I))
    G = 1
    for i in range(np.size(I)):
        #print(alpha[i],"\t", alpha_g)
        if( np.abs(alpha_g) > alpha[i]):
            G = D * np.sin(alpha[i] * np.pi/180) / beam_width
        else:
            G = 1
        #print(G)
        I_corr[i] = I[i] / G
        I_corr[0] = I[0]
    return I_corr


def parrat_neu(a_i, delta2, delta3, sigma2, sigma3, z2, beta2, beta3):
    n2 = 1 - delta2 + 1j*beta2
    n3 = 1 - delta3 + 1j*beta3

    a_i = a_i * np.pi/180 #conversion to rad


    kz1 = k * np.sqrt(n1**2-np.cos(a_i)**2)
    kz2 = k * np.sqrt(n2**2-np.cos(a_i)**2)
    kz3 = k * np.sqrt(n3**2-np.cos(a_i)**2)


    r12 = (kz1 - kz2) / (kz1 + kz2) * np.exp(-2 * sigma2**2 * kz1 * kz2 )
    r23 = (kz2 - kz3) / (kz2 + kz3) * np.exp(-2 * sigma3**2 * kz2 * kz3 )


    x2 = np.exp(-2j * kz2 * z2) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)

    R_parr = np.abs(x1)**2

    return R_parr


loc_1 = 23
loc_2 = 30

beam_width = z1[loc_2] - z1[loc_1]

gauss_max = 838646.2462

alpha_g1 = just_rock1_2theta[6]
alpha_g2 = just_rock1_2theta[-4]
alpha_ges = 1/2 *(alpha_g1 - alpha_g2)

### Geometriefaktor###
print("\n\n---Geometriefaktor---")


refl_I_corr = geometrie_corr(refl_I, refl_2theta, alpha_ges, beam_width, dicke)
diffus_I_corr = geometrie_corr(diffus_I, diffus_2theta, alpha_ges, beam_width, dicke)


corr_data = refl_I_corr - diffus_I_corr








###Parrat###
print(f"\n\n---Paratt---")
print(f"Fit-Params direkt aus der V44.py nehmen")



#old params 

#n1 = 1 #Luft
#delta2 = 8*10**(-7)
#delta3 = 4*10**(-6)
#sigma1 = 3*10**(-9) # m
#sigma2 = 13*10**(-11) # m 
#z2 = 8.55*10**(-8) # m   #verändert die Frequenz

########################## 2 entspricjht poly und 3 sil
sil_del = 7.6 * 10**(-6)
poly_del = 3.5 * 10**(-6)

#better ones

#try 1
#n1 = 1 #Luft
#delta2 = 1.5*10**(-7)
#delta3 = 7.6*10**(-6)
#sigma2 = 0.7*10**(-6) # m
#sigma3 =  4.8*10**(-10) # m 
#z2 = 6.5*10**(-8) # m   #verändert die Frequenz
#beta2 = 7.5 * 10**(-7)
#beta3 = 35 * 10**(-8)

#better one
#n1 = 1 #Luft
#delta2 = 2*10**(-7) #poly
#delta3 = 7.5*10**(-6) # sil
#sigma2 = 0.8*10**(-10) # m // luft
#sigma3 = 5*10**(-10) # m  // sil #lloking good
#z2 = 8.65*10**(-8) # m   #verändert die Frequenz
#beta2 = 6.3 * 10**(-7)
#beta3 = 70 * 10**(-8)

n1 = 1 #Luft
delta2 = 7*10**(-7)
delta3 = 7*10**(-6)
sigma2 = 4*10**(-10) # m
sigma3 =  9*10**(-10) # m 
z2 = 6*10**(-8) # m   #verändert die Frequenz
beta2 = 8 * 10**(-7)
beta3 = 40 * 10**(-8)

n1 = 1 #Luft
delta2 = 4*10**(-7)
delta3 = 11*10**(-6)
sigma2 = 4*10**(-10) # m
sigma3 =  6*10**(-10) # m 
z2 = 8.55*10**(-8) # m   #verändert die Frequenz
beta2 = delta2/70 #3 * 10**(-7)
beta3 = delta3/20  #4 * 10**(-8)

#delta2 = poly_del
#delta3 = sil_del

delta2plus = delta2  #+ 0.1 * delta2
delta3plus = delta3  #+ 0.1 * delta2 
sigma1plus = sigma2  #+ 0.1 * sigma2
sigma2plus = sigma3  #+ 0.1 * sigma3
z2plus     = z2      #+ 0.1 * z2
beta2plus  = beta2   + 0.1 * beta2
beta3plus  = beta3   #+ 0.1 * beta3


delta2minus = delta2  #- 0.1 * delta2
delta3minus = delta3  #- 0.1 * delta2 
sigma1minus = sigma2  #- 0.1 * sigma2
sigma2minus = sigma3  #- 0.1 * sigma3
z2minus     = z2      #- 0.1 * z2
beta2minus  = beta2   - 0.1 * beta2
beta3minus  = beta3   #- 0.1 * beta3

#test mit hteo

#delta2 = poly_del
#delta3 = sil_del
#z2 = schichtdicke

#params_par, cov_par = curve_fit(parrat_rau, refl_2theta[thresh1:thresh2+ 60], corr_data[thresh1:thresh2+ 60]/gauss_max, p0 = (12*10**(-6), 10*10**(-6), 6*10**(-10), 12*10**(-10), 8.55*10**(-8) ))
#cov_par = np.sqrt(np.diag(cov_par))
#uparam2 = unp.uarray(params_par, cov_par)
#print(uparam2)

x = np.linspace(refl_2theta[0], refl_2theta[-1],10000)
plt.figure()
plt.plot(x, parrat_neu(x, delta2, delta3, sigma2, sigma3, z2, beta2, beta3), label = "Paratt-Fit (händisch)")
plt.plot(x, parrat_neu(x, delta2plus, delta3plus, sigma1plus, sigma2plus, z2plus, beta2plus, beta3plus), alpha = 0.5, label = "Paratt-Fit (+)")
plt.plot(x, parrat_neu(x, delta2minus, delta3minus, sigma1minus, sigma2minus, z2minus, beta2minus, beta3minus), alpha = 0.5, label = "Paratt-Fit (-)")
#plt.plot(refl_2theta[thresh1:thresh2 + 60], parrat_rau(refl_2theta[thresh1:thresh2 + 60],*params_par), label = "Paratt-Fit")
plt.plot(refl_2theta, corr_data/ gauss_max, label="korrigierte Daten")
plt.yscale('log')
plt.ylabel(r"Reflektivität")
plt.xlabel(r"$\alpha_{i}$ $/$ Grad")
#plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
#plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/paratt.pdf")

