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

#Silizium

sil_rho = 20 * 10**(-6) # per metre**2
sil_del = 7.6 * 10**(-6)
sil_mu = 8600 #per metre
sil_alp = 0.174 # degree

#Polyester

poly_rho = 9.5 * 10**(-6) # per metre**2
poly_del = 3.5 * 10**(-6)
poly_mu = 400 #per metre
poly_alp = 0.153 # degree

#####RECHNUNGEN#######














###Funktionen###

def rel_abw(theo,a):
    c = (theo - a)/theo
    print(f"Relative Abweichung in Prozent: {noms(c) * 100 :.4f} \pm {stds(c) * 100 :.5f}\n")


def gauss(x, mu, sigma, A, B):
    return A / ( np.sqrt(2 * np.pi * sigma) ) * np.exp( -(( x -mu )**2 / (2*sigma**2)) + B ) 

def geometrie_corr(I, alpha, alpha_g, beam_width):
    I_corr = np.zeros(np.size(I))
    G = 1
    for i in range(np.size(I)):
        #print(alpha[i],"\t", alpha_g)
        if( np.abs(alpha_g) > alpha[i]):
            G = np.sin(alpha[i] * np.pi/180) / beam_width
        else:
            G = 1
        #print(G)
        I_corr[i] = I[i] * G
    return I_corr

def parrat():
    
    return











###Gauß-Fit###

param1, cov1 = curve_fit(gauss, detectorscan, detectorscan_I)

cov1 = np.sqrt(np.diag(cov1))
uparam1 = unp.uarray(param1, cov1)

FWHM = 2 * np.sqrt(2 * np.log(2)) * uparam1[1]
print(f"\n---Gauss-Fit--- \nmu = {noms(uparam1[0]):.4f} \pm {stds(uparam1[0]):.4f} \nsigma = {noms(uparam1[1]):.4f} \pm {stds(uparam1[1]):.4f}")
print(f"A = {noms(uparam1[2]):.4f} \pm {stds(uparam1[2]):.4f} \nB = {noms(uparam1[3]):.4f} \pm {stds(uparam1[3]):.4f}  \n")
print(f"I_max Fit = {gauss(param1[0] , *param1):.4f}")
print(f"FWHM = {noms(FWHM):.4f} \pm {stds(FWHM):.4f}")

x1 = np.linspace(-0.5,0.5,1000)
FWHM_plot_1 = np.linspace(param1[0] - noms(FWHM)/2, param1[0] + noms(FWHM)/2, 100)
FWHM_plot_2 = np.ones((100)) * 1/2 * gauss(param1[0] , *param1)

plt.figure()
plt.plot(detectorscan, detectorscan_I,"x",label="Messwerte")
plt.plot(x1,gauss(x1,*param1), label="Gauss-Fit" )
plt.plot(FWHM_plot_1, FWHM_plot_2, label="Halbwertsbreite")
#plt.xscale('log')
plt.ylabel(r"Intensität")
plt.xlabel(r"$\alpha_{i}$ $/$ Grad")
#plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
#plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/gauss.pdf")






###Erster Z-Scan (Strahlbreite)###
print("\n\n---Erster Z-Scan (Strahlbreite)---")

loc_1 = 23
loc_2 = 30

beam_width = z1[loc_2] - z1[loc_1]

print(f"Beginn & Ende = {z1[loc_1]} & {z1[loc_2]} in mm\nStrahlbreite = {beam_width} mm")
print(f"Geometriewinkel über breite = {np.arcsin(((beam_width) / dicke) ) * 180/np.pi } in Grad")
x2 = np.linspace(z1[loc_1], z1[loc_2],100)
height_1 = np.ones(100) * z1_I[loc_1] / 2

plt.figure()
plt.plot(z1, z1_I,"x",label="Messwerte")
plt.bar(z1[loc_1], z1_I[loc_1] + 20000, width=0.01, bottom=-10000, color = "orange")
plt.bar(z1[loc_2], z1_I[loc_1] + 20000, width=0.01, bottom=-10000, color = "orange")
plt.plot(x2, height_1, color = "orange", label = "Strahlbreite")
#plt.yscale('log')
plt.ylabel(r"Intensität")
plt.xlabel(r"h $/$ mm")
#plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
#plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/z1_scan.pdf")







###Rocking-Scan1 (Dreieck)###
print("\n\n---Rocking-Scan1 (Dreieck)---")

alpha_g1 = just_rock1_2theta[6]
alpha_g2 = just_rock1_2theta[-4]
alpha_ges = 1/2 *(alpha_g1 - alpha_g2)

print(f"Geometriewinkel = {alpha_g1} & {alpha_g2} in Grad\nGemittelt = {alpha_ges}")

plt.figure()
plt.plot(just_rock1_2theta, just_rock1_I,"x",label="Messwerte")
plt.plot((alpha_g1, alpha_g2), (just_rock1_I[6],just_rock1_I[-4]),"rx",label="Geometriewinkel")
#plt.yscale('log')
plt.ylabel(r"Intensität")
plt.xlabel(r"$\alpha_{i}$ $/$ Grad")
#plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
#plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/dreieck.pdf")







### Geometriefaktor###
print("\n\n---Geometriefaktor---")


refl_I_corr = geometrie_corr(refl_I, refl_2theta, alpha_ges, beam_width)






###Reflektivitäts & Diffuser Scan###
print("\n\n---Reflekitvitäts & Diffuser Scan---")

plt.figure()
plt.plot(diffus_2theta, diffus_I,"x",label="Diffuserscan")
plt.plot(refl_2theta, refl_I, "x", label="Refklektivitäätsscan")
plt.plot(refl_2theta, refl_I - diffus_I, "x", label="korrigierte Daten")
#plt.plot(refl_2theta, refl_I_corr, "x", label="Reflektivitätsscan korrigiert")
#plt.yscale('log')
plt.ylabel(r"Intensität")
plt.xlabel(r"$\alpha_{i}$ $/$ Grad")
#plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
#plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/refl.pdf")




