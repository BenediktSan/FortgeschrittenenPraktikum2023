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

def rel_abw(theo,a):
    c = (theo - a)/theo
    print(f"Relative Abweichung in Prozent: {noms(c) * 100 :.4f} \pm {stds(c) * 100 :.5f}\n")


def gauss(x, mu, sigma, A ):
    return A / ( np.sqrt(2 * np.pi * sigma) ) * np.exp( -(( x -mu )**2 / (2*sigma**2)) ) 

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


def dispersion(k, rho_e):
    r_0 = const.value("classical electron radius")
    return 2* np.pi * rho_e * r_0/k**2


#alter rauer parrat

#def parrat_rau(a_i,delta2,delta3,sigma1,sigma2,z2):
#    n2 = 1. - delta2
#    n3 = 1. - delta3
#
#    a_i = a_i * np.pi/180
#
#    kz1 = k * np.sqrt(np.abs(n1**2 - np.cos(a_i)**2))
#    kz2 = k * np.sqrt(np.abs(n2**2 - np.cos(a_i)**2))
#    kz3 = k * np.sqrt(np.abs(n3**2 - np.cos(a_i)**2))
#
#    r12 = (kz1 - kz2) / (kz1 + kz2) * np.exp(-2 * kz1 * kz2 * sigma1**2)
#    r23 = (kz2 - kz3) / (kz2 + kz3) * np.exp(-2 * kz2 * kz3 * sigma2**2)
#
#    x2 = np.exp(-2j * kz2 * z2) * r23
#    x1 = (r12 + x2) / (1 + r12 * x2)
#    R_parr = np.abs(x1)**2
#
#    return R_parr

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


#def theorie_sil(alpha):
#    return (np.abs((k * np.sin(alpha)- k*np.sqrt(n**2-np.cos(alpha)**2))/(k * np.sin(alpha)+ k*np.sqrt(n**2-np.cos(alpha)**2))))**2


def theorie_sil(alpha, alpha_crit, mu, lam):
    k = 2* np.pi /lam
    beta = mu/(2*k)
    R = ( alpha - np.sqrt( alpha**2 - alpha_crit**2 + 2 * 1j  * beta) ) / ( alpha + np.sqrt( alpha**2 - alpha_crit**2 + 2 * 1j  * beta) )
    R = R * np.conjugate(R)
    return R 

###Gauß-Fit###

param1, cov1 = curve_fit(gauss, detectorscan, detectorscan_I, p0 = (-0.0144, 0.04, 33.6))

cov1 = np.sqrt(np.diag(cov1))
uparam1 = unp.uarray(param1, cov1)

FWHM = 2 * np.sqrt(2 * np.log(2)) * uparam1[1]
gauss_max = gauss(param1[0], *param1)
print(f"\n---Gauss-Fit--- \nmu = {noms(uparam1[0]):.4f} \pm {stds(uparam1[0]):.4f} \nsigma = {noms(uparam1[1]):.4f} \pm {stds(uparam1[1]):.4f}")
print(f"A = {noms(uparam1[2]):.4f} \pm {stds(uparam1[2]):.4f} \n")
print(f"I_max Fit = {gauss_max:.4f}")
print(f"FWHM = {noms(FWHM):.4f} \pm {stds(FWHM):.4f}")

x1 = np.linspace(-0.5,0.5,1000)
FWHM_plot_1 = np.linspace(param1[0] - noms(FWHM)/2, param1[0] + noms(FWHM)/2, 100)
FWHM_plot_2 = np.ones((100)) * 1/2 * gauss_max

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
rel_abw(np.arcsin(((beam_width) / dicke))* 180/np.pi, np.abs(alpha_ges))

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


normierung_zeit = (refl_2theta[1] - refl_2theta[0]) * 10**3


refl_I_corr = refl_I / (normierung_zeit)
diffus_I_corr = diffus_I /(normierung_zeit)


refl_I_corr = geometrie_corr(refl_I_corr, refl_2theta, alpha_ges, beam_width, dicke)
diffus_I_corr = geometrie_corr(diffus_I_corr, diffus_2theta, alpha_ges, beam_width, dicke)



corr_data = refl_I_corr - diffus_I_corr





###Reflektivitäts & Diffuser Scan###
print("\n\n---Reflekitvitäts & Diffuser Scan---")

alpha_crit_theo = np.sqrt(2 * sil_del) * 180/np.pi # in Grad #nur hier fürs plotten der THeo. kommt später nochmal


plt.figure()
plt.plot(diffus_2theta, diffus_I/ gauss_max, label="normierter Diffuserscan")
plt.plot(refl_2theta , refl_I/ gauss_max, label="normierter Refklektivitätsscan")
plt.plot(refl_2theta, corr_data/ gauss_max, label="korrigierte Daten")
plt.plot(refl_2theta, theorie_sil(refl_2theta , alpha_crit_theo, sil_mu, lam) , label = "Theoriekurve Silizium")
plt.yscale('log')
plt.ylabel(r"Reflektivität")
plt.xlabel(r"$\alpha_{i}$ $/$ Grad")
#plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
#plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/refl1.pdf")









###Darstellung des relevanten Bereichs
print("\n\n---Peak Abstände & Schichtdicke---")

thresh1 = 40        #Grenzen um nur bestimmten BEreih zu betrachten (Nach Anfangsbums und vorm ausfasern der Daten)
thresh2 = 160

peaks, peak_heights = find_peaks(corr_data[thresh1:thresh2]/ gauss_max, height = 0.00044)
diff = np.zeros(np.size(peaks) -1)

for i in range(0,np.size(peaks) - 1):
    diff[i] = refl_2theta[peaks[i+1]] - refl_2theta[peaks[i]]   #ohne thresh1, wegen gleichbleibendem Abstand der Punkte

print(diff)
diff_mittel = np.sum(diff)/np.size(diff)
schichtdicke = lam/(2*diff_mittel * np.pi/180)
print(f"\ngemittelter Peak Abstand {diff_mittel} Grad\nSchichtdicke = {schichtdicke}")


print(f"\nWerte an den Grenzen. alpha \in  {refl_2theta[thresh1]} & {refl_2theta[thresh2]}")

plt.figure()
plt.plot(refl_2theta[thresh1:thresh2], corr_data[thresh1:thresh2]/ gauss_max, label="korrigierte Daten")
plt.plot(refl_2theta[peaks + thresh1], corr_data[peaks + thresh1]/ gauss_max,"rx", label = "Maxima")
plt.yscale('log')
plt.ylabel(r"Reflektivität")
plt.xlabel(r"$\alpha_{i}$ $/$ Grad")
#plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
#plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/refl2.pdf")







###Dispersion###
print(f"\n\n---Dispersion---")

print("???????")





###Kritischer Winkel###
print(f"\n\n---Kritscher Winkel---")

alpha_crit_theo = np.sqrt(2 * sil_del) * 180/np.pi # in Grad
print(f"\nTheorie kritischer Winkel = {alpha_crit_theo}")


alpha_crit, crit_height = find_peaks(refl_I[:thresh1 +10]/ gauss_max)
print(f"kritische Winkel = {refl_2theta[ alpha_crit]}")
rel_abw(alpha_crit_theo, refl_2theta[ alpha_crit[1]])



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

########################## 2 entspricht poly und 3 sil

n1 = 1 #Luft
delta2 = 10*10**(-7)
delta3 = 8.15*10**(-6)
sigma2 = 9*10**(-10) # m
sigma3 = 7.8*10**(-10) # m 
z2 = 8.8*10**(-8) # m   #verändert die Frequenz
beta2 = 3 * 10**(-10)
beta3 = delta3/50  #4 * 10**(-8)

print(f"beta_Si = {beta3}")

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
plt.plot(x, parrat_neu(x, delta2, delta3, sigma2, sigma3, z2, beta2, beta3),color ="teal", label = "Paratt-Fit (händisch)")
#plt.plot(refl_2theta[thresh1:thresh2 + 60], parrat_rau(refl_2theta[thresh1:thresh2 + 60],*params_par), label = "Paratt-Fit")
plt.plot(refl_2theta, corr_data/ gauss_max, color = "orange", label="korrigierte Daten")
plt.yscale('log')
plt.ylabel(r"Reflektivität")
plt.xlabel(r"$\alpha_{i}$ $/$ Grad")
#plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
#plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/paratt.pdf")





print(f"\n\n---Abweichungen---")

print("Vergleich dicken. zuerst bestimmt = theo")
rel_abw(schichtdicke, z2)

print("Delta. erst Poly dann Si")
rel_abw(poly_del, delta2)
print(sil_del,"\t", delta3)
rel_abw(sil_del, delta3)




#plt.figure()
#plt.plot(diffus_2theta, diffus_I/ gauss_max, label="normierter Diffuserscan")
#plt.plot(refl_2theta , refl_I/ gauss_max, label="normierter Refklektivitäätsscan")
#plt.plot(refl_2theta, corr_data/ gauss_max, label="korrigierte Daten")
#plt.plot(refl_2theta, theorie_sil(refl_2theta * np.pi /180), label = "Theoriekurve Silizium")
#plt.plot(refl_2theta[thresh1:], parrat_rau(refl_2theta[thresh1:], delta2, delta3, sigma1, sigma2, z2), label = "Paratt-Fit (händisch)")
#plt.yscale('log')
#plt.ylabel(r"Reflektivität")
#plt.xlabel(r"$\alpha_{i}$ $/$ Grad")
##plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
##plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
#plt.tight_layout()
#plt.legend()
#plt.savefig("build/plots/refl1.pdf")


