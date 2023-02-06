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

###invertierter Linearverstärker #1
#Spannungen peak to peak

lin1_R1 = 1000 #ohm
lin1_R2 = 10000 #ohm

lin1_f = np.array([15, 30.5, 50, 80.5, 100, 150, 250, 500, 750, 1000, 1500, 2000, 3500, 5000, 7500, 10000, 
                   15000, 20000, 30000, 40000, 50000, 100000, 250000, 500000, 1000000]) #Hz

lin1_Ue = np.array([205, 209, 205, 205, 205, 209, 205, 209, 209, 209, 209, 209, 209, 
                    209, 209, 209, 209, 209, 213, 217, 217, 215, 215, 213, 209]) * 10**(-3)#Volt

lin1_Ua = np.array([2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14,
                    2.10, 2.08, 2, 1.68, 1.33, 0.98, 0.6, 0.35, 0.24, 0.24]) #Volt

lin1_phase = np.array([180, 179, 178, 180, 180, 180, 179, 180, 181, 180, 176, 175, 170, 165, 160, 150, 140, 
                       100, 100, 88, 70, 47, np.NaN, np.NaN, np.NaN]) #rad


print(np.size(lin1_f), np.size(lin1_Ue), np.size(lin1_Ua), np.size(lin1_phase))

###invertierter Linearverstärker #2
#Spannungen peak to peak

lin2_R1 = 1000 #ohm
lin2_R2 = 100000 #ohm

lin2_f = np.array([15, 30, 50, 100, 200, 500, 750, 1000, 1500, 2500, 5000, 7500, 10000, 
            20000, 30000, 50000, 100000, 250000, 500000, 1000000]) #Hertz

lin2_Ue = np.array([205, 205, 205, 205, 205, 205, 205, 205, 205, 205, 209, 209, 
                    211, 215, 215, 215, 213, 213, 213, 209]) * 10**(-3) #Volt

lin2_Ua = np.array([19.8, 19.8, 19.8, 19.8, 19.8, 19.8, 19.7, 19.7, 19.6, 19, 
                    17, 14, 11.1, 5.2, 3.2, 1.37, 0.64, 0.27, 0.17, 0.14])  #Volt

lin2_phase = np.array([180, 179, 179, 179, 179, 184, 174, 171, 168, 161, 140, 126,
                       115, 80, 70, 55, 40, 30, np.NaN, np.NaN]) #rad

print(np.size(lin2_f), np.size(lin2_Ue), np.size(lin2_Ua), np.size(lin2_phase))


###invertierter Linearverstärker #3
#Spannungen peak to peak

lin3_R1 = 1000 #ohm
lin3_R2 = 68000#ohm

lin3_f = np.array([15, 30, 50, 100, 200, 500, 750, 1000, 1500, 2500, 5000, 7500, 10000, 
            20000, 30000, 50000, 100000, 250000, 500000, 1000000]) #Hertz

lin3_Ue = np.array([205, 205, 205, 205, 205, 205, 205, 205, 205, 205, 205, 209, 
                    209, 213, 213, 215, 213, 213, 213, 207]) * 10**(-3) #Volt
                    
lin3_Ua = np.array([13.9, 14, 14, 14, 14, 14, 13.9, 13.9, 13.8, 13.6, 13, 11.8,
                    10.2, 5.5, 3.3, 1.6, 1.1, 0.27, 0.155, 0.12])  #Volt

lin3_phase = np.array([179, 179, 179, 179, 179, 176, 176, 174, 171, 166,
                       152, 139, 126, 95, 75, 47, 43, 26, 15, np.NaN]) #rad

print(np.size(lin3_f), np.size(lin3_Ue), np.size(lin3_Ua), np.size(lin3_phase))


###Umkehr Integrator

int_R1 = 10000 #Ohm
int_C = 100 * 10**(-9) #farad

int_f = np.array([5, 10, 20, 35, 50, 75, 100, 125, 150, 175, 200, 225, 250, 300, 500]) #Hertz

int_Ue = np.array([1.07, 1.09, 1.08, 1.09, 1.08, 1.08, 1.08, 1.09, 1.09, 1.09, 1.09, 1.09, 1.09, 1.09, 1.09 ])  #Volt
                    
int_Ua = np.array([23.5, 13.4, 7.6, 4.9, 3.3, 2.4, 1.9, 1.6, 1.4, 1.2, 1.1, 1.1, 0.9, 0.8, 0.7])  #Volt

int_phase = np.array([96, 94, 93, 92, 90, 89, 89, 91, 89, 90, 87, 83, 90, 80, 80]) #rad


print(np.size(int_f), np.size(int_Ue), np.size(int_Ua), np.size(int_phase))


###invertierter Differenzierer

diff_R1 = 100000 #Ohm
diff_C = 22 * 10**(-9) #farad

diff_f = np.array([20, 40, 60, 100, 125, 175, 350, 500, 750, 1000, 1500]) #Hertz

diff_Ue = np.array([2.01, 2.01, 2.01, 2.05, 2.05, 2.05, 2.09, 2.09, 2.09, 2.09, 2.09 ])  #Volt
                    
diff_Ua = np.array([0.54, 0.96, 1.37, 2.2, 2.7, 3.74, 7.8, 11.3, 16.7, 25.3, 28.3])  #Volt

diff_phase = np.array([90, 90, 90, 91, 91, 90, 90, 90, 90, 90, 150 ]) #rad


print(np.size(diff_f), np.size(diff_Ue), np.size(diff_Ua), np.size(diff_phase))




###nicht invertierter Schmitt Trigger

schmitt_R1 = 10000 #Ohm
schmitt_R2 = 100000 #Ohm

schmitt_Ue = 14 #Volt # max spannung ds schmitt

#Kipppunkte

schmitt_x = np.array([0.396, 1.71, -1.04, -2.35]) * 10**(-3) #second

schmitt_y = np.array([1.86250, -1.6625, -1.675, 1.8625])    #volt


###Signalverstärker

signal_R1 = 10000 #Ohm
signal_R2 = 100000 #Ohm
signal_R3 = 1000
signal_C = 1e-6 #farad




#####Funktionen######


def lin(x,m,n):
    return m*x + n

def hyperbel(x, A):
    return 1/(A*x) 

def weirdfunc(x, A, B):
    if B < 0:
        return (x)**B * 1/A 
    elif B > 0:
        return (x)**B * A

    return 0 

def rel_abw(theo,a):
    c = (theo - a)/theo
    print(f"Relative Abweichung in Prozent: {noms(c) * 100 :.4f} \pm {stds(c) * 100 :.5f}\n")


def printer4(a,b,c,d,name):
    print(f"### {name} ###")
    table1 ={'Messreihe 1': a, 'Messreihe 2': b,  'Messreihe 3': c, 'Messreihe 4': d }
    print("\n", tabulate(table1, tablefmt = "latex_raw"))  

def printer2(a,b,name):
    print(f"### {name} ###")
    table1 ={'Messreihe 1': a, 'Messreihe 2': b }
    print("\n", tabulate(table1, tablefmt = "latex_raw"))  

 

#####RECHNUNGEN#######


size_label = 15

plt.rc('axes', labelsize=size_label)


#### Invertierter Linearverstärker ####
print("\n\n\t#### Invertierter Linearverstärker ####\n")


def inv_lin(R1, R2, f, Ue, Ua, thresh, name):

    print(f"#### {name} ####")



    V_theo = R2/R1
    print(f"V_theo = {V_theo}")

    #Plateaumittel

    mittel = np.mean( (Ua/Ue)[:thresh])
    std = np.std((Ua/Ue)[:thresh])
    umittel = ufloat(mittel, std)       #Leerlaufverstärkung
    print(f"\n\tLeerlaufverstärkung übers Plateau:\nmean = {noms(umittel):.4f} \pm {stds(umittel):.4f} s\n")

    #Flankenfit
    
    params1, cov1 = curve_fit(lin, np.log(f[thresh:]), np.log((Ua/Ue)[thresh:]))
    cov1 = np.sqrt(np.diag(abs(cov1)))
    uparams1 = unp.uarray(params1, cov1)

    print(f"\tFlanke:\nm = {noms(uparams1[0]):.4f} \pm {stds(uparams1[0]):.4f} s\nn = {noms(uparams1[1]):.4f} \pm {stds(uparams1[1]):.4f}\n")
    rel_abw(V_theo, mittel)

    #other stuff


    f_grenz = unp.exp((unp.log((umittel)/ np.sqrt(2))  - uparams1[1] ) / uparams1[0])
    band_prod = f_grenz * umittel

    print(f"\tother stuff:\nf_grenz = {noms(f_grenz):.4f} \pm {stds(f_grenz):.4f} 1/s")
    print(f"Bandbreitenprod = {noms(band_prod):.4f} \pm {stds(band_prod):.4f} 1/s\n")



    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.plot(f[:thresh], (Ua/Ue)[:thresh],"x",label="Messwerte für Mittelwertberechnung")
    plt.plot(f[thresh:], (Ua/Ue)[thresh:],"x",label="Messwerte")
    plt.plot(f[(thresh-2):-1], np.exp(lin(np.log(f[(thresh-2):-1]),*params1)), label="Flankenfit")
    plt.axhline(y = noms(mittel), xmax = np.log(f[thresh + 2]) / np.log(f[-1]), label = "Mittelwert")
    plt.axhline(y = noms(mittel) / np.sqrt(2),xmax = np.log(f[thresh + 4]), color = "red", label = r"V / $\sqrt{2}$")
    #plt.axvline(x = noms(f_grenz),color = "red")
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r"$V=\frac{U_e}{U_a}$")
    plt.xlabel(r"$f  $ / $ Hz$")
    #plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
    #plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
    plt.tight_layout()
    plt.legend()
    plt.savefig("build/plots/" + name + ".pdf")
    
    return


#invertierter Linearverstärker

#threshholds über ablesen

thresh1 = 17
thresh2 = 10
thresh3 = 11


inv_lin(lin1_R1, lin1_R2, lin1_f, lin1_Ue, lin1_Ua, thresh1, 'lin1')
inv_lin(lin2_R1, lin2_R2, lin2_f, lin2_Ue, lin2_Ua, thresh2, 'lin2')
inv_lin(lin3_R1, lin3_R2, lin3_f, lin3_Ue, lin3_Ua, thresh3, 'lin3')

plt.figure()
plt.rc('axes', labelsize=size_label)
plt.plot(lin1_f, lin1_phase,"x",label=r"Phasenverschiebung V $\approx$ 10 ")
plt.plot(lin2_f, lin2_phase,"x",label=r"Phasenverschiebung V $\approx$ 100 ")
plt.plot(lin3_f, lin3_phase,"x",label=r"Phasenverschiebung V $\approx$ 67 ")
plt.xscale('log')
plt.ylabel(r"$Phase$ / $rad$")
plt.xlabel(r"$f  $ / $ Hz$")
#plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
#plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/lin_phase.pdf")















#### Umkehr Integrator ####
print("\n\n\n\t#### Umkehr Integrator & Differenzierer ####\n")

def intdiff(R, C, f, Ue, Ua, name):

    print(f"#### {name} ####")


    #Umrechnung von f zu omega
    omega = 2 * np.pi * f 

    RC_theo = R * C


    #Flankenfit
    
    params1, cov1 = curve_fit(lin, np.log(omega), np.log(Ua/Ue))
    cov1 = np.sqrt(np.diag(abs(cov1)))
    uparams1 = unp.uarray(params1, cov1)



    
    print(f"\nm = {noms(uparams1[0]):.4f} \pm {stds(uparams1[0]):.4f} s\nn = {noms(uparams1[1]):.4f} \pm {stds(uparams1[1]):.4f}")
    m_theo = 0

    if name == "integrator":
        m_theo = -1
        RC = 1/unp.exp(uparams1[1])
    elif name == "differenzierer":
        m_theo = 1
        RC = unp.exp(uparams1[1])
    else:
        print("\nFalscher name")
    
    print(f"RC_theo = {RC_theo:.4f} s")
    print(f"RC = {noms(RC):.4f} \pm {stds(RC):.4f} s")
    rel_abw(RC_theo, RC)

    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.plot(omega, (Ua/Ue), "x", label="Messwerte")
    #plt.plot(omega, weirdfunc(omega, noms(RC), noms(uparams1[0])), label = "RC-test")
    plt.plot(omega, weirdfunc(omega, noms(RC_theo), m_theo), label = "Theorierkurve")
    plt.plot(omega, np.exp(lin(np.log(omega),*params1)), label="Fit-Funktion")
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r"$V=\frac{U_e}{U_a}$")
    plt.xlabel(r"$\omega  $ / $ Hz$")
    #plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
    #plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
    plt.tight_layout()
    plt.legend()
    plt.savefig("build/plots/" + name + ".pdf")
    
    return

intdiff(int_R1, int_C, int_f, int_Ue, int_Ua,"integrator")
intdiff(diff_R1, diff_C, diff_f, diff_Ue, diff_Ua,"differenzierer")















#### nicht-invertierende-Schmitt-Trigger ####
print("\n\n\n\t#### nicht-invertierende-Schmitt-Trigger ####\n")

def schmitt(x, y, R1, R2, Ue):
    
    kipp_theo = R1/R2 * Ue

    x_pos = np.array([x[0], x[-1]])
    y_pos = np.array([y[0], y[-1]])
    kipp_pos = ufloat(np.mean(y_pos), np.std(y_pos))

    x_neg = np.array([x[1], x[2]])
    y_neg = np.array([y[1], y[2]])
    kipp_neg = ufloat(np.mean(y_neg), np.std(y_neg))

    print(f"\nKipp_theo = {kipp_theo} V")
    print(f"\nkipp_pos = {kipp_pos} V")
    rel_abw(kipp_theo, kipp_pos)
    print(f"\nkipp_neg = {kipp_neg} V")
    rel_abw(kipp_theo, np.absolute(kipp_neg))

    return



schmitt(schmitt_x, schmitt_y, schmitt_R1, schmitt_R2, schmitt_Ue)











#### Signalgenerator ####
print("\n\n\n\t#### Signalgenerator ####\n")

def signal(R1, R2, R3, C):

    U_mess = 5.3 / 2
    f_mess = 1.64e3

    U_max = 27.5 / 2
    f_theo = R2 / (4 * C * R1 * R3 )
    Amp_theo = U_max * R1/R2

    print(f"\nf_theo = {f_theo:.4f} 1/s")
    rel_abw(f_theo, f_mess)
    print(f"Amp_theo = {Amp_theo:.4f} 1/s")
    rel_abw(U_max, U_mess)

    #print(f"\nm = {noms(uparams1[0]):.4f} \pm {stds(uparams1[0]):.4f} s\nn = {noms(uparams1[1]):.4f} \pm {stds(uparams1[1]):.4f}")

    return





signal(signal_R1, signal_R2, signal_R3, signal_C)

########Grafiken########


plt.figure()
plt.rc('axes', labelsize=size_label)
plt.plot(diff_f, (diff_Ua/diff_Ue), "x", label="Messwerte")
plt.ylabel(r"$V=\frac{U_e}{U_a}$")
plt.xlabel(r"$f  $ / $ Hz$")
#plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
#plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/hehehe.pdf")




######Messdaten###########
print("\n\n\n\t#########MESSDATEN#######")

#printer4(lin1_f, lin1_Ue, lin1_Ua, lin1_phase, "lin1")
#printer4(lin2_f, lin2_Ue, lin2_Ua, lin2_phase, "lin2")
#printer4(lin3_f, lin3_Ue, lin3_Ua, lin3_phase, "lin3")

#printer4(int_f, int_Ue, int_Ua, int_phase, "int")
#printer4(diff_f, diff_Ue, diff_Ua, diff_phase, "diff")

#printer2(schmitt_x, schmitt_y, "Schmitt")