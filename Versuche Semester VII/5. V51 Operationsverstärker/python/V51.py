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


print(np.size(lin1_f), np.size(lin1_Ue), np.size(lin1_Ua))

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

print(np.size(lin2_f), np.size(lin2_Ue), np.size(lin2_Ua))


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

print(np.size(lin3_f), np.size(lin3_Ue), np.size(lin3_Ua))



#####RECHNUNGEN#######

########Grafiken########


#plt.figure()
#plt.plot(f,U/(10),"x",label="Messwerte")
#plt.plot(x,func(x,omega0) )
##plt.xscale('log')
#plt.xlabel(r"$\Omega=\frac{\nu}{\nu_0}$")
#plt.ylabel(r"$\frac{U_Br}{U_0}$")
##plt.xticks([5*10**3,10**4,2*10**4,4*10**4],[r"$5*10^3$", r"$10^4$", r"$2*10^4$", r"$4*10^4$"])
##plt.yticks([0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2],[r"$0$",r"$\frac{\pi}{8}$", r"$\frac{\pi}{4}$",r"$\frac{3\pi}{8}$", r"$\frac{\pi}{2}$"])
#plt.tight_layout()
#plt.legend()
#plt.savefig("build/plots/plot1.pdf")