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


#Laden der Daten aus der Datei "echo_gradient.csv" 
#Die erste Spalte enthält die Zeiten in Sekunden, die zweite Spalte  
#den Realteil und die dritte Spalte den Imaginärteil 
times, real, imag= np.genfromtxt("python/data/echo_gradient.csv", skip_header= 5, unpack = True)



#######DIFFUSIONSKONSTANTE ÜBER BESSEL#########
print("\n\t##### Diffusionskonstante über Bessel  ######")

peaks_real, peak_heights_real = find_peaks( np.real(real), height = 1.386)

print(f"peak height: {peak_heights_real} {peaks_real}")
max_height =np.real(real[peaks_real[0]])
print(f"max height / 2: {max_height/2 }")
limit = 0.002513


mask1 = np.real(real) >= max_height/2 - limit 
mask2 = np.real(real) <= max_height/2 + limit

half = np.real(real[mask1 * mask2])
print(f"half: {half} ")
#print(max_height/2 - half)

FWHM = (times[mask1 * mask2])[1] - (times[mask1 * mask2])[0] #in ms
uFWHM = ufloat(FWHM, times[1] - times[0])

print(f"Halbwertsbreite: {FWHM} \pm ({times[1] - times[0]}) ms")

d = 4.2*10**(-3) # in m
gyro = 2.67*10**8 #für Protonen T/s

G_fwhm = 4 * 2.2 / ( d * gyro * uFWHM)

print(f"Gradient über Halbwertsbreite: {noms(G_fwhm)} \pm {stds(G_fwhm) }")

plt.figure()
plt.plot(times, real, label = "Realteil") 
plt.plot(times, imag, label = "Imaginärteil" )
plt.plot(times[mask1 * mask2], real[mask1 * mask2], label = "Halbwertsbreite")
plt.xlabel(r"$\tau$ $/$ $ms$") 
plt.ylabel(r"$U$ $/$ $V$")
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/signal.pdf")


#Suchen des Echo-Maximums und alle Daten davor abschneiden 
start = np.argmax(real) 
times = times[start:] 
real = real[start:] 
imag = imag[start:] 
#Phasenkorrektur - der Imaginärteil bei t=0 muss = 0 sein 
phase = np.arctan2(imag[0], real[0]) 
#Daten in komplexes Array mit Phasenkorrektur speichern 
compsignal = (real*np.cos(phase)+imag*np.sin(phase))+ (-real*np.sin(phase)+imag*np.cos(phase))*1j 
#Offsetkorrektur, ziehe den Mittelwert der letzten 512 Punkte von allen Punkten ab 
compsignal = compsignal - compsignal[-512:-1].mean() 
#Der erste Punkt einer FFT muss halbiert werden 
compsignal[0] = compsignal[0]/2.0 
#Anwenden einer Fensterfunktion (siehe z. Bsp. #https://de.wikipedia.org/wiki/Fensterfunktion ) 
#Hier wird eine Gaußfunktion mit sigma = 100 Hz verwendet 
apodisation = 100.0*2*np.pi 
compsignal = compsignal*np.exp(-1.0/2.0*((times-times[0])*apodisation)**2) 
#Durchführen der Fourier-Transformation 
fftdata = np.fft.fftshift(np.fft.fft(compsignal)) 
#Generieren der Frequenzachse 
freqs = np.fft.fftshift(np.fft.fftfreq(len(compsignal), times[1]-times[0])) 
#Speichern des Ergebnisses als txt 
np.savetxt("build/plots/echo_gradient_fft.txt", np.array([freqs, np.real(fftdata), np.imag(fftdata)]).transpose()) 
#Erstellen eines Plots 








#######MINE########



peaks, peak_heights = find_peaks( -1 * np.real(fftdata), height = 18)
print(np.real(fftdata)[peaks])

#print(peaks)
save_peaks = peaks
save_peaks[0] += 2
#save_peaks[1] -=1
#print(save_peaks)
#print(peaks)


lower = 5050
upper = 5200
plt.figure()
plt.plot(freqs[lower:upper], np.real(fftdata)[lower:upper], label = "Fouriertransformierte Daten") 
plt.plot(freqs[peaks], np.real(fftdata)[peaks], label = "Peakbreite" )
plt.xlabel(r"$f$ $/$ $Hz$")
plt.ylabel(r"")
plt.tight_layout()
plt.legend()
plt.savefig("build/plots/echo_gradient.pdf")

#plt.axis([0,17500,-7,30])

####G berechnen

d = 4.2*10**(-3) # in m
gyro = 2.67*10**8 #für Protonen T/s
d_f = 0 #Durchmesser Spektrum 



d_f = freqs[save_peaks[1]] - freqs[save_peaks[0]] #Durchmesser Spektrum 
ud_f = ufloat(d_f, freqs[1] - freqs[0])

G = (2 * np.pi * ud_f)/(gyro * d)
print(f"\nVerteilungsdicke d_f = {d_f:.4f} \pm {stds(ud_f)} \nGradientenstärke G = {noms(G)} \pm {stds(G)}\n")

