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


#Laden der Daten aus der Datei "echo_gradient.csv" 
#Die erste Spalte enthält die Zeiten in Sekunden, die zweite Spalte  
#den Realteil und die dritte Spalte den Imaginärteil 
times, real, imag= np.genfromtxt("python/data/echo_gradient.csv", skip_header= 5, unpack = True)


plt.figure()
plt.plot(times, real, label = "Realteil") 
plt.plot(times, imag, label = "Imaginärteil" )
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

lower = 5050
upper = 5200


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


peaks, peak_heights = find_peaks( -1 * np.real(fftdata), height = 18)
print(np.real(fftdata)[peaks])

d_f = freqs[peaks[1]] - freqs[peaks[0]] #Durchmesser Spektrum 

G = (2 * np.pi * d_f)/(gyro * d)
print(f"\nVerteilungsdicke d_f = {d_f:.4f}\nGradientenstärke G = {G:.4f}")
