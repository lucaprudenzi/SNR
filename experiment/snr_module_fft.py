import math as m
import numpy as np
from tabulate import tabulate

# Constants
c = 2.99e10
freqin = 30 # Ligo inband frequency
G = 6.67e-8
Fs = 4096 # number of sample

## Masses functions
def SolarMass(mass):
    M_sun = 1.989e33 # Sun mass [g]
    return  mass*M_sun

def ChirpMass(mass1, mass2):
    mass1 = SolarMass(mass1)
    mass2 = SolarMass(mass2)
    return (mass1*mass2)**(3./5.)/(mass1+mass2)**(1./5.)

# symmetric mass ratio
def SMR(mass1, mass2):
    mass1 = SolarMass(mass1)
    mass2 = SolarMass(mass2)
    return (mass1*mass2)/(mass1+mass2)**2

## Distances functions
# distance in Mpc
def Mpc(distance):
    return distance*3.1e24

# proper distance
def dp(dl, z):
    dl = Mpc(dl)
    return dl/(1.+z)

# tau See Note 3 Chap 4 Maggiore
def Tau(Mc):
    tcoal = (1./np.pi)**(8/3)*(5./256.)*(G*Mc/c**3)**(-5./3)* \
            (1/freqin)**(8./3)
    N = Fs*tcoal 
    t = np.linspace(0, tcoal-0.01, N)
    tau = tcoal-t
    return tau

def PHI(Mc, t):
    Phi0 = 0 # Phase at coalescence 
    return -2*(5*G*Mc/c**3)**(-5./8.)*t**(5./8.)+Phi0

def Freq(Mc, t):
    return 1./m.pi*(5./256.*1./t)**(3./8.)*(G*Mc/c**3)**(-5./8.)

## Detector pattern functions
# F+
def Fplus(theta, phi, psi):
    return 0.5*(1+m.cos(theta)**2)*m.cos(2*phi)* \
            m.cos(2*psi)-m.cos(theta)*m.sin(2*phi)*m.sin(2*psi)

# Fx
def Fcross(theta, phi, psi):
    return 0.5*(1+m.cos(theta)**2)*m.cos(2*phi)* \
            m.sin(2*psi)+m.cos(theta)*m.sin(2*phi)*m.cos(2*psi)

## GW Amplitude
def Hplus(mass1, mass2, dl, z, iota):
    r = dp(dl, z) # proper distance
    Mc = ChirpMass(mass1, mass2)
    tau = Tau(Mc)
    Phi = PHI(Mc, tau)
    freq = Freq(Mc, tau)

    hplus = 4./r*(G*Mc/c**2)**(5./3.)*(m.pi*freq/c)**(2./3.)* \
            (1+np.cos(iota)**2)/2.*np.cos(Phi)
    return hplus

def Hcross(mass1, mass2, dl, z, iota):
    r = dp(dl, z) # proper distance
    Mc = ChirpMass(mass1, mass2)
    tau = Tau(Mc)
    Phi = PHI(Mc, tau)
    freq = Freq(Mc, tau)

    hcross = 4./r*(G*Mc/c**2)**(5./3.)*(m.pi*freq/c)**(2./3.)* \
            np.cos(iota)*np.sin(Phi)
    return hcross

def H(mass1, mass2, dl, z, iota, theta, phi, psi):
    hplus = Hplus(mass1, mass2, dl, z, iota)
    hcross = Hcross(mass1, mass2, dl, z, iota)
    fplus = Fplus(theta, phi, psi)
    fcross = Fcross(theta, phi, psi)

    return fplus*hplus+fcross*hcross

## FFT
def Hplusfft(mass1, mass2, dl, z, iota):
    hplus = Hplus(mass1, mass2, dl, z, iota)
    Mc = ChirpMass(mass1, mass2)
    tau = Tau(Mc)
    
    window = np.blackman(hplus.size)
    windowedhplus = hplus*window
    hplusfft = np.fft.rfft(windowedhplus)/Fs
    freq = np.fft.rfftfreq(len(windowedhplus))*Fs
    return freq, hplusfft, Fs

def Hcrossfft(mass1, mass2, dl, z, iota):
    hcross = Hcross(mass1, mass2, dl, z, iota)
    Mc = ChirpMass(mass1, mass2)
    tau = Tau(Mc)
    
    window = np.blackman(hcross.size)
    hcrosswindow = hplus*window
    hcrossfft = np.fft.rfft(hcrosswindow)/Fs
    freq = np.fft.rfftfreq(len(windowedhplus))*Fs
    return freq, hcrossfft, Fs

def Hfft(mass1, mass2, dl, z, iota, theta, phi, psi):
    freq, hplusfft, samplefreq = Hplusfft(mass1, mass2, dl, z, iota)
    freq, hcrossfft, samplefreq = Hcrossfft(mass1, mass2, dl, z, iota)
    fplus = Fplus(theta, phi, psi)
    fcross = Fcross(theta, phi, psi)
    
    return fplus*hplusfft+fcross*hcrossfft

# Amplitude spectral density
def ASD(freqfinal):
    # frequency interval
    freq = np.linspace(freqin,freqfinal,1e3)
    # power spectral density
    psd = (1.e-22*(18./(0.1+freq))**2)**2+0.7e-23**2+ \
            ((freq/2000.)*4.e-23)**2
    # ampliude spectral density
    asd = np.sqrt(psd)
    return freq, asd


def SNR(mass1, mass2, dl, z, iota, theta, phi, psi):
    Mc = ChirpMass(mass1, mass2)
    time = Tau(Mc)
    freq = Freq(Mc, time)
    freqfinal = freq[-1] #last element
    freq, asd = ASD(freqfinal)

    hfft = Hfft(mass1, mass2, dl, z, iota, theta, phi, psi)
    fraction = 4*(np.abs(hfft))**2/asd**2
    
    result = trapezoidal(fraction, freq)
    snr = np.sqrt(result)
    print (tabulate([['Mass1',mass1],['Mass2',mass2],['d_l',dl],['z',z],\
                     ['iota',iota],['theta',theta],['phi',phi],\
                     ['psi',psi],['SNR',snr]],\
                    headers=['Quantities','Values']))
    return snr

# Compute the integral
def trapezoidal(func, freq):
    n = len(func)
    dx = (freq[-1]-freq[0])/(n+1)
    s=0.0

    s=func[0]+func[-1]
    for i in range(1,n-1):
        s+=2*func[i]
    return s*dx/2.

## Plot functions
def HPlot(mass1, mass2, dl, z, iota, theta, phi, psi):
    h = H(mass1, mass2, dl, z, iota, theta, phi, psi)
    Mc = ChirpMass(mass1, mass2)
    time = Tau(Mc)

    return time, h

def FPlot(mass1, mass2):
    Mc = ChirpMass(mass1, mass2)
    time = Tau(Mc)
    freq = Freq(Mc, time)

    return time, freq

def HfftPlot(mass1, mass2, dl, z, iota, theta, phi, psi):
    hfft = Hfft(mass1, mass2, dl, z, iota, theta, phi, psi)
    Mc = ChirpMass(mass1, mass2)
    time = Tau(Mc)
    freq = Freq(Mc, time)

    return freq, hfft
