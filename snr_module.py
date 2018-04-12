import math as m
import numpy as np
from tabulate import tabulate

# Constants
c = 2.99e10
freqin = 30 # Ligo inband frequency
G = 6.67e-8
N = 4096
tbefore = 0.01 # time before divergence of linearized quantities

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
    t = np.linspace(0, tcoal-tbefore, N)
    tau = tcoal-t
    return tau

def Freq(Mc, t):
    return 1./m.pi*(5./256.*1./t)**(3./8.)*(G*Mc/c**3)**(-5./8.)

## Detector pattern functions (Schutz11)
# F+
def Fplus(theta, phi, psi):
    return 0.5*(1+m.cos(theta)**2)*m.cos(2*phi)* \
            m.cos(2*psi)-m.cos(theta)*m.sin(2*phi)*m.sin(2*psi)

# Fx
def Fcross(theta, phi, psi):
    return 0.5*(1+m.cos(theta)**2)*m.cos(2*phi)* \
            m.sin(2*psi)+m.cos(theta)*m.sin(2*phi)*m.cos(2*psi)

## GW Amplitude
def PHI(Mc, t):
    Phi0 = 0 # Phase at coalescence 
    return -2*(5*G*Mc/c**3)**(-5./8.)*t**(5./8.)+Phi0

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

## Fourier transforms (Maggiore)
def PSIplus(Mc, freq):
    Phi0 = 0
    return  2*np.pi*freq-Phi0-np.pi/4.+3./4.*\
            (G*Mc/c**3*8*np.pi*freq)**(-5./3.)

def PSIcross(Mc, freq):
    Psiplus = PSIplus(Mc, freq)
    
    return  Psiplus+np.pi/2.

def Hplusft(mass1, mass2, dl, z, iota):
    Mc = ChirpMass(mass1, mass2)
    r = dp(dl, z)
    tau = Tau(Mc)
    freq = Freq(Mc, tau)
    #
    freq = np.linspace(freqin, freq[-1], N)
    #
    Psiplus = PSIplus(Mc, freq)
    A = 1./m.pi**(2./3.)*(5./24.)**(1./2.)
    hplusft = A*np.exp(1j*Psiplus)*c/r*(G*Mc/c**3)*\
            freq**(-7./6.)*(1+np.cos(iota)**2)/2.
    
    return hplusft

def Hcrossft(mass1, mass2, dl, z, iota):
    Mc = ChirpMass(mass1, mass2)
    r = dp(dl, z)
    tau = Tau(Mc)
    freq = Freq(Mc, tau)
    #
    freq = np.linspace(freqin, freq[-1],N)
    #
    Psicross = PSIcross(Mc, freq)
    A = 1./m.pi**(2./3)*(5./24.)**(1./2.)
    hcrossft = A*np.exp(1j*Psicross)*c/r*(G*Mc/c**3)*\
            freq**(-7./6.)*np.cos(iota)

    return hcrossft

def Hft(mass1, mass2, dl, z, iota, theta, phi, psi):
    hplusft = Hplusft(mass1, mass2, dl, z, iota)
    hcrossft = Hcrossft(mass1, mass2, dl, z, iota)
    fplus = Fplus(theta, phi, psi)
    fcross = Fcross(theta, phi, psi)

    return fplus*hplusft+fcross*hcrossft

# Amplitude spectral density (hard-coded model from Ligo tutorial)
def ASD(mass1, mass2):
    Mc = ChirpMass(mass1, mass2)
    tau = Tau(Mc)
    # same interval of the Fourier transform of h 
    freq = Freq(Mc, tau)
    #
    freq = np.linspace(freqin, freq[-1],N)
    #
    # power spectral density
    psd = (1.e-22*(18./(0.1+freq))**2)**2+0.7e-23**2+ \
            ((freq/2000.)*4.e-23)**2
    # ampliude spectral density
    asd = np.sqrt(psd)
    
    return freq, asd

def SNR(mass1, mass2, dl, z, iota, theta, phi, psi):
    freq, asd = ASD(mass1, mass2)
    hft = Hft(mass1, mass2, dl, z, iota, theta, phi, psi)

    # check that the signal is above the sensibility curve
    factor1 = 2*(np.abs(hft))*np.sqrt(freq)
    factor2 = asd
    mask = factor1 > factor2
    factor1 = factor1[mask]
    factor2 = factor2[mask]
    freq = freq[mask]

    # integral computation
    fraction = factor1**2/factor2**2
    integral = trapezoidal(fraction, freq)
    snr = np.sqrt(integral)

    print (tabulate([['Mass1 (Solar masses)',mass1],\
                     ['Mass2 (Solar masses)',mass2],\
                     ['d_l (Mpc)',dl],['z',z],\
                     ['iota (between n and L)',iota],\
                     ['theta (from z-axis)', theta],\
                     ['phi (from x-arm)',phi],\
                     ['psi (binary axes orientation)',psi],\
                     ['frequency in band', freqin],\
                     ['SNR (inspiral)',snr]], \
                    headers=['Quantities','Values']))

    return snr

# Compute the integral
def trapezoidal(func, freq):
    n = len(func)
    print(n)
    if n==0:
        return 0
    freqlog = np.log(freq)

    s=func[0]*(freqlog[1]-freqlog[0])/2.
    s+=func[-1]*(freqlog[n-1]*freqlog[n-2])/2.
    for i in range(1,n-1):
        s+=func[i]*(freqlog[i+1]-freqlog[i])
    
    return s

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

def HftPlot(mass1, mass2, dl, z, iota, theta, phi, psi):
    hft = Hft(mass1, mass2, dl, z, iota, theta, phi, psi)
    Mc = ChirpMass(mass1, mass2)
    time = Tau(Mc)
    freq = Freq(Mc, time)
    #
    freq = np.linspace(freqin, freq[-1],N)
    #

    return freq, hft

def ASDPlot():
    freqfinal = 1000
    # frequency interval
    freq = np.linspace(freqin,freqfinal,N)
    # power spectral density
    psd = (1.e-22*(18./(0.1+freq))**2)**2+0.7e-23**2+ \
            ((freq/2000.)*4.e-23)**2
    # ampliude spectral density
    asd = np.sqrt(psd)
    return freq, asd


