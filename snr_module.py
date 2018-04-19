import math as m
import numpy as np

# Constants
c = 3e10
M_sun = 1.989e33 # Sun mass [g]
freqin = 1 # Ligo inband frequency
G = 6.67e-8
N = 150000 # points in time array
# choose between H1, L1 and V1

# detector = 'L1'
detector = 'H1'
# detector = 'V1'

def SolarMass(mass):
    return  mass*M_sun

def ChirpMass(mass1, mass2, z):
    mass1 = SolarMass(mass1)
    mass2 = SolarMass(mass2)
    return (mass1*mass2)**(3./5.)/(mass1+mass2)**(1./5.)*(1+z)

## Distances functions
# distance in Mpc
def Mpc(distance):
    return distance*3.1e24

# tau See Note 3 Chap 4 Maggiore
def Tau(Mc, mass1, mass2, z):
    mass1 = SolarMass(mass1)*(1+z)
    mass2 = SolarMass(mass2)*(1+z)
    # at 10Rs 
    #tbefore = 10*5./16*G*(mass1+mass2)**2/(c**3*(mass1*mass2)/(mass1+mass2))
    #print(tbefore)
    fisco = 2*1./(6*m.sqrt(6)*2*m.pi)*(c**3/(G*(mass1+mass2)))
    #print(fisco)
    #print(1./fisco)

    tcoal = (1./np.pi)**(8/3)*(5./256.)*(G*Mc/c**3)**(-5./3)* \
            (1./freqin)**(8./3)
    tcoalfin = (1./np.pi)**(8/3)*(5./256.)*(G*Mc/c**3)**(-5./3)* \
            (1./fisco)**(8./3)

    #print(tcoalfin)
    t = np.linspace(0, tcoal-tcoalfin, N)
    tau = tcoal-t
    return tau

def Freq(Mc, t, mass1 = SolarMass(15), mass2 = SolarMass(7)):
       return 1./m.pi*(5./256.*1./t)**(3./8.)*(G*Mc/c**3)**(-5./8.)

## Detector pattern functions (Schutz11)
# F+
def Fplus(theta, phi, psi, detector='H1'):
    
    if detector=='H1':
        latitude, longitude, rotation = H1position()
    elif detector=='L1':
        latitude, longitude, rotation = L1position()
    elif detector=='V1':
        latitude, longitude, rotation = V1position()
    
    # theta has 0 on the z-axis, latitude on the x-y plane
    theta = theta - (np.pi/2-latitude)
    # phi - longitudine = x-axis points toward south
    # rotation = x-axis rotation from South
    phi = phi - longitude - rotation
    # psi trasforms in the same way as phi, because projection on
    # the sky plane doesn't change angles
    psi = psi - longitude - rotation

    return 0.5*(1+m.cos(theta)**2)*m.cos(2*phi)* \
            m.cos(2*psi)-m.cos(theta)*m.sin(2*phi)*m.sin(2*psi)

# Fx
def Fcross(theta, phi, psi, detector='H1'):
 
    if detector=='H1':
        latitude, longitude, rotation = H1position()
    elif detector=='L1':
        latitude, longitude, rotation = L1position()
    elif detector=='V1':
        latitude, longitude, rotation = V1position()

    # theta has 0 on the z-axis, latitude on the x-y plane
    theta = theta - (np.pi/2-latitude)
    # phi - longitudine = x-axis points toward south
    # rotation = x-axis rotation from South
    phi = phi - longitude - rotation
    # psi trasforms in the same way as phi, because projection on
    # the sky plane doesn't change angles 
    psi = psi - longitude - rotation

    return 0.5*(1+m.cos(theta)**2)*m.cos(2*phi)* \
            m.sin(2*psi)+m.cos(theta)*m.sin(2*phi)*m.cos(2*psi)

## GW Amplitude (Maggiore)
def PHI(Mc, t):
    Phi0 = 0 # Phase at coalescence 
    return -2*(5*G*Mc/c**3)**(-5./8.)*t**(5./8.)+Phi0

def Hplus(mass1, mass2, dl, z, iota):
    Mc = ChirpMass(mass1, mass2, z)
    tau = Tau(Mc, mass1, mass2, z)
    Phi = PHI(Mc, tau)
    freq = Freq(Mc, tau)
    dl = Mpc(dl)

    hplus = 4./dl*(G*Mc/c**2)**(5./3.)*(m.pi*freq/c)**(2./3.)* \
            (1+np.cos(iota)**2)/2.*np.cos(Phi)
    
    return hplus

def Hcross(mass1, mass2, dl, z, iota):
    Mc = ChirpMass(mass1, mass2, z)
    tau = Tau(Mc, mass1, mass2, z)
    Phi = PHI(Mc, tau)
    freq = Freq(Mc, tau)
    dl = Mpc(dl)

    hcross = 4./dl*(G*Mc/c**2)**(5./3.)*(m.pi*freq/c)**(2./3.)* \
            np.cos(iota)*np.sin(Phi)
    
    return hcross

def H(mass1, mass2, dl, z, iota, theta, phi, psi):
    hplus = Hplus(mass1, mass2, dl, z, iota)
    hcross = Hcross(mass1, mass2, dl, z, iota)
    fplus = Fplus(theta, phi, psi, detector)
    fcross = Fcross(theta, phi, psi, detector)
    
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
    Mc = ChirpMass(mass1, mass2, z)
    tau = Tau(Mc, mass1, mass2, z)
    freq = Freq(Mc, tau)
    dl = Mpc(dl)
    #
    freq = np.logspace(np.log10(freqin), np.log10(freq[-1]), N)
    #
    Psiplus = PSIplus(Mc, freq)
    A = 1./m.pi**(2./3.)*(5./24.)**(1./2.)
    hplusft = A*np.exp(1j*Psiplus)*c/dl*(G*Mc/c**3)**(5./6.)*\
            freq**(-7./6.)*(1+np.cos(iota)**2)/2.
    
    return hplusft

def Hcrossft(mass1, mass2, dl, z, iota):
    Mc = ChirpMass(mass1, mass2, z)
    tau = Tau(Mc, mass1, mass2, z)
    freq = Freq(Mc, tau)
    dl = Mpc(dl)
    #
    freq = np.logspace(np.log10(freqin), np.log10(freq[-1]),N)
    #
    Psicross = PSIcross(Mc, freq)
    A = 1./m.pi**(2./3)*(5./24.)**(1./2.)
    hcrossft = A*np.exp(1j*Psicross)*c/dl*(G*Mc/c**3)**(5./6.)*\
            freq**(-7./6.)*np.cos(iota)

    return hcrossft

def Hft(mass1, mass2, dl, z, iota, theta, phi, psi):
    hplusft = Hplusft(mass1, mass2, dl, z, iota)
    hcrossft = Hcrossft(mass1, mass2, dl, z, iota)
    fplus = Fplus(theta, phi, psi, detector)
    fcross = Fcross(theta, phi, psi, detector)

    return fplus*hplusft+fcross*hcrossft

def PSDInput(freq):
    # power spectral density (from Ligo tutorial)
    psd = (1.e-22*(18./(0.1+freq))**2)**2+0.7e-23**2+ \
            ((freq/2000.)*4.e-23)**2
    return psd

# Amplitude spectral density (hard-coded model from Ligo tutorial)
def ASD(mass1, mass2, z):
    Mc = ChirpMass(mass1, mass2, z)
    tau = Tau(Mc, mass1, mass2, z)
    # same interval of the Fourier transform of h 
    freq = Freq(Mc, tau)
    #
    freq = np.logspace(np.log10(freqin), np.log10(freq[-1]),N)
    #
    psd = PSDInput(freq)
    # ampliude spectral density
    asd = np.sqrt(psd)
    
    return freq, asd

def SNR(mass1, mass2, dl, z, iota, theta, phi, psi, table=0):
    freq, asd = ASD(mass1, mass2, z)
    hft = Hft(mass1, mass2, dl, z, iota, theta, phi, psi)

    # check that the signal is above the sensibility curve
    factor1 = 2*(np.abs(hft))*np.sqrt(freq)
    factor2 = asd
    mask = factor1 > factor2
   
    #
    factor1 = 2*(np.abs(hft)) # this is the function to evaluate 
    #                         # in the integral
    # reduce the evaluation interval only to inband part
    factor1 = factor1[mask]
    factor2 = factor2[mask]
    freq = freq[mask]
    fraction = factor1**2/factor2**2 
    
    # Integral computation (two equivalent methods)
    integral = np.trapz(freq*fraction, np.log(freq))
    #integral = np.trapz(fraction, freq)
    
    snr = np.sqrt(integral)
    if table == 1:
        from tabulate import tabulate
        print (tabulate([['Mass1 (Solar masses)',mass1],\
                     ['Mass2 (Solar masses)',mass2],\
                     ['d_l (Mpc)',dl],['z',z],\
                     ['iota (between n and L)',round(iota*180./np.pi,2)],\
                     ['theta (from z-axis)', round(theta*180./np.pi,2)],\
                     ['phi (from x-arm)',round(phi*180./np.pi,2)],\
                     ['psi (binary axes orientation)',\
                      round(psi*180/np.pi,2)],\
                     ['frequency in band', freqin],\
                     ['Detector', detector],\
                     ['SNR (inspiral)',round(snr,2)]], \
                    headers=['Quantities','Values']))
    else:
        print('SNR = ', snr)

    return snr

## Plot functions
def HPlot(mass1, mass2, dl, z, iota, theta, phi, psi):
    h = H(mass1, mass2, dl, z, iota, theta, phi, psi)
    Mc = ChirpMass(mass1, mass2, z)
    time = Tau(Mc, mass1, mass2, z)
    
    return time, h

def FPlot(mass1, mass2, z):
    Mc = ChirpMass(mass1, mass2, z)
    time = Tau(Mc, mass1, mass2, z)
    freq = Freq(Mc, time)

    return time, freq

def HftPlot(mass1, mass2, dl, z, iota, theta, phi, psi):
    hft = Hft(mass1, mass2, dl, z, iota, theta, phi, psi)
    Mc = ChirpMass(mass1, mass2, z)
    time = Tau(Mc, mass1, mass2, z)
    freq = Freq(Mc, time)
    #
    freq = np.logspace(np.log10(freqin), np.log10(freq[-1]),N)
    #

    return freq, hft

def ASDPlot():
    freqfinal = 2000
    # frequency interval
    freq = np.linspace(freqin,freqfinal,N)
    # power spectral density
    psd = (1.e-22*(18./(0.1+freq))**2)**2+0.7e-23**2+ \
            ((freq/2000.)*4.e-23)**2
    # ampliude spectral density
    asd = np.sqrt(psd)
    return freq, asd

# Detector positions

def H1position():
    latitude = (49+27./60.+19./3600.)*m.pi/180
    longitude = 2*m.pi - (119+44./60.+28./3600)*m.pi/180.
    # x-arm respect to South position
    rotation = (180+36)*m.pi/180.

    return latitude, longitude, rotation

def L1position():
    latitude = (30+33./60.+46./3600.)*m.pi/180
    longitude = 2*m.pi - (90+46./60.+27./3600)*m.pi/180
    # x-arm orientation respect to South position
    rotation = (270+18)*m.pi/180.

    return latitude, longitude, rotation

def V1position():
    latitude = (43+37./60.+53./3600.)*m.pi/180.
    longitude = (10+30./60.+16./3600)*m.pi/180.
    # x-arm orientation respect to South position
    rotation = (180-19)*m.pi/180.

    return latitude, longitude, rotation
