import numpy as np
import matplotlib.pyplot as plt
import math

## CONSTANTS
A = np.pi**(-2./3)*(5./24)**(1./2)
c = 2.99792458e10 # [cm s-1]
G  = 6.67259e-8 # [erg cm g-2]
a = 1/(1+0) # 1/(1+z0) scale factor today
M_sun = 1.989e33 # Sun mass [g]
pc = 3.1e18
f_i = 30 # Ligo in-band frequency

## INTEGRATION PARAMETERS
tau_0 = 0

## PARAMETERS
# DISTANCES
z = 0.09 # redshift
r = 420e6*pc  # distance of binary from us [cm]
d_l = (1+z)*a*r  # luminosity distance (distance with redsift) [cm]

# MASSES
m_p = 36.2*M_sun # primary mass [g]
m_s = 29.1*M_sun # secondary mass [g]
nu = (m_p*m_s)/(m_p+m_s)**2 # symmetric mass ratio [dimensionless] [0.01 - 0.25]
Mc = (m_p*m_s)**(3./5.)/(m_p+m_s)**(1./5.) # chirp mass [g]

# ANGLES
# iota = nL : angle between the line of sight n and the unit vector parallel to the orbital angular momentum L. iota = 0 is a head-on observation
iota = 0
phi_0 = 0 # used in phi_plus

# ANTENNA PARAMETERS
theta =  0 # angle from z axis
phi = 0 # angle on xy plane
# F_for and F_plus are referred to detectores own axis (psi = 0)
F_cross = 1./2.*(1+(np.cos(theta))**2)*np.cos(2.*phi) # maximum value is 1
F_plus = np.cos(theta)*np.sin(2.*phi) # maximum value is 1

## FREQUENCIES DOMAIN
# bisogna usare la f fino alla frequenza di coalescenza, mentre  dopo quella f il contributo del segnale e' nullo
f_signal = c*c*c/(6*math.sqrt(6)*2*math.pi*G*(m_p+m_s))
vector_f = np.linspace(10, f_signal, 1e3)
print(f_signal)
f = np.logspace(1,3,1e3)
# time
t_coal = 5./256/(np.pi**(8./3))*(c**3/(G*Mc))**(5./3)*f_i**(-8./3)
# Terms for PN computation
tau0 = 5./256./np.pi/f_i*(np.pi*G*(m_p+m_s)/c**3*f_i)**(-5./3)/nu
tau1 = 5./192./np.pi/f_i*(np.pi*G*(m_p+m_s)/c**3*f_i)**(-1)/nu*(743./336.+11./4.*nu)
tau15 = 1./8./f_i*(np.pi*G*(m_p+m_s)/c**3*f_i)**(-2./3.)/nu
tau2 = 5./128./np.pi/f_i*(np.pi*G*(m_p+m_s)/c**3*f_i)**(-1./3)/nu*(3058673./1016064.+5429./1008.*nu+617./144.*nu**2)

## H_TILDE
PSI_plus = 2*np.pi*f*(t_coal+r/c)-phi_0-np.pi/4+2*np.pi*f_i\
*(3./5*tau0*(f/f_i)**(-5./3)+tau1*(f/f_i)**(-1)\
-3./2*tau15*(f/f_i)**(-2./3)+3*tau2*(f/f_i)**(-1./3))
PSI_cross = PSI_plus + np.pi/2

h_tilde_plus = (5./6)**(1./2)*c/(2*np.pi**(2./3)*r)*(G*Mc/c**3)**(5./6)*f**(-7./6)*np.exp(1j*PSI_plus)*(1+np.cos(iota)**2)/2

h_tilde_cross = (5./6)**(1./2)*c/(2*np.pi**(2./3)*d_l)*(G*Mc/c**3)**(5./6)*f**(-7./6)*np.exp(1j*PSI_plus)*np.cos(iota)

h_tilde = F_cross*h_tilde_cross+F_plus*h_tilde_plus

## ASD
# usa quella del tutorial di Ligo
f0_aligo = 215
S0_aligo = 10**(-49)
F_aligo = f/f0_aligo # dimensionless frequency
Snoise = (F_aligo**(-4.14)-5*F_aligo**(-2)+111*(1-F_aligo**2+0.5*F_aligo**4)/(1+0.5*F_aligo**2))*S0_aligo

# SNR
# vedi 7.51: bisogna fare 2xla radice dell'integrale da 0 fino alla frequenza di merger del rapporto tra il modulo quadro di h tilde e S

## PLOT
plt.loglog(vector_f, np.absolute(h_tilde)*f**(1./2),label=r'$m_1$='+str(int(m_p/M_sun))+'$M_{\odot}, m_2$='+str(int(m_s/M_sun))+'$M_{\odot}$')
# plt.axis([99.,100,10**(-22),10**(-21)])
plt.legend()
plt.xlabel('$f$ [Hz]')
plt.ylabel(r'$|h|\sqrt{f}$ and $\sqrt{S_{noise}}$')
plt.loglog(f,np.sqrt(Snoise))

plt.savefig('asd.pdf',format='pdf')
plt.show()
