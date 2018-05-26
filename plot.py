import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import snr_module as snr
import numpy as np
import string

# Parameters
geoparamfile = open("geoparams.txt", "r")
geoparams = [list(map(float, line.strip().split(' '))) for line in geoparamfile]
mass1= geoparams[0][0] # in Solar Masses
mass2 = geoparams[1][0] # in Solar Masses
dl = geoparams[2][0]
z = geoparams[3][0]
iota = geoparams[7][0] # angle between line of sight and angular momentum of the binary
iota = iota*np.pi/180
detector = geoparams[8][0]

detectorname = ''

if (detector == 1):
    detectorname = 'Hanford'
if (detector == 2):
    detectorname = 'Livingstone'
if (detector == 3):
    detectorname = 'Virgo'

geoparamfile.close()

# change last value of geoparams.txt in 1, 2 or 3 to choose H1, L1, V1
detanglesfile = open("detangles.txt")
detangles = [list(map(float, line.strip().split(' '))) for line in detanglesfile]

# Angles seen by detector 
theta = detangles[0][0] # angle from the z axis of interferometer to source position
phi = detangles[1][0] # angle from x-arm of the interferometer to source position
psi = detangles[2][0] # angle from the projected x-arm of the interferometer on binary plane
                      # and semijar axis of binary

theta = theta*np.pi/180.
phi = phi*np.pi/180.
psi = psi*np.pi/180.

detanglesfile.close()

time, h = snr.HPlot(mass1,mass2,dl,z,iota,theta,phi,psi)
time, freq = snr.FPlot(mass1, mass2, z)

freq1, hft  = snr.HftPlot(mass1,mass2,dl,z,iota, theta, phi, psi)
freq2, asd = snr.ASDPlot()

# Print a table with data and SNR
# print table with results (needs tabulate package) 
snr_result = snr.SNR(mass1,mass2,dl,z,iota, theta, phi, psi, detectorname, 1)

# print only snr (no additional package needed)
# snr_result = snr.SNR(mass1,mass2,dl,z,iota,theta,phi,psi)

## Plots
fig = plt.figure(figsize=(12,5))
gs = gridspec.GridSpec(2, 4)
gs.update(wspace=.7)
fig.suptitle('mass1: %.1f, mass2: %.1f, dl: %.1fMpc, z: %.1f, iota: %.1f, \
             theta: %.1f, phi: %.1f, psi: %.1f'
             %(mass1,mass2,dl,z,iota*180/np.pi,theta*180/np.pi,phi*180/np.pi,psi*180/np.pi), fontsize=11)
# remove vertical gap between subplot
plt.subplots_adjust(hspace=.0)
# First left hand plot
ax0 = plt.subplot(gs[0,:2])
line0, = ax0.plot(-time, h, color='r', label='h')
plt.ylabel('Strain h')
plt.grid(True)
# Second left hand plot
ax1 = plt.subplot(gs[1,:2], sharex=ax0)
line1, = ax1.plot(-time, freq, color='g', label='freq')
plt.xlabel('Time [s]')
plt.ylabel('Frequency [Hz]')

ax1.legend((line0, line1), ('h', 'Frequency [Hz]'), loc='upper left')
plt.grid(True)
plt.setp(ax0.get_xticklabels(), visible=False)
# remove last tick on below plot
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)

# Right hand plot
ax2 = plt.subplot(gs[0:,2:])
# sqrt(f) factor needed because of the loglog plot
plt.loglog(freq1, 2*np.abs(hft)*np.sqrt(freq1), \
           label=r'$2\|\mathrm{hft}\|\sqrt{\mathrm{freq}}$', color='c')
plt.loglog(freq2, asd, label='ASD', color='m')
plt.xlabel('Frequency [Hz]')
ax2.set_ylabel(r'$2\|hft\|\sqrt{\mathrm{freq}}\mathrm{,}\sqrt{S_{noise}}[\mathrm{strain}/\sqrt{\mathrm{Hz}}]$')
plt.legend()

plt.show()
