import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import snr_module as snr
import numpy as np

# Parameters
mass1 =  # Solar masses
mass2 =  # Solar masses
dl =  # Mpc
z = 
iota =  # angle between line of sight and angular momentum of the binary
theta = # angle from the z axis, perpendicular the interferometer
phi =  # angle from x-arm of the interferometer
psi =  # orientation of the axes of the binary respect to the axes of the
        # interferometers

time, h = snr.HPlot(mass1,mass2,dl,z,iota,theta,phi,psi)
time, freq = snr.FPlot(mass1, mass2)

freq1, hft  = snr.HftPlot(mass1,mass2,dl,z,iota,theta,phi,psi)
freq2, asd = snr.ASDPlot()

# Print a table with data and SNR
snr.SNR(mass1,mass2,dl,z,iota,theta,phi,psi)
 
## Plots
fig = plt.figure(figsize=(12,5))
gs = gridspec.GridSpec(2, 4)
gs.update(wspace=.7)
fig.suptitle('mass1: %.1f, mass2: %.1f, dl: %.1fMpc, z: %.1f, iota: %.1f, theta: %.1f, phi: %.1f, psi: %.1f' \
             %(mass1,mass2,dl,z,iota,theta,phi,psi), fontsize=11)
# remove vertical gap between subplot
plt.subplots_adjust(hspace=.0)
# First left hand plot
ax0 = plt.subplot(gs[0,:2])
line0, = ax0.plot(-time, h, color='r', label='h')
plt.ylabel('h')
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
plt.subplot(gs[0:,2:])
# sqrt(f) factor needed because of the loglog plot
plt.loglog(freq1, 2*np.abs(hft)*np.sqrt(freq1), \
           label='2|hft|sqrt(f)', color='c')
plt.loglog(freq2, asd, label='sqrt(S_noise)', color='m')
plt.xlabel('Frequency [Hz]')
plt.ylabel('2|hft|sqrt(f) and sqrt(S_noise) [1/sqrt(Hz)]')
plt.legend()

plt.show()