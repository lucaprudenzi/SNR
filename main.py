import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import snr_module as snr
import numpy as np

# Parameters
mass1 = 7 # Solar masses
mass2 = 14 # Solar masses
dl = 440 # Mpc
z = 0.09 
iota = 0 # angle between line of sight and angular momentum of the binary
theta = 0 # angle from the z axis, perpendicular the interferometer
phi = 0 # angle from x-arm of the interferometer
psi = 0 # orientation of the axes of the binary respect to the axes of the
        # interferometers

time, h = snr.HPlot(mass1,mass2,dl,z,iota,theta,phi,psi)
time, freq1 = snr.FPlot(mass1, mass2)

freq2, hft  = snr.HftPlot(mass1,mass2,dl,z,iota,theta,phi,psi)
freq3, asd = snr.ASDPlot()

# Print a table with data and SNR
snr.SNR(mass1,mass2,dl,z,iota,theta,phi,psi)
 
## Plots
fig = plt.figure(figsize=(12,5))
gs = gridspec.GridSpec(2, 4)
gs.update(wspace=.7)
fig.suptitle('mass1: %.1f, mass2: %.1f, dl: %.1fMpc'%(mass1,mass2,dl), fontsize=10)
# remove vertical gap between subplot
plt.subplots_adjust(hspace=.0)
# Left hand plot
ax0 = plt.subplot(gs[0,:2])
line0, = ax0.plot(-time, h, color='r', label='h')
plt.ylabel('h')
plt.grid(True)

ax1 = plt.subplot(gs[1,:2], sharex=ax0)
line1, = ax1.plot(-time, freq1, color='g', label='freq')
plt.xlabel('Time')
plt.ylabel('Frequency')

ax1.legend((line0, line1), ('h', 'Frequency'), loc='upper left')

plt.grid(True)
plt.setp(ax0.get_xticklabels(), visible=False)
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)

# Right hand plot
plt.subplot(gs[0:,2:])
#sqrt(f) factor needed because of the loglog plot
plt.loglog(freq2, 2*np.abs(hft)*np.sqrt(freq2), \
           label='2abs(hfft)sqrt(freq)', color='c')
plt.loglog(freq3, asd, label='sqrt(S_noise)', color='m')
plt.xlabel('Frequency')
plt.legend()

plt.show()
