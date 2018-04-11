import matplotlib.pyplot as plt
import snr_module_fft as snr
import numpy as np

#time, h = snr.HPlot(29,32,420,0.09,0,0,0,0)
#time, f = snr.FPlot(29, 32)
# plt.plot(-time, h)
# plt.plot(-time, f)

#freq, hfft = snr.HfftPlot(29,32,420,0.09,0,0,0,0)
freq1, asd = snr.ASD(1000)
# sqrt(f) factor needed because of the loglog plot 
# plt.loglog(freq,2*np.abs(hfft)*np.sqrt(freq))
plt.loglog(freq1, asd)

#snr.SNR(29,32,420,0.09,0,0,0,0)
freq, hfft, fs = snr.Hplusfft(29,32,420,0.09,0)

plt.loglog(freq, 2*np.abs(hfft)*np.sqrt(freq))
plt.axis([10, fs/2.0, 1e-24, 1e-18])
plt.xlabel('Freq (Hz)')
plt.ylabel('Strain / Hz')
plt.show()
