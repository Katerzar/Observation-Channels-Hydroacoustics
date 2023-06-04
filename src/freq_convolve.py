import params
from scipy.fft import fft,ifft
import numpy as np
from matplotlib import pyplot as plt

FFT = [ [0]*params.n for i in range(params.N)]
for i in range(params.N):
    y_signal = params.signals_array_FFT[i]
    FFT[i] = fft(y_signal, params.n)



n = 2**(int(np.log2(params.N))+1)



FFT_signals_fk = [0]*n
for t in range(params.N):
    FFT_signals_fk[t] = FFT[t][params.fk_index]


signals_fk_FFT = fft(FFT_signals_fk,n)

t_diff_L_coeff = params.t_diff_arr[params.m_index:params.m_index+params.M]
L_coeff = np.exp(1j * 2 * np.pi * params.frequency_FFT * np.asarray(t_diff_L_coeff))
Ln_coeff = [0]*n
for j in range(len(L_coeff)):
    Ln_coeff[j] = L_coeff[j]
Ln_coeff_FFT = fft(Ln_coeff,n)


Y = [signals_fk_FFT[i]*Ln_coeff_FFT[i] for i in range(n)]
#Y = np.convolve(FFT_signals_fk,L_coeff)

print("Y_mult")
print("len(Y_mult) = "+str(len(Y)))
print(Y)

Y_ifft = ifft(Y,n)
print("Y_ifft")
print("len(Y_ifft) = "+str(len(Y_ifft)))
x_ifft = np.linspace(0,n-1,n)
alpha_x = [round((params.alpha[i] * 180) / np.pi) for i in range(params.iter)]
plt.plot(alpha_x,np.abs(Y_ifft[params.M-1:params.N]))
plt.vlines(x =round((params.phi_signal * 180) / np.pi), ymin = min(np.abs(Y_ifft[params.M-1:params.N])), ymax = max(np.abs(Y_ifft[params.M-1:params.N])),
           colors = 'purple',
           label = 'vline_multiple - full height')
plt.title("График сигнального отклика системы формирования веера ПКН (угол прихода сигнала = "+str(round((params.phi_signal * 180) / np.pi))+")\n Параметры: Количество ПЭ N = "+str(params.N)+", частота = "+str(params.frequency_FFT)+"\n Количество направлений фазирования = "+str(params.N-params.M+1))
plt.show()
