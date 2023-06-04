import params
from scipy.fft import fft
import numpy as np
from matplotlib import pyplot as plt

FFT = [ [0]*params.n for i in range(params.N)]
for i in range(params.N):
    y_signal = params.signals_array_FFT[i]
    FFT[i] = fft(y_signal, params.n)



FFT_signals_fk = [0]*params.N
for t in range(params.N):
    FFT_signals_fk[t] = FFT[t][params.fk_index]

fig, (ax11,ax) = plt.subplots(2,1)
ax11.plot(params.x[:600], params.signals_array_FFT[49][:600])
ax11.set_title('Сигнал')
ax11.set_xlabel('Время, с')
ax.plot(params.fk, np.abs(FFT[49]))
ax.set_title('Амплитудный спектр')
ax.set_xlabel('Частота, Гц')
fig.tight_layout()

plt.show()

t_diff_L_coeff = [0]*params.M
L_coeff = [0]*params.M
V = [0]*params.M
V_sum = [0]*params.n
W = [0]*params.iter
t_diff_L_coeff = params.t_diff_arr[params.m_index:params.m_index+params.M]
L_coeff = np.exp(1j * 2 * np.pi * params.frequency_FFT * np.asarray(t_diff_L_coeff))

for i in range(params.iter):
    V = [FFT_signals_fk[i+r]*L_coeff[r] for r in range(params.M)]
    V_sum = 0
    for j in range(params.M):
        V_sum+=V[j]
    W[i] = np.abs(V_sum)


print(W)
print(len(W))
alpha_x = [round((params.alpha[i] * 180) / np.pi) for i in range(params.iter)]

print(alpha_x)
plt.plot(alpha_x,W)
plt.vlines(x =round((params.phi_signal * 180) / np.pi), ymin = min(W), ymax = max(W),
           colors = 'purple',
           label = 'vline_multiple - full height')
plt.title("График сигнального отклика системы формирования веера ПКН (угол прихода сигнала = "+str(round((params.phi_signal * 180) / np.pi))+")\n Параметры: Количество ПЭ N = "+str(params.N)+", частота = "+str(params.frequency_FFT)+"\n Количество направлений фазирования = "+str(params.N-params.M+1))
plt.show()

