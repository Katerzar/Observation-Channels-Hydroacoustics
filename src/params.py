import numpy as np
import functions
N = 90            # число приёмных элементов
M = 30               # число приёмных элементов в рабочем пятне
c = 1500             # скорость звука в воде
frequency = 5000     # частота
a_XH = 0.7           # коэффициент направленности
freq_min = 3000      # нижний порог частоты, ГЦ
freq_max = 7500      # верхний порог частоты, ГЦ
duration = 0.1       # длительность сигнала, с
n = 256              # Число точек БПФ для свертки во временной области

k = 0
phi_signal = k*np.pi/3 # Угол прихода сигнала


freq_d = 2.5*freq_max                # Частота дискретизации f_d = 2.5f_max -? 2500 (sample_rate)
phi_diff = ((2 * np.pi) / N)         # Угловое расстояние между ПЭ
delta_time = 1 / freq_d  # 40 мкс    # шаг дискретизации
lambda_signal = c/freq_max           # Длина волны сигнала на верхней частоте
d = lambda_signal*0.5                # lambda = c/freq_max = 1500/4000 = 0.375
                                     # Расстояние между приёмными элементами
R = round((d*N)/(2*np.pi),2)         # радиус антенны

phi = [-np.pi+(np.pi / N) + phi_diff * i for i in range(N)]  # создание массива угловых положений ПЭ
x_coord = [round(R * np.cos(phi[i]),3) for i in range(N)]    # создание массива х-координат ПЭ
y_coord = [round(R * np.sin(phi[i]),3) for i in range(N)]    # создание массива у-координат ПЭ
alpha = [-np.pi + (np.pi / N) + (phi_diff*(M-1)/2) + phi_diff*i for i in range(N-M+1)]
                                     # создание массива угловых направлений рабочих пятен (N-M+1 направление)
t_diff_arr = [((-R*(np.cos(phi[m]-phi_signal))) / c) for m in range (N)]  #создание массива временных задержек для каждого приёмного элемента в зависимости от угла прихода волны
x_XH_array = [max(0, (a_XH + np.cos(phi[m]-phi_signal)) / (a_XH + 1)) for m in range(N)] # создание массива характеристической направленности в зависимости от угла прихода волны

m_index = 0 # индекс номера элемента, соответствующего началу рабочего пятна при заданном угле прихода сигнала
for i in range(N-M+1):
    if np.abs(alpha[i] - phi_signal)<0.001:
        m_index = i
        break

x = np.linspace(0, duration, int(freq_d * duration), endpoint=False) # массив моментов времени
print(len(x))

signals_array_FFT = [functions.create_signal(x, frequency, t_diff_arr[i], x_XH_array[i]) for i in range(N)]# двумерный массив сигналов для каждого ПЭ [N][len(x)]

df = freq_d / n                     # шаг дискретизации спектра
fk = df * np.linspace(0, n - 1, n)  # массив частот с n элементами и шагом df
fk_index = round(frequency/df)      # индекс частоты из массива частот спектра,наиболее близкой к заданной частоте сигнала
frequency_FFT = fk[fk_index]        # частота из массива частот спектра, наиболее близкая к заданной частоте сигнала
iter = N-M+1
