import numpy as np
from matplotlib import pyplot as plt
from scipy.fft import fft, ifft
import params

def create_antenna():
    x = params.x_coord
    y = params.y_coord
    plt.scatter(x, y, s=20, marker='o', c='red')

    xx_tail = max(x) + 0.05
    yx_tail = 0
    xx_head = min(x) - 0.05
    yx_head = 0
    x_dx = -xx_head + xx_tail
    x_dy = -yx_head + yx_tail
    plt.arrow(xx_head, yx_head, x_dx, x_dy, width=0.01)

    xy_tail = 0
    yy_tail = min(y) - 0.05
    xy_head = 0
    yy_head = max(y) + 0.05
    y_dx = -xy_head + xy_tail
    y_dy = -yy_head + yy_tail
    plt.arrow(xy_head, yy_head, y_dx, y_dy, width=0.01)

    plt.axis('square')
    plt.title("Расположение N = " + str(params.N) + " приёмных элементов на КАР")
    plt.show()

def generate_sine_wave(x, freq, sample_rate, duration):

    frequencies = x * freq
    # 2pi для преобразования в радианы
    y = np.sin((2 * np.pi) * frequencies)
    return x, y


def generate_sine_wave_diff_time(freq, sample_rate, duration, time_diff, x_XH):
    x = np.linspace(0, duration, sample_rate * duration, endpoint=False)
    frequencies = (x - time_diff) * freq*x_XH
    # 2pi для преобразования в радианы
    y = np.sin((2 * np.pi) * frequencies)
    return x, y


def print_phi_angle(phi,M):
    f_phi = open('phi_angle.txt', 'w')
    for i in range(M):
        f_phi.write(str(i + 1))
        f_phi.write(": ")
        f_phi.write(str(round((phi[i] * 180) / np.pi)))
        f_phi.write("\n")
    f_phi.close()

def print_alpha_angle(alpha, M):
    f_alpha = open('alpha_angle.txt', 'w')
    for i in range(M):
        f_alpha.write(str(i + 1))
        f_alpha.write(": ")
        f_alpha.write(str(round((alpha[i] * 180) / np.pi)))
        f_alpha.write("\n")
    f_alpha.close()


def print_coord(x_coord, y_coord,M):
    f_coord = open('coord.txt', 'w')
    for i in range(M):
        f_coord.write(str(i + 1))
        f_coord.write(": ")
        f_coord.write(str(x_coord[i]))
        f_coord.write("\t")
        f_coord.write(str(y_coord[i]))
        f_coord.write("\n")
    f_coord.close()



def print_XH(x_XH,N):
    f_XH = open('XH.txt', 'w')
    for m in range(N):
       f_XH.write(str(m + 1))
       f_XH.write(": ")
       f_XH.write(str(x_XH[m]))
       f_XH.write("\n")
    f_XH.close()





def print_L_coef(L_coeff):
    f_L_coeff = open('L_coeff.txt', 'w')
    for i in range(M):
        f_L_coeff.write(str(i + 1))
        f_L_coeff.write(": ")
        f_L_coeff.write(str(L_coeff[i]))
        f_L_coeff.write("\n")
    f_L_coeff.close()

def print_t_diff_arr(t_diff_arr,N,phi,phi_s,R):
    f_t_diff_arr = open('t_diff_arr.txt', 'w')
    f_t_diff_arr.write("Угол прихода сигнала: ")
    f_t_diff_arr.write(str(round(phi_s * 180) / np.pi))
    f_t_diff_arr.write("\n")

    for m in range(N):
        f_t_diff_arr.write(str(m + 1))
        f_t_diff_arr.write(": ")
        f_t_diff_arr.write(str(t_diff_arr[m]*1000))
        f_t_diff_arr.write("\t")
        f_t_diff_arr.write("phi[")
        f_t_diff_arr.write(str(m + 1))
        f_t_diff_arr.write("] = ")
        f_t_diff_arr.write(str(round((phi[m] * 180) / np.pi)))
        f_t_diff_arr.write("\n")

    f_t_diff_arr.close()

def create_sine_wave(x,freq):
    y = np.sin((2 * np.pi) * x * freq) #+ np.cos((100 * np.pi) * x * freq)
    return y

def create_signal(x,freq,t_diff,x_XH):
    y = []
    for i in range(len(x)):
        time = x[i]-t_diff
        #if time < 0:
        #    time = 0
        y.append(create_sine_wave(time,freq)*x_XH)
    #noise = np.random.normal(0, 5, len(x))
    noise = [0]*len(x)
    y_res = y + noise
    return np.array(y_res)

def create_t_diff_array(N,phi_s,phi,c,R):

    t_diff_arr = [((-R*(np.cos(phi[m]-phi_s))) / c) for m in range (N)]
    return t_diff_arr



def create_signals_array(x, frequency, t_diff_arr, x_XH_array, N):
    signals_array = [create_signal(x, frequency, t_diff_arr[i], x_XH_array[i]) for i in range(N)]
    return signals_array



def create_FFT_plots(n,freq_d,signals,M):
    df = freq_d / n
    fk = df * np.linspace(0, n - 1, n)  # массив частот с n элементами и шагом df
    yfs = create_sine_wave(x, frequency) # создание волны sin(2*pi*f*t)
    y_signal2 = signals[1] #создание сигнала на 2-м приёмнике
    Fs = fft(yfs, n) #БПФ для синусовой волны
    Fs_signal2 = fft(y_signal2, n) #БПФ для сигнала на 2-м приёмнике
    plt.plot(fk, np.abs(Fs))
    plt.plot(fk, np.abs(Fs_signal2))
    plt.show()
    plt.clf()

    for i in range(M):
        y_signal = signals[i]
        Fs = fft(y_signal, n)
        plt.plot(fk, np.abs(Fs))
    plt.show()
