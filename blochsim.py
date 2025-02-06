import numpy as np
import matplotlib.pyplot as plt
import ctypes
import os
import time as tm

lib_path = os.path.abspath("./UM_Blochsim/UM_Blochsim.so")
lib = ctypes.CDLL(lib_path)

lib.c_blochsim_eul.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_double]
lib.c_blochsim_rk4.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_double]

def blochsim_eul(B, T1, T2, dt, plot=False, dsample=1, timer=False):
    if not isinstance(B, np.ndarray):
        raise TypeError("Input must be a NumPy array.")
    if B.dtype != np.float64:
        raise TypeError("B must have dtype np.float64.")
    if B.ndim != 2:
        raise ValueError("B must be 2 dimensional")
    if B.shape[1] != 3:
        raise ValueError("The first dimension of B must have length 3 (Simulating 3D Space), got", B.shape[1])
    
    ntime = np.shape(B)[0]
    M = np.zeros(np.shape(B), dtype=np.float64)

    B_ptr = B.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    M_ptr = M.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    start_t = tm.time()
    lib.c_blochsim_eul(B_ptr, M_ptr, T1, T2, ntime, np.float64(dt))
    end_t = tm.time()

    if timer:
        print(f"Time: {end_t - start_t}")

    if plot:
        plotter(M, dt, dsample)

    return M


def blochsim_rk4(B, T1, T2, dt, plot=False, dsample=1, timer=False):
    if not isinstance(B, np.ndarray):
        raise TypeError("Input must be a NumPy array.")
    if B.dtype != np.float64:
        raise TypeError("B must have dtype np.float64.")
    if B.ndim != 2:
        raise ValueError("B must be 2 dimensional")
    if B.shape[1] != 3:
        raise ValueError("The first dimension of B must have length 3 (Simulating 3D Space), got", B.shape[1])
    
    ntime = np.shape(B)[0]
    M = np.zeros((ntime, 3), dtype=np.float64)

    B_ptr = B.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    M_ptr = M.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    start_t = tm.time()
    lib.c_blochsim_rk4(B_ptr, M_ptr, T1, T2, ntime, np.float64(dt))
    end_t = tm.time()

    if timer:
        print(f"Time: {end_t - start_t}")
    
    if plot:
        plotter(M, dt, dsample)

    print(M.shape)

    return M


def plotter(M, dt, dsample):
        time = np.arange(M.shape[0]) * dt
        time = time[::dsample]
        plt.plot(time[ : ], M[::dsample, 0], label = 'x')
        plt.plot(time[ : ], M[::dsample, 1], label = 'y')
        plt.plot(time[ : ], M[::dsample, 2], label = 'z')
        plt.xlabel("Time [ms]")
        plt.ylabel("Magnetization")
        plt.title("Simulated Magnetization [T/volume]")
        plt.legend()
        plt.show()
