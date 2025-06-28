import numpy as np
import matplotlib.pyplot as plt
import ctypes
import os
import time as tm


lib_dir = os.path.dirname(os.path.abspath(__file__)) + "/UM_Blochsim.so"
lib = ctypes.CDLL(lib_dir)
    
lib.c_blochsim_eul.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_double]
lib.c_blochsim_rk4.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_double]
lib.c_blochsim_ljn.argtypes = [
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int), ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, 
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double  
]
lib.c_blochsim_ljn_dyntime.argtypes = [
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int), ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double, 
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double 
]

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

    return M

def blochsim_ljn(B, s, M_start, R1f_app, R2f_app, R1s_app, dt, ks, kf, f, M0_f, M0_s, crusher_inds=np.array([]), absorp=0.0, s_sat=0.0, plot=False, dsample=1, timer=False):
    if not isinstance(B, np.ndarray):
        raise TypeError("Input must be a NumPy array.")
    if B.dtype != np.float64:
        raise TypeError("B must have dtype np.float64.")
    if B.ndim != 2:
        raise ValueError("B must be 2 dimensional")
    if B.shape[1] != 3:
        raise ValueError("The first dimension of B must have length 3 (Simulating 3D Space), got", B.shape[1])
    if len(s.shape) != 1:
        raise ValueError("The first dimension of s must have length 1, got", s.shape[1])
    if B.shape[0] != s.shape[0]:
        raise ValueError("B and s must have the same length in dimesnion 0, got lengths: B: ", B.shape[0], " and s: ", s.shape[0])
    if len(M_start) != 4:
        raise ValueError("The starting magnetization must have dimension (1, 4), 3 flow components and one semisoid component, got: ", M_start.shape)
    
    ntime = np.shape(B)[0]
    M = np.zeros((ntime, 4), dtype=np.float64)

    B_ptr = B.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    M_ptr = M.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    M_start_ptr = M_start.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    s_ptr = s.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    
    # Crusher Pointer
    if crusher_inds.size == 0:
        crush_ptr = ctypes.POINTER(ctypes.c_int)()
        num_crushers = 0
    else:
        crush_ptr = crusher_inds.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        num_crushers = np.size(crusher_inds)

    start_t = tm.time()
    lib.c_blochsim_ljn(M_ptr, B_ptr, s_ptr, M_start_ptr, crush_ptr, num_crushers, ntime, np.float64(dt), f, ks, kf, R1f_app, R2f_app, R1s_app, M0_s, M0_f, absorp, s_sat)
    end_t = tm.time()

    if timer:
        print(f"Time: {end_t - start_t}")
    
    if plot:
        plotter_4D(M, dt, dsample)

    return M

def blochsim_ljn_dyntime(B, s, M_start, time_vec, R1f_app, R2f_app, R1s_app, ks, kf, f, M0_f, M0_s, crusher_inds=np.array([]), absorp=0.0, s_sat=0.0, plot=False, timer=False):
    if not isinstance(B, np.ndarray):
        raise TypeError("Input must be a NumPy array.")
    if B.ndim != 2:
        raise ValueError("B must be 2 dimensional")
    if B.shape[1] != 3:
        raise ValueError("The first dimension of B must have length 3 (Simulating 3D Space), got", B.shape[1])
    if len(s.shape) != 1:
        raise ValueError("The first dimension of s must have length 1, got", s.shape[1])
    if (B.shape[0] != s.shape[0]) or (B.shape[0] != time_vec.shape[0]):
        raise ValueError("B, s, and time_vec must have the same length in dimesnion 0, got lengths: B: ", B.shape[0], " s: ", s.shape[0], " time_vec: ", time_vec.shape[0])
    if len(M_start) != 4:
        raise ValueError("The starting magnetization must have dimension (1, 4), 3 flow components and one semisoid component, got: ", M_start.shape)
    
    ntime = np.shape(B)[0]
    M = np.zeros((ntime, 4), dtype=np.float64)

    B_ptr = B.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    M_ptr = M.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    M_start_ptr = M_start.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    s_ptr = s.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    time_ptr = time_vec.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    
    # Crusher Pointer
    if crusher_inds.size == 0:
        crush_ptr = ctypes.POINTER(ctypes.c_int)()
        num_crushers = 0
    else:
        crush_ptr = crusher_inds.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        num_crushers = np.size(crusher_inds)

    start_t = tm.time()
    lib.c_blochsim_ljn_dyntime(M_ptr, B_ptr, s_ptr, M_start_ptr, time_ptr, crush_ptr, num_crushers, ntime, f, ks, kf, R1f_app, R2f_app, R1s_app, M0_s, M0_f, absorp, s_sat)
    end_t = tm.time()

    if timer:
        print(f"Time: {end_t - start_t}")
    
    if plot:
        plt.plot(time_vec, M[:, 0], label = 'x free')
        plt.plot(time_vec, M[:, 1], label = 'y free')
        plt.plot(time_vec, M[:, 2], label = 'z free')
        plt.plot(time_vec, M[:, 3], label = 'semisolid')
        plt.xlabel("Time [ms]")
        plt.ylabel("Magnetization")
        plt.title("Simulated Magnetization [T/volume]")
        plt.legend()
        plt.show()

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

def plotter_4D(M, dt, dsample):
        time = np.arange(M.shape[0]) * dt
        time = time[::dsample]
        plt.plot(time[ : ], M[::dsample, 0], label = 'x free')
        plt.plot(time[ : ], M[::dsample, 1], label = 'y free')
        plt.plot(time[ : ], M[::dsample, 2], label = 'z free')
        plt.plot(time[ : ], M[::dsample, 3], label = 'semisoid')
        plt.xlabel("Time [ms]")
        plt.ylabel("Magnetization")
        plt.title("Simulated Magnetization [T/volume]")
        plt.legend()
        plt.show()
