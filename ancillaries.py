import numpy as np

constants = {
    'T_CMB':2.7255,
    'h':6.62607004e-34,
    'k_B':1.38064852e-23,
    'c':2.99792e8
}

units = {
    'n':1.e-9, 
    'u':1.e-6, 
    'm':1.e-3, 
    '1': 1.,
    'k':1.e3, 
    'M':1.e6, 
    'G':1.e9,
    'K_RJ': lambda x: 3.072387e4*x**2,
    'K_CMB': lambda x:  (0.017608676*x/(np.exp(0.017608676*x)-1))**2*np.exp(0.017608676*x)*3.072387e4*x**2,
    'Jysr': lambda x: 1.
}
