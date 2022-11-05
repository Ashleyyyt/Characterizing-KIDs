#imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv as I0
from scipy.special import kv as K0
import scipy.constants as const

#Define function
def coth(x):
    return np.cosh(x)/np.sinh(x)

def csch(x):
    return 1/np.sinh(x)

#Define constants 
TC = 1.5 
Delta_0 = (3.5*const.Boltzmann*TC)/2
sigma_n = 6.0e7          # Normal stae conductvity if superconducting film
Thick = 20e-9            # Thickness of superconducting fil           
n = 18e28     
f = 500e6
w = 2 * np.pi * f
me = const.m_e
miu_0 = 4*np.pi*10**-7


#Varying range of temperature
T = np.linspace(0, 1.5, num=500)
frac = T/TC
ns = n*(1-(frac)**4)
nn =  n - ns

#Plot
ratio_ns = ns/n
ratio_nn = nn/n
plt.plot(frac, ratio_ns, label="$n_s$", color="red")
plt.plot(frac, ratio_nn, label="$n_n$", color="black")
plt.ylabel("$n_x/n$")
plt.xlabel("$T/T_C$")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True)
plt.rcParams['figure.dpi'] = 1000
