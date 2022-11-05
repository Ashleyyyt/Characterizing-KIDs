import sys
fundir = '/Users/simondoyle/OneDrive - Cardiff University/Programs/Python/Functions/Simulation_functions/'
if fundir not in sys.path:
    sys.path.append(fundir)

import numpy as np
from numpy import log10 as log10
from numpy import pi as pi
from numpy import tan as tan
from numpy import cos as cos
from numpy import sin as sin
from numpy import exp as exp
from numpy import sqrt as sqrt
from numpy import arctan as Atan
from numpy import sinh as sinh
from numpy import tanh as tanh
from matplotlib import pyplot as plt

from scipy.constants import h as h
from scipy.constants import hbar as hbar
from scipy.constants import k as kB
from scipy.constants import mu_0 as mu0
from scipy.constants import electron_volt as eV

from scipy.special import iv as I0
from scipy.special import kv as K0

plt.close('all')




BW=10000



Sample_length=0.1
Spec_dense=0.001   #x/(root_Hz)

Sigma=Spec_dense*(sqrt(BW))

Sig_F=200


ave=100

points=2.0*BW*Sample_length
points=int(points)
TS=np.linspace(0,Sample_length, points)


AC_sig=(0.5*sin(2*pi*Sig_F*TS))+0.5
DC_sig=(AC_sig*0)+1.0
NO_sig=(AC_sig*0)

noise = np.random.normal(0,Sigma,points)

Noisy_AC_sig=(0.5*sin(2*pi*Sig_F*TS))+noise
Noisy_DC_sig=(AC_sig*0)+1.0+noise
Noisy_NO_sig=(AC_sig*0)+noise




from scipy import signal


FS=points/Sample_length

for x in range (0,ave):
    
        noise = np.random.normal(0,Sigma,points)


        if x == 0:
            f, Pxx_den = signal.periodogram(noise, FS)
            
            PS_ave=Pxx_den
        
        if x != 0:
            f, Pxx_den = signal.periodogram(noise, FS)
            PS_ave=PS_ave+Pxx_den


PS_ave=(PS_ave/ave)

PS_ave=sqrt(PS_ave)




for x in range (0,ave):
    
        noise = np.random.normal(0,Sigma,points)
        Noisy_sig=(0.5*sin(2*pi*Sig_F*TS))+0.5+noise

        if x == 0:
            f, Pxx_den = signal.periodogram(Noisy_sig, FS)
            
            Sig_PS_ave=Pxx_den
        
        if x != 0:
            f, Pxx_den = signal.periodogram(Noisy_sig, FS)
            Sig_PS_ave=Sig_PS_ave+Pxx_den
            

Sig_PS_ave=(Sig_PS_ave/ave)

Sig_Spec_ave=sqrt(Sig_PS_ave)






fig1,ax1 = plt.subplots(1,1)

ax1.plot(TS, DC_sig, label='DC signal')
ax1.plot(TS, NO_sig, label='NO signal')
ax1.plot(TS, AC_sig, label='Modulated signal')
ax1.get_xaxis().get_major_formatter().set_useOffset(False)
plt.ylabel('Amplitude', fontsize=16)
plt.xlabel('Time', fontsize=16)
plt.legend(loc=1)
ax1.plot(TS, AC_sig,'b')
plt.show()


fig2,ax2 = plt.subplots(1,1)
ax2.plot(TS, Noisy_AC_sig, label='Modulated signal')
ax2.plot(TS, Noisy_DC_sig, label='DC signal')
ax2.plot(TS, Noisy_NO_sig, label='NO signal')
ax2.plot(TS, Noisy_AC_sig,'b')
ax2.get_xaxis().get_major_formatter().set_useOffset(False)
plt.ylabel('Amplitude', fontsize=16)
plt.xlabel('Time', fontsize=16)
plt.legend(loc=1)
plt.show()


fig3,ax3 = plt.subplots(1,1)
ax3.plot(f,PS_ave)
ax3.get_xaxis().get_major_formatter().set_useOffset(False)
plt.ylabel('Amplitude', fontsize=16)
plt.xlabel('freq', fontsize=16)
plt.legend(loc=1)
plt.show()


fig4,ax4 = plt.subplots(1,1)
ax4.plot(f,Sig_Spec_ave)
ax4.get_xaxis().get_major_formatter().set_useOffset(False)

plt.ylabel('Noise / Signal Amplitude', fontsize=16)
plt.xlabel('freq', fontsize=16)
plt.legend(loc=1)
plt.show()

TS=np.linspace(0,Sample_length, points)


F=np.linspace(0.,1e6, 1000)
W=2.*pi*F

R=10e6
C=1e-12
tau=R*C

out1=   1./(sqrt(1.+((W*tau)**2)))

R=1e6
C=1e-12
tau=R*C

out2=   1./(sqrt(1.+((W*tau)**2)))




fig5,ax5 = plt.subplots(1,1)
ax5.plot(F,out1, label=r"$R=1M\Omega$ C=1pF $\tau$ =1 $\mu s$")
ax5.plot(F,out2, label=r"$R=10M\Omega$ C=1pF $\tau$ =10 $\mu s$")
ax5.get_xaxis().get_major_formatter().set_useOffset(False)
plt.ylabel('Filter function', fontsize=16)
plt.xlabel('freq/Hz', fontsize=16)
plt.legend(loc=1)
plt.show()


fig11,ax11 = plt.subplots(1,1)


#ax11.plot(TS, AC_sig, label='Modulated signal No noise')
ax11.plot(TS, AC_sig+1*noise, label='Modulated signal with noise')
ax11.get_xaxis().get_major_formatter().set_useOffset(False)
plt.ylabel('Amplitude', fontsize=16)
plt.xlabel('Time', fontsize=16)
plt.legend(loc=1)
ax1.plot(TS, AC_sig,'b')
plt.show()
