#imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.special import iv as I0
from scipy.special import kv as K0

#Define Global Variables
L_geo = 55.6e-9
Z0 = 50.0
F0_base = 0.95e9        #At lowest Temp
squares= 27223
c_couple = 1.5e-14
TC = 1.5
Delta_0 = (3.5*const.Boltzmann*TC)/2
sigma_n = 6.0e7          # Normal stae conductvity if superconducting film
Thick = 20e-9            # Thickness of superconducting fil                
w = 2 * np.pi * F0_base
me = const.m_e
miu_0 = 4*np.pi*10**-7
pi = np.pi

#Main code
def main():
    #Define temperature range with step 0.01K
    step = 0.1
    temp = np.arange(0.2, 0.3   , step)
    
    #Find sigma1 and sigma 2 and Lint
    sigma1, sigma2 = find_sigma1_sigma2(sigma_n ,Thick, TC, Delta_0, w, temp)
    Lint = find_Lint_square(Thick, w, sigma2) * squares
    
    #Find lk
    Lk = find_lk(Thick, w, sigma2)
    
    #Find Res
    sigma12Ratio = sigma1/sigma2
    Res = Lk*w*sigma12Ratio *squares
    
    #IDC for Lowest Temp (0.2K)
    Ltot_lowest = Lint[0] + L_geo
    IDC = find_IDC(w, Ltot_lowest, c_couple)
    
    #Find S21
    Sweep_points = 20000
    BW = 5e6 
    I_raw = np.zeros((Sweep_points, len(temp)), dtype="float")
    Q_raw = np.copy(I_raw)
    Phase = np.copy(Q_raw)
    S21_Volt = np.copy(I_raw)
    for i in range(0, len(Lint)):
        Sweep, S21_Volt[:,i], Phase[:,i], I_raw[:,i], Q_raw[:,i],_,_,_,_,_ = Capacitive_Res_Sim(F0_base, c_couple, Z0, L_geo, Lint[i], Res[i], BW, Sweep_points, IDC)
        plt.plot(Sweep/1e9, S21_Volt[:,i], label=str("{:.2f}".format(temp[i])))
    
    #Graph labels and title
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, title="Temperature / K")
    plt.xlabel('Frequency / GHz', fontsize=13)
    
    plt.ylabel('Output Amplitude', fontsize=13)
    #plt.title("S21 Amplitude For Varying Temperatures")
    plt.xlim(0.9490, 0.9505)
    plt.locator_params(nbins=6)
    plt.savefig("S21 Plot with Resistance")
    plt.rcParams['figure.dpi'] = 300
    plt.figure()
    
    #Q vs I plots
    for i in range(0, len(Lint)):
        plt.plot(I_raw[:,i], Q_raw[:,i], linewidth=1,label=str("{:.2f}".format(temp[i])))
    
    #Minimum S21 at lowest temp
    S21_Base = min(S21_Volt[:,0])
    I_Base = np.zeros(len(temp), dtype="float")
    Q_Base = np.copy(I_Base)
    
    #Obtain F0_base and I and Q values for Lowest Temp
    for i in range(0, len(S21_Volt[:,0])):
        if S21_Base == S21_Volt[i,0]:
            F0_Base = Sweep[i]
    
    #Plot I and Q values at F0_Base
    for i in range(0, len(temp)):
        for j in range(0, len(Sweep)):
            if F0_Base == Sweep[j]: 
                I_Base[i] = I_raw[j,i]
                Q_Base[i] = Q_raw[j,i]
                plt.plot(I_Base[i], Q_Base[i], markersize=4, marker="x", color='black')
    
    #labels
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, title="Temperature / K")
    plt.xlabel('I / V', fontsize=13)
    plt.ylabel('Q / V', fontsize=13)
    #plt.title("Q vs I Plot for Varying Temperature")
    plt.savefig("Q vs I plot for varying temp")
    plt.figure()        
    
    #Finding F0 for the different Temperatures
    F0 = np.zeros(len(temp))
    for i in range(0, len(temp)):
        S21_min = min(S21_Volt[:,i])
        for j in range(0, len(Sweep)):
            if S21_min == S21_Volt[j,i]:
                F0[i] = Sweep[j]
    
    #Plotting F0 vs Temp
    plt.plot(temp, F0/1e9, color='k', linewidth="1", label="Minimum Of S21")
    plt.xlabel('Temperature / K', fontsize=13)
    plt.ylabel('F0 / GHz', fontsize=13)
    plt.rcParams['figure.dpi'] = 300
    #plt.title("F0 vs Temperature")
    
    #Finding dI/dF and dQ/dF for lowest temperature
    #Using numerical derivatives
    step = abs((Sweep[0]-Sweep[-1])/Sweep_points)
    for i in range(0, len(Sweep)):
        if Sweep[i] == F0_Base:
            didf = (I_raw[i+1,0] - I_raw[i-1,0])/(2*step)
            dqdf = (Q_raw[i+1,0] - Q_raw[i-1,0])/(2*step)

    #Use Magic Formula
    di = np.zeros(len(temp))
    dq = np.copy(di)
    di = abs(I_Base - I_Base[0])
    dq = abs(Q_Base - Q_Base[0])
    dF0 = Magic_Formula(di, dq, didf, dqdf)
    
    #Find F0 for different temp
    F0_Magic = F0_Base - abs(dF0)
    plt.plot(temp, F0_Magic/1e9, label="dF0 Formula")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True)
    plt.ticklabel_format(useOffset=False)
    plt.rcParams['figure.dpi'] = 1000
    plt.xlim(0.20, 0.22)
    plt.ylim(0.949980, 0.95)
    plt.savefig("Magic Formula plot")
    
#KID Simulating Function
def Capacitive_Res_Sim(F0, C_couple, Z0, L_geo, L_int, Res, Sweep_BW, Sweep_points, Capacitance):
    """ Help file here"""    
    j=complex(0,1)    
    Cc=C_couple
    F_min=F0-(Sweep_BW/2.0)
    F_max=F0+(Sweep_BW/2.0)
    Sweep=np.linspace(F_min, F_max, Sweep_points)
    W=Sweep*2.0*pi
    W0=2.0*pi*F0
    L=L_geo+L_int
    C=Capacitance
    Zres= 1.0/((1./((j*W*L)+Res))+(j*W*C))   # Impedance of resonator section    
    Zc=1.0/(j*W*Cc)  #impedance of coupler
    ZT=Zres+Zc
    YT=1.0/ZT
    S21 = 2.0/(2.0+(YT*Z0))
    I_raw=S21.real
    Q_raw=S21.imag
    shift=((1.0-min(I_raw))/2.0)+min(I_raw)
    I_cent=I_raw-shift
    Q_cent=Q_raw
    Phase=Atan(abs(Q_cent/I_cent))
    QU=(W0*L)/Res
    QL=(C*2)/(W0*(Cc**2)*Z0)
    S21_Volt=abs(S21) 
    I_offset=shift
    return  (Sweep, S21_Volt, Phase, I_raw, Q_raw, I_cent, Q_cent, QU, QL, I_offset)

#Function to find sigma1 and sigma2
def find_sigma1_sigma2(sigma_n ,Thick, TC, Delta_0, w, T):
    #An interpolation formula for delta_T 
    delta_T = Delta_0*np.tanh(1.74*np.sqrt((TC/T)-1))
    
    #Define constants to simplify eqn
    multiplying_constant = delta_T/(const.hbar * w)
    e_const_1 = - Delta_0/(const.Boltzmann*T)
    e_const_2 = (const.hbar*w)/(2*const.Boltzmann*T)

    #Parts of the sigma1 Ratio
    A = 2*multiplying_constant
    B = np.exp(e_const_1)
    C = K0(0, e_const_2)
    D = 2*(np.sinh(e_const_2))

    #Find Sigma 1 and Sigma 2
    sigma1Ratio = A * B * C * D
    sigma2Ratio = np.pi*multiplying_constant*(1 - (2*np.exp(e_const_1)*np.exp(-e_const_2)*I0(0,e_const_2)))
    sigma2 = sigma2Ratio * sigma_n
    sigma1 = sigma1Ratio * sigma_n
    return sigma1, sigma2
 
def find_lk(Thick, w, sigma2):
    #Depth
    lower_fraction = miu_0*sigma2*w
    Lambda_T_MB = (1/lower_fraction)**0.5
    fraction = Thick/(2*Lambda_T_MB)
    
    #Terms for lk
    A = (miu_0*Lambda_T_MB)/4
    B = coth(fraction)
    C = fraction*(csch(fraction))**2
    
    #R vs T
    lk = A*(B+C)
    return lk
   
def find_Lint_square(Thick, w, sigma2):
    #Depth
    lower_fraction = miu_0*sigma2*w
    Lambda_T_MB = (1/lower_fraction)**0.5
    
    #Internal Inductance
    fraction = Thick/(2*Lambda_T_MB)
    L_int = (miu_0*Lambda_T_MB/2)*coth(fraction)
    return L_int

#Define coth and csch
def coth(x):
    return np.cosh(x)/np.sinh(x)

def csch(x):
    return 1/np.sinh(x)

def Atan(x):
    return np.arctan(x)

#Find IDC function
def find_IDC(w0, Ltot, Cc):
    IDC = 1/((w0**2)*Ltot) - Cc
    return IDC

def Magic_Formula(di, dq, didf, dqdf):
    return (di*didf + dq*dqdf)/(didf**2 + dqdf**2)
    
main()