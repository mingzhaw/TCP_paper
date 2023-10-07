import numpy as np
from math import exp, log
import random
from scipy import stats

# Definition of model parameters
alpha_37 = 0.386
alpha_41 = 1.73*alpha_37

beta_37 = alpha_37/17.9
beta_41 = 0.41*beta_37

mu = 0.5/60 #in min^-1
S = 392.08
H = 147907.8


# Definition of the extended LQ model
def LQ_T(T, t_int):
    alpha=alpha_37*exp((T-37)*log(alpha_41/alpha_37)*exp(-mu*t_int)/(41-37))
    beta=beta_37*exp((T-37)*log(beta_41/beta_37)*exp(-mu*t_int)/(41-37))
    
    return alpha, beta


# Definition of treatment conditions
N_fract = 30     # Number of RT sessions
D = 2            # Dose per fraction
N_HT = 4         # Number of HT sessions
t_int = 30       # Time interval in min
T_min=39         # Minimum temperature
T_max=43         # Maximum temperature

N_patients=10000 # Number of studied patients

volume=180       # Tumor volume in cm^3

# Definition of vectors to save the results
Temp = np.zeros(N_patients)           # Mean temperature vector
Temp_min = np.ones(N_patients)*T_max  # Minimum temperature vector
Temp_max = np.zeros(N_patients)       # Maximum temperature vector
TCP = np.zeros(N_patients)            # TCP vector (without direct HT cell killing) 
TCP_kill = np.zeros(N_patients)       # TCP vector (with direct HT cell killing)


# Calculation of the baseline TCP (only RT)
N=volume/(2e-9)
N=N*exp(-N_fract*alpha_37*2-N_fract*beta_37*2*2)
baseline = np.exp(-N)

print("Baseline TCP (only RT): ",baseline)
print()


random.seed(10)   # Random generator seed to obtaion always the same results

# Loop to calculate the TCP values
for n in range(0, N_patients):
    # Definition of initial number of cells
    N=volume/(2e-9)
    N_kill=volume/(2e-9)
    
    for i in range(0, N_fract):
        if i<N_HT:
            # Random selection of temperature
            temp_achi = random.uniform(T_min, T_max)
            
            # Record mean temperature
            Temp[n] += temp_achi/N_HT
            
            # Record maximum temperature
            if temp_achi>Temp_max[n]:
                Temp_max[n]=temp_achi
            
            # Record minimum temperature
            if temp_achi<Temp_min[n]:
                Temp_min[n]=temp_achi
            
            # Calculation of direct HT cell killing
            HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*3600)
            
            # Calculation of the LQ paramaters after HT session
            alpha_T, beta_T = LQ_T(temp_achi, t_int)
            
            # Calculation of survival after HT+RT fraction
            N=N*exp(-alpha_T*D-beta_T*D*D)
            N_kill=N_kill*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
        else:
            # Calculation of survival after RT fraction
            N=N*exp(-alpha_37*D-beta_37*D*D)
            N_kill=N_kill*exp(-alpha_37*D-beta_37*D*D)            

    
    # Calculation of TCP after the HT+RT treatment
    TCP[n]=exp(-N)   
    TCP_kill[n]=exp(-N_kill)


# Save results
np.savetxt('results/TCP_T.dat',(np.c_[Temp_min, Temp, Temp_max, TCP, TCP_kill]), header="Temp_min Temp Temp_max TCP TCP_kill", comments="")

# Calculation of correlations

# No killing
print("No killing")
res = stats.spearmanr(Temp_min, TCP)
spearman_T = res.correlation
spearman_p_T = res.pvalue
    
kendall_T, kendall_p_T = stats.kendalltau(Temp_min, TCP)

print("Min T Spearman's rho: ", spearman_T)
print("Min T Kendall's tau: ", kendall_T)
print("p-value: ",spearman_p_T," ",kendall_p_T)
print()

res = stats.spearmanr(Temp, TCP)
spearman_T = res.correlation
spearman_p_T = res.pvalue
    
kendall_T, kendall_p_T = stats.kendalltau(Temp, TCP)

print("Mean T Spearman's rho: ", spearman_T)
print("Mean T Kendall's tau: ", kendall_T)
print("p-value: ",spearman_p_T," ",kendall_p_T)
print()

res = stats.spearmanr(Temp_max, TCP)
spearman_T = res.correlation
spearman_p_T = res.pvalue
    
kendall_T, kendall_p_T = stats.kendalltau(Temp_max, TCP)

print("Max T Spearman's rho: ", spearman_T)
print("Max T Kendall's tau: ", kendall_T)
print("p-value: ",spearman_p_T," ",kendall_p_T)
print()

# Killing
print("Killing")
res = stats.spearmanr(Temp_min, TCP_kill)
spearman_T = res.correlation
spearman_p_T = res.pvalue
    
kendall_T, kendall_p_T = stats.kendalltau(Temp_min, TCP_kill)

print("Min T Spearman's rho: ", spearman_T)
print("Min T Kendall's tau: ", kendall_T)
print("p-value: ",spearman_p_T," ",kendall_p_T)
print()

res = stats.spearmanr(Temp, TCP_kill)
spearman_T = res.correlation
spearman_p_T = res.pvalue
    
kendall_T, kendall_p_T = stats.kendalltau(Temp, TCP_kill)

print("Mean T Spearman's rho: ", spearman_T)
print("Mean T Kendall's tau: ", kendall_T)
print("p-value: ",spearman_p_T," ",kendall_p_T)
print()

res = stats.spearmanr(Temp_max, TCP_kill)
spearman_T = res.correlation
spearman_p_T = res.pvalue
    
kendall_T, kendall_p_T = stats.kendalltau(Temp_max, TCP_kill)

print("Max T Spearman's rho: ", spearman_T)
print("Max T Kendall's tau: ", kendall_T)
print("p-value: ",spearman_p_T," ",kendall_p_T)
print()
