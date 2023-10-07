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
def LQ_T(T, t_int, alpha_37, alpha_41, beta_37, beta_41):
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
Temp = np.zeros(N_patients)                             # Mean temperature vector
Temp_min = np.ones(N_patients)*T_max                    # Minimum temperature vector
Temp_max = np.zeros(N_patients)                         # Maximum temperature vector
TCP = np.zeros(N_patients)                              # TCP vector without heterogeneity 
TCP_5 = np.zeros(N_patients)                            # TCP vector without heterogeneity (c.v. = 5%)
TCP_50 = np.zeros(N_patients)                           # TCP vector without heterogeneity (c.v. = 50%)
a_5 = np.random.normal(0.386, 0.386*0.05,N_patients)    # Diverse alpha_37 values vector (c.v. = 5%)
a_50 = np.random.normal(0.386, 0.386*0.5,N_patients)    # Diverse alpha_37 values vector (c.v. = 50%)
V = np.random.normal(201.6, 119.7,N_patients)           # Diverse tumor volume vector



# Calculation of the baseline TCP (only RT)
N=volume/(2e-9)
N=N*exp(-N_fract*alpha_37*2-N_fract*beta_37*2*2)
baseline = np.exp(-N)

print("Baseline TCP (only RT): ",baseline)
print()


for j in range(0,N_patients):
    # Check that the radiosensitivity is within the defined limits
    while a_5[j]<0.2 or a_5[j]>0.74:
        a_5[j] = np.random.normal(0.386, 0.386*0.05)
        
    while a_50[j]<0.2 or a_50[j]>0.74:
        a_50[j] = np.random.normal(0.386, 0.386*0.5)
    
    # Check that the voulme is not negative
    while V[j]<0:
        V[j] = np.random.normal(201.6, 119.7)
        

random.seed(10)   # Random generator seed to obtaion always the same results


# Loop to calculate the TCP values
for n in range(0, N_patients):
    # Definition of initial number of cells
    N=volume/(2e-9)
    N_5=V[n]/(2e-9)
    N_50=V[n]/(2e-9)
    
    for i in range(0,N_fract):
        if i<N_HT:
            # Random selection of temperature
            temp_achi = random.uniform(T_min, T_max)
            
            # Record mean temperature
            Temp[n] += temp_achi/4
            
            # Calculation of direct HT cell killing
            HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*3600)
            
            # Calculation of survival after HT+RT fraction without heterogeneity
            alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.386, 1.73*0.386, 0.386/17.9, 0.41*0.386/17.9)
            N=N*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
            # Calculation of survival after HT+RT fraction with heterogeneity (c.v. = 5%)
            alpha_T, beta_T = LQ_T(temp_achi, t_int, a_5[n], 1.73*a_5[n], a_5[n]/17.9, 0.41*a_5[n]/17.9)
            N_5=N_5*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
            # Calculation of survival after HT+RT fraction with heterogeneity (c.v. = 50%)
            alpha_T, beta_T = LQ_T(temp_achi, t_int, a_50[n], 1.73*a_50[n], a_50[n]/17.9, 0.41*a_50[n]/17.9)
            N_50=N_50*exp(-alpha_T*D-beta_T*D*D)*HT_killing
        else:
            # Calculation of survival after RT fraction
            N=N*exp(-0.386*D-beta_37*D*D)
            N_5=N_5*exp(-a_5[n]*D-a_5[n]*D*D/17.9)
            N_50=N_50*exp(-a_50[n]*D-a_50[n]*D*D/17.9)
            

    # Calculation of TCP after the HT+RT treatment
    TCP[n]=exp(-N)   
    TCP_5[n]=exp(-N_5)
    TCP_50[n]=exp(-N_50)


# Save results    
np.savetxt('results/TCP_T_het.dat',(np.c_[Temp, TCP, TCP_5, TCP_50]), header="Temp TCP TCP_5 TCP_50", comments="")


# Calculation of correlations
res = stats.spearmanr(Temp, TCP)
spearman_T = res.correlation
spearman_p_T = res.pvalue
    
kendall_T, kendall_p_T = stats.kendalltau(Temp, TCP)

print("No heterogeneity")
print("Mean T Spearman's rho: ", spearman_T)
print("Mean T Kendall's tau: ", kendall_T)
print("p-value: ",spearman_p_T," ",kendall_p_T)
print()

res = stats.spearmanr(Temp, TCP_5)
spearman_T = res.correlation
spearman_p_T = res.pvalue
    
kendall_T, kendall_p_T = stats.kendalltau(Temp, TCP_5)

print("cv = 5%")
print("Mean T Spearman's rho: ", spearman_T)
print("Mean T Kendall's tau: ", kendall_T)
print("p-value: ",spearman_p_T," ",kendall_p_T)
print()

res = stats.spearmanr(Temp, TCP_50)
spearman_T = res.correlation
spearman_p_T = res.pvalue
    
kendall_T, kendall_p_T = stats.kendalltau(Temp, TCP_50)

print("cv = 50%")
print("Mean T Spearman's rho: ", spearman_T)
print("Mean T Kendall's tau: ", kendall_T)
print("p-value: ",spearman_p_T," ",kendall_p_T)
print()
