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

N_patients=int(1e6) # Number of studied patients

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

print("Baseline TCP (only RT, no heterogeneity): ",baseline)
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


# Definition of vectors to save the results
T=np.linspace(39.5,42.5,4*5)     # Vector of mean temperature intervals 
T_mean=np.zeros(len(T)-1)        # Mean temperature vector
TCP_mean=np.zeros(len(T)-1)      # Mean TCP vector without heterogeneity  
TCP_mean_5=np.zeros(len(T)-1)    # Mean TCP vector without heterogeneity (c.v. = 5%)
TCP_mean_50=np.zeros(len(T)-1)   # Mean TCP vector without heterogeneity (c.v. = 50%)
count=np.zeros(len(T)-1)         # Number of patients in a certain mean temperature interval


# Loop to go through all mean temperature intervals
for i in range(0,len(T)-1):
    T_mean[i]=(T[i]+T[i+1])/2
    # Loop to calculate the mean TCP value
    for n in range(0,N_patients):
        # Condition to check if the simulated patient is in a certain mean temperature interval
        if Temp[n]>T[i] and Temp[n]<=T[i+1]:
            # Add TCP of simulated patient to, eventually, calculate the mean TCP
            TCP_mean[i]+=TCP[n]
            TCP_mean_5[i]+=TCP_5[n]
            TCP_mean_50[i]+=TCP_50[n]
            
            # Increase counter to, eventually, calculate the mean TCP
            count[i]+=1
    
    # Eliminate intervals without patients to avoid diving by 0
    if count[i]>0:
        TCP_mean[i]=TCP_mean[i]/count[i]
        TCP_mean_5[i]=TCP_mean_5[i]/count[i]
        TCP_mean_50[i]=TCP_mean_50[i]/count[i]

        
# Eliminate intervals without patients
T_mean = T_mean[count!=0]
TCP_mean = TCP_mean[count!=0]
TCP_mean_5 = TCP_mean_5[count!=0]
TCP_mean_50 = TCP_mean_50[count!=0]


# Calculation of the mean TCP baseline (only RT) 
TCP_mean_baseline_5=0
TCP_mean_baseline_50=0


for n in range(0,N_patients):
    N=volume/(2e-9)
    N=N*exp(-N_fract*a_5[n]*2-N_fract*a_5[n]*2*2/17.9)
    TCP_mean_baseline_5+=exp(-N)/N_patients
    
    N=volume/(2e-9)
    N=N*exp(-N_fract*a_50[n]*2-N_fract*a_50[n]*2*2/17.9)
    TCP_mean_baseline_50+=exp(-N)/N_patients
    
    
print("Baseline mean TCP (only RT, c.v. = 5%): ",TCP_mean_baseline_5)
print()

print("Baseline mean TCP (only RT, c.v. = 50%): ",TCP_mean_baseline_50)
print()
    

# Save results
np.savetxt('results/TCP_mean_T_het.dat',(np.c_[T_mean, TCP_mean, TCP_mean_5, TCP_mean_50]), header="Temp TCP TCP_5 TCP_50", comments="")
