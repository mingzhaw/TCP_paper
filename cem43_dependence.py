import numpy as np
from math import exp, log
import random

# Definition of model parameters
alpha_37 = 0.386
alpha_41 = 1.73*alpha_37

beta_37 = alpha_37/17.9
beta_41 = 0.41*beta_37

S = 392.08
H = 147907.8
mu = 0.5/60 #in min-1


# Definition of the extended LQ model
def LQ_T(T, t_int):
    alpha=alpha_37*exp((T-37)*log(alpha_41/alpha_37)*exp(-mu*t_int)/(41-37))
    beta=beta_37*exp((T-37)*log(beta_41/beta_37)*exp(-mu*t_int)/(41-37))
    
    return alpha, beta

# Definition of treatment conditions
N_fract = 30          # Number of RT sessions
D = 2                 # Dose per fraction
T_min=37              # Minimum temperature
T_max=43              # Minimum temperature
N_patients=int(1e6)   # Number of studied patients

volume=180       # Tumor volume in cm^3

# Definition of vectors to save the results
cem43 = np.zeros(4*N_patients)                  # Total CEM43 vector
treatments = np.zeros(4*N_patients)             # Number of HT sessions vector
TCP_10 = np.zeros(4*N_patients)                 # TCP vector (without direct HT cell killing and t_int = 10 min) 
TCP_kill_10 = np.zeros(4*N_patients)            # TCP vector (with direct HT cell killing and t_int = 10 min)  
TCP_60 = np.zeros(4*N_patients)                 # TCP vector (without direct HT cell killing and t_int = 60 min) 
TCP_kill_60 = np.zeros(4*N_patients)            # TCP vector (with direct HT cell killing and t_int = 60 min)
TCP_120 = np.zeros(4*N_patients)                # TCP vector (without direct HT cell killing and t_int = 120 min) 
TCP_kill_120 = np.zeros(4*N_patients)           # TCP vector (with direct HT cell killing and t_int = 120 min)
TCP_240 = np.zeros(4*N_patients)                # TCP vector (without direct HT cell killing and t_int = 240 min) 
TCP_kill_240 = np.zeros(4*N_patients)           # TCP vector (with direct HT cell killing and t_int = 240 min)



# Calculation of the baseline TCP (only RT)
N=volume/(2e-9)
N=N*exp(-N_fract*alpha_37*2-N_fract*beta_37*2*2)
baseline = np.exp(-N)

print("Baseline TCP (only RT): ",baseline)
print()

random.seed(10)   # Random generator seed to obtaion always the same results

# Loop to calculate the TCP values
for n in range(0, N_patients):
    # Loop to simulate diverse number of HT sessions
    for n_ht in range(3,7):
        # Definition of initial number of cells
        N_10=volume/(2e-9)
        N_kill_10=volume/(2e-9)
        N_60=volume/(2e-9)
        N_kill_60=volume/(2e-9)
        N_120=volume/(2e-9)
        N_kill_120=volume/(2e-9)
        N_240=volume/(2e-9)
        N_kill_240=volume/(2e-9)
    
        # Record number of HT sessions
        treatments[n+(n_ht-3)*N_patients] = n_ht
    
        for i in range(0,N_fract):
            if i<n_ht:
                # Random selection of temperature
                temp_achi = random.uniform(T_min, T_max)
                
                # Record total CEM43 value
                cem43[n+(n_ht-3)*N_patients] += (0.25**(43-temp_achi))*60
                
                # Calculation of direct HT cell killing
                HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*3600)

                # Calculation of survival after HT+RT fraction (t_int = 10 min)
                alpha_T, beta_T = LQ_T(temp_achi, 10)
                N_10=N_10*exp(-alpha_T*D-beta_T*D*D)
                N_kill_10=N_kill_10*exp(-alpha_T*D-beta_T*D*D)*HT_killing
                
                # Calculation of survival after HT+RT fraction (t_int = 60 min)
                alpha_T, beta_T = LQ_T(temp_achi, 60)
                N_60=N_60*exp(-alpha_T*D-beta_T*D*D)
                N_kill_60=N_kill_60*exp(-alpha_T*D-beta_T*D*D)*HT_killing
                
                # Calculation of survival after HT+RT fraction (t_int = 120 min)
                alpha_T, beta_T = LQ_T(temp_achi, 120)
                N_120=N_120*exp(-alpha_T*D-beta_T*D*D)
                N_kill_120=N_kill_120*exp(-alpha_T*D-beta_T*D*D)*HT_killing
                
                # Calculation of survival after HT+RT fraction (t_int = 240 min)
                alpha_T, beta_T = LQ_T(temp_achi, 240)
                N_240=N_240*exp(-alpha_T*D-beta_T*D*D)
                N_kill_240=N_kill_240*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            else:        
                # Calculation of survival after RT fraction
                N_10=N_10*exp(-alpha_37*D-beta_37*D*D)
                N_kill_10=N_kill_10*exp(-alpha_37*D-beta_37*D*D)
                N_60=N_60*exp(-alpha_37*D-beta_37*D*D)
                N_kill_60=N_kill_60*exp(-alpha_37*D-beta_37*D*D)
                N_120=N_120*exp(-alpha_37*D-beta_37*D*D)
                N_kill_120=N_kill_120*exp(-alpha_37*D-beta_37*D*D)
                N_240=N_240*exp(-alpha_37*D-beta_37*D*D)
                N_kill_240=N_kill_240*exp(-alpha_37*D-beta_37*D*D)

        # Calculation of TCPs after the HT+RT treatment (t_int = 10 min)
        TCP_10[n+(n_ht-3)*N_patients]=exp(-N_10)     
        TCP_kill_10[n+(n_ht-3)*N_patients]=exp(-N_kill_10)  

        # Calculation of TCPs after the HT+RT treatment (t_int = 60 min)
        TCP_60[n+(n_ht-3)*N_patients]=exp(-N_60)     
        TCP_kill_60[n+(n_ht-3)*N_patients]=exp(-N_kill_60)  

        # Calculation of TCPs after the HT+RT treatment (t_int = 120 min)
        TCP_120[n+(n_ht-3)*N_patients]=exp(-N_120)     
        TCP_kill_120[n+(n_ht-3)*N_patients]=exp(-N_kill_120)  

        # Calculation of TCPs after the HT+RT treatment (t_int = 240 min)
        TCP_240[n+(n_ht-3)*N_patients]=exp(-N_240)     
        TCP_kill_240[n+(n_ht-3)*N_patients]=exp(-N_kill_240)  
    
    

# Definition of CEM43 intervals    
x=np.linspace(0,1,11)
x=np.append(x,np.linspace(2,150,149))    
    
# Definition of vectors to save the results
TCP_mean_10 = np.zeros(4*len(x))        # Mean TCP vector (without direct HT cell killing and t_int = 10 min)
TCP_mean_60 = np.zeros(4*len(x))        # Mean TCP vector (without direct HT cell killing and t_int = 60 min)
TCP_mean_120 = np.zeros(4*len(x))       # Mean TCP vector (without direct HT cell killing and t_int = 120 min)
TCP_mean_240 = np.zeros(4*len(x))       # Mean TCP vector (without direct HT cell killing and t_int = 240 min)
TCP_mean_kill_10 = np.zeros(4*len(x))   # Mean TCP vector (with direct HT cell killing and t_int = 10 min)
TCP_mean_kill_60 = np.zeros(4*len(x))   # Mean TCP vector (with direct HT cell killing and t_int = 60 min)
TCP_mean_kill_120 = np.zeros(4*len(x))  # Mean TCP vector (with direct HT cell killing and t_int = 120 min)
TCP_mean_kill_240 = np.zeros(4*len(x))  # Mean TCP vector (with direct HT cell killing and t_int = 240 min)
n = np.zeros(4*len(x))                  # Number of patients in a certain CEM43 interval
cem = np.zeros(4*len(x))                # Mean CEM43
N_HT = np.zeros(4*len(x))               # Number of HT sessions


for n_p in range(0, N_patients):
    # Loop to cover diverse number of HT sessions
    for n_ht in range(3,7):
        # Loop to go through all CEM43 intervals
        for i in range(1,len(x)):
            # Definition of mean CEM43 value and number of HT sessions
            cem[i+(n_ht-3)*len(x)] = (x[i-1] + x[i])/2
            N_HT[i+(n_ht-3)*len(x)] = n_ht
            
            #Condition to check if the simulated patient is in a certain CEM43 interval
            if cem43[n_p+(n_ht-3)*N_patients]>=x[i-1] and cem43[n_p+(n_ht-3)*N_patients]<x[i]:
                # Add TCP of simulated patient to, eventually, calculate the mean TCP
                TCP_mean_10[i+(n_ht-3)*len(x)] += TCP_10[n_p+(n_ht-3)*N_patients]
                TCP_mean_kill_10[i+(n_ht-3)*len(x)] += TCP_kill_10[n_p+(n_ht-3)*N_patients]
                
                TCP_mean_60[i+(n_ht-3)*len(x)] += TCP_60[n_p+(n_ht-3)*N_patients]
                TCP_mean_kill_60[i+(n_ht-3)*len(x)] += TCP_kill_60[n_p+(n_ht-3)*N_patients]
                
                TCP_mean_120[i+(n_ht-3)*len(x)] += TCP_120[n_p+(n_ht-3)*N_patients]
                TCP_mean_kill_120[i+(n_ht-3)*len(x)] += TCP_kill_120[n_p+(n_ht-3)*N_patients]
                
                TCP_mean_240[i+(n_ht-3)*len(x)] += TCP_240[n_p+(n_ht-3)*N_patients]
                TCP_mean_kill_240[i+(n_ht-3)*len(x)] += TCP_kill_240[n_p+(n_ht-3)*N_patients]
                
                # Increase counter to, eventually, calculate the mean TCP
                n[i+(n_ht-3)*len(x)]+=1
     
    
        # Set mean TCP(CEM43 = 0) as baseline (only RT)
        TCP_mean_10[(n_ht-3)*len(x)] = baseline
        TCP_mean_60[(n_ht-3)*len(x)] = baseline
        TCP_mean_120[(n_ht-3)*len(x)] = baseline
        TCP_mean_240[(n_ht-3)*len(x)] = baseline
        TCP_mean_kill_10[(n_ht-3)*len(x)] = baseline
        TCP_mean_kill_60[(n_ht-3)*len(x)] = baseline
        TCP_mean_kill_120[(n_ht-3)*len(x)] = baseline
        TCP_mean_kill_240[(n_ht-3)*len(x)] = baseline
        
        N_HT[(n_ht-3)*len(x)] = n_ht
        n[(n_ht-3)*len(x)] = 1


# Eliminate intervals without patients to avoid diving by 0
TCP_mean_10 = TCP_mean_10[n!=0]
TCP_mean_kill_10 = TCP_mean_kill_10[n!=0]
TCP_mean_60 = TCP_mean_60[n!=0]
TCP_mean_kill_60 = TCP_mean_kill_60[n!=0]
TCP_mean_120 = TCP_mean_120[n!=0]
TCP_mean_kill_120 = TCP_mean_kill_120[n!=0]
TCP_mean_240 = TCP_mean_240[n!=0]
TCP_mean_kill_240 = TCP_mean_kill_240[n!=0]
cem = cem[n!=0]
N_HT = N_HT[n!=0]
n = n[n!=0]

# Loop to calculate the mean TCP value
for i in range(0,len(n)):
    TCP_mean_10[i]=TCP_mean_10[i]/n[i]
    TCP_mean_kill_10[i]=TCP_mean_kill_10[i]/n[i]
    
    TCP_mean_60[i]=TCP_mean_60[i]/n[i]
    TCP_mean_kill_60[i]=TCP_mean_kill_60[i]/n[i]
    
    TCP_mean_120[i]=TCP_mean_120[i]/n[i]
    TCP_mean_kill_120[i]=TCP_mean_kill_120[i]/n[i]
    
    TCP_mean_240[i]=TCP_mean_240[i]/n[i]
    TCP_mean_kill_240[i]=TCP_mean_kill_240[i]/n[i]
    

# Save results
np.savetxt('results/cem43.dat',(np.c_[cem, N_HT, TCP_mean_10, TCP_mean_kill_10, TCP_mean_60, TCP_mean_kill_60, TCP_mean_120, TCP_mean_kill_120, TCP_mean_240, TCP_mean_kill_240]), 
           header="cem N_HT TCP_10 TCP_kill_10 TCP_60 TCP_kill_60 TCP_120 TCP_kill_120 TCP_240 TCP_kill_240", comments="")
