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


# Definition of the extended LQ model
def LQ_T(T, t_int, mu):
    alpha=alpha_37*exp((T-37)*log(alpha_41/alpha_37)*exp(-mu*abs(t_int))/(41-37))
    beta=beta_37*exp((T-37)*log(beta_41/beta_37)*exp(-mu*abs(t_int))/(41-37))
    
    return alpha, beta


# Define vector for studied time intervals
t_int = np.linspace(-4,4,19)

# Definition of vectors to save the results
T = np.zeros(3*len(t_int))               # Temperature vector
interval = np.zeros(3*len(t_int))        # Time interval vector
survival = np.zeros(3*len(t_int))        # Survival vector with mu = 0.027 h^-1
survival_high = np.zeros(3*len(t_int))   # Survival vector with mu = 0.5 h^-1


# Calculation of the baseline survival (only RT)
baseline = exp(-alpha_37*2-beta_37*2*2)
print('Baseline survival (only RT): ',baseline)
print()


# Loop to study different time intervals
for i in range(0,len(t_int)):
    # Simulate survival at diverse temperatures with mu = 0.027 h^-1
    temp_achi=39.
    HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*60*60)
    alpha_T, beta_T = LQ_T(temp_achi, t_int[i], 0.027)
    survival[3*i] = exp(-alpha_T*2-beta_T*2*2)*HT_killing
    
    
    temp_achi=41.
    HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*60*60)
    alpha_T, beta_T = LQ_T(temp_achi, t_int[i], 0.027)
    survival[3*i+1] = exp(-alpha_T*2-beta_T*2*2)*HT_killing
    
    
    temp_achi=43.
    HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*60*60)
    alpha_T, beta_T = LQ_T(temp_achi, t_int[i], 0.027)
    survival[3*i+2] = exp(-alpha_T*2-beta_T*2*2)*HT_killing
    
    
    # Simulate survival at diverse temperatures with mu = 0.5 h^-1
    temp_achi=39.
    HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*60*60)
    alpha_T, beta_T = LQ_T(temp_achi, t_int[i], 0.5)
    survival_high[3*i] = exp(-alpha_T*2-beta_T*2*2)*HT_killing
    
    
    temp_achi=41.
    HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*60*60)
    alpha_T, beta_T = LQ_T(temp_achi, t_int[i], 0.5)
    survival_high[3*i+1] = exp(-alpha_T*2-beta_T*2*2)*HT_killing
    
    
    temp_achi=43.
    HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*60*60)
    alpha_T, beta_T = LQ_T(temp_achi, t_int[i], 0.5)
    survival_high[3*i+2] = exp(-alpha_T*2-beta_T*2*2)*HT_killing
    
    
    # Record time intervals and temperatures
    interval[3*i] = t_int[i]
    interval[3*i+1] = t_int[i]
    interval[3*i+2] = t_int[i]
    T[3*i] = 39
    T[3*i+1] = 41
    T[3*i+2] = 43
    

# Save results
np.savetxt('results/t_int.dat',(np.c_[interval, T, survival, survival_high]), 
           header="t_int T survival survival_high", comments="")
