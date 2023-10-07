import numpy as np
from math import exp, log
import random

# Definition of model parameters
alpha_37 = 0.386
alpha_41 = 1.73*alpha_37

beta_37 = alpha_37/17.9
beta_41 = 0.41*beta_37

mu = 0.5/60 #in min-1
S = 392.08
H = 147907.8


# Definition of the extended LQ model
def LQ_T(T, t_int):
    alpha=alpha_37*exp((T-37)*log(alpha_41/alpha_37)*exp(-mu*t_int)/(41-37))
    beta=beta_37*exp((T-37)*log(beta_41/beta_37)*exp(-mu*t_int)/(41-37))
    
    return alpha, beta


# Definition of treatment conditions
N_fract=61                           # Number of RT fractions
SF2 = exp(-alpha_37*2-beta_37*2*2)   # Survival fraction after 2 Gy (without HT)
volume=180                           # Tumor volume in cm^3
N=volume/(2e-9)                      # Number of initial cells
N_patients = 10000                   # Number of simulated patients


# Definition of vectors to save the results
TCP = np.zeros(N_fract)                                  # Mean TCP vector (no heterogeneity) 
TCP_het_alpha = np.zeros(N_fract)                        # Mean TCP vector (radiosensitivity heterogeneity)
TCP_het_vol = np.zeros(N_fract)                          # Mean TCP vector (volume heterogeneity)
TCP_het_both = np.zeros(N_fract)                         # Mean TCP vector (radiosensitivity and volume heterogeneity)
D = np.zeros(N_fract)                                    # Total delivered dose vector
alpha = np.random.normal(0.386, 0.386*0.5,N_patients)    # Diverse alpha_37 values vector
V = np.random.normal(201.6, 119.7,N_patients)            # Diverse tumor volume vector


# Calculation of TCP (without heterogeneity)
for i in range(0,N_fract):
    TCP[i] = exp(-N*(SF2**i))
    D[i] = 2*i


# Calculation of TCP (with heterogeneity)
for j in range(0,N_patients):
    # Check that the radiosensitivity is within the defined limits
    while alpha[j]<0.2 or alpha[j]>0.74:
        alpha[j] = np.random.normal(0.386, 0.386*0.5)
    
    # Check that the volume is not negative
    while V[j]<0:
        V[j] = np.random.normal(201.6, 119.7)
        
    # Redefinition the initial number of cells
    N_vol=V[j]/(2e-9)
    
    # Redefinition the Survival fraction after 2 Gy (without HT)
    SF2_alpha = np.exp(-2*alpha[j]-4*alpha[j]/17.9)
    
    for i in range(0,N_fract):
        TCP_het_alpha[i] += exp(-N*(SF2_alpha**i))/N_patients
        TCP_het_vol[i] += exp(-N_vol*(SF2**i))/N_patients
        TCP_het_both[i] += exp(-N_vol*(SF2_alpha**i))/N_patients

# Save results
np.savetxt('results/dose_control_het.dat',(np.c_[D, TCP, TCP_het_alpha, TCP_het_vol, TCP_het_both]), 
           header="D TCP TCP_het_alpha TCP_het_vol TCP_het_both", comments="")
