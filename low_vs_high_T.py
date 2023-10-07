###################################################################################################################
#                                                                                                                 #
#    Four different schedules were used to study the impact of low/high temperature HT sessions:                  #
#                                                                                                                 #
#    1.  A baseline schedule with all the HT sessions at a medium temperature (41 °C).                            #
#    2.  A schedule with one session at low temperature (39 °C) and the rest at a medium temperature (41 °C).     #
#    3.  A schedule with one session at a high temperature (43 °C) and the rest at a medium temperature (41 °C).  #
#    4.  A schedule with a session at high temperature (43 °C) and the rest at a low temperature (39 °C).         #
#                                                                                                                 #
###################################################################################################################


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
N_fract = 30                         # Number of RT sessions
volume=180                           # Tumor volume in cm^3
N=volume/(2e-9)                      # Number of initial cells
SF2 = exp(-alpha_37*2-beta_37*2*2)   #Survival after to 2 Gy (without HT)


# Definition of vectors to save the results
N_HT = np.linspace(3,6,4)                # Number of HT sessions
TCP_10_1 = np.zeros(len(N_HT))           # TCP vector for case 1 (without direct HT cell killing and t_int = 10 min)
TCP_kill_10_1 = np.zeros(len(N_HT))      # TCP vector for case 1 (with direct HT cell killing and t_int = 10 min)
TCP_10_2 = np.zeros(len(N_HT))           # TCP vector for case 2 (without direct HT cell killing and t_int = 10 min)
TCP_kill_10_2 = np.zeros(len(N_HT))      # TCP vector for case 2 (with direct HT cell killing and t_int = 10 min)
TCP_10_3 = np.zeros(len(N_HT))           # TCP vector for case 3 (without direct HT cell killing and t_int = 10 min)
TCP_kill_10_3 = np.zeros(len(N_HT))      # TCP vector for case 3 (with direct HT cell killing and t_int = 10 min)
TCP_10_4 = np.zeros(len(N_HT))           # TCP vector for case 4 (without direct HT cell killing and t_int = 10 min)
TCP_kill_10_4 = np.zeros(len(N_HT))      # TCP vector for case 4 (with direct HT cell killing and t_int = 10 min)

TCP_120_1 = np.zeros(len(N_HT))          # TCP vector for case 1 (without direct HT cell killing and t_int = 120 min)
TCP_kill_120_1 = np.zeros(len(N_HT))     # TCP vector for case 1 (with direct HT cell killing and t_int = 120 min)
TCP_120_2 = np.zeros(len(N_HT))          # TCP vector for case 2 (without direct HT cell killing and t_int = 120 min)
TCP_kill_120_2 = np.zeros(len(N_HT))     # TCP vector for case 2 (with direct HT cell killing and t_int = 120 min)
TCP_120_3 = np.zeros(len(N_HT))          # TCP vector for case 3 (without direct HT cell killing and t_int = 120 min)
TCP_kill_120_3 = np.zeros(len(N_HT))     # TCP vector for case 3 (with direct HT cell killing and t_int = 120 min)
TCP_120_4 = np.zeros(len(N_HT))          # TCP vector for case 4 (without direct HT cell killing and t_int = 120 min)
TCP_kill_120_4 = np.zeros(len(N_HT))     # TCP vector for case 4 (with direct HT cell killing and t_int = 120 min)



# Calculation of direct HT cell killing and radiosensitization effect at 39 °C (t_int = 10 min)
t_int = 10         # Time interval in min
T=39.
HT_killing_39=np.exp(-(2.05e10)*(T+273.15)*exp( S/2 - H/(2*(273.15+T)))*3600)
alpha_T, beta_T = LQ_T(T, t_int)
SF_39 = exp(-alpha_T*2-beta_T*2*2)


# Calculation of direct HT cell killing and radiosensitization effect at 41 °C (t_int = 10 min)
T=41.
HT_killing_41=np.exp(-(2.05e10)*(T+273.15)*exp( S/2 - H/(2*(273.15+T)))*3600)
alpha_T, beta_T = LQ_T(T, t_int)
SF_41 = exp(-alpha_T*2-beta_T*2*2)


# Calculation of direct HT cell killing and radiosensitization effect at 43 °C (t_int = 10 min)
T=43.
HT_killing_43=np.exp(-(2.05e10)*(T+273.15)*exp( S/2 - H/(2*(273.15+T)))*3600)
alpha_T, beta_T = LQ_T(T, t_int)
SF_43 = exp(-alpha_T*2-beta_T*2*2)


# Loop to consider diverse number of HT sessions with t_int = 10 min
for i in range(0,len(N_HT)):
        #Case 1 (All sessions at 41 degrees)
        TCP_10_1[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*(SF_41**N_HT[i]))
        TCP_kill_10_1[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*((SF_41*HT_killing_41)**N_HT[i]))
        
        #Case 2 (One session at 39 degrees and the rest at 41 degrees)
        TCP_10_2[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*SF_39*(SF_41**(N_HT[i]-1)))
        TCP_kill_10_2[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*(SF_39*HT_killing_39)*((SF_41*HT_killing_41)**(N_HT[i]-1)))
        
        #Case 3 (One session at 43 degrees and the rest at 41 degrees)
        TCP_10_3[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*SF_43*(SF_41**(N_HT[i]-1)))
        TCP_kill_10_3[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*(SF_43*HT_killing_43)*((SF_41*HT_killing_41)**(N_HT[i]-1)))
        
        #Case 4 (One session at 43 degrees and the rest at 39 degrees)
        TCP_10_4[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*SF_43*(SF_39**(N_HT[i]-1)))
        TCP_kill_10_4[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*(SF_43*HT_killing_43)*((SF_39*HT_killing_39)**(N_HT[i]-1)))
        
        
        
        
# Calculation of direct HT cell killing and radiosensitization effect at 39 °C (t_int = 120 min)
t_int = 120         # Time interval in min
T=39.
HT_killing_39=np.exp(-(2.05e10)*(T+273.15)*exp( S/2 - H/(2*(273.15+T)))*3600)
alpha_T, beta_T = LQ_T(T, t_int)
SF_39 = exp(-alpha_T*2-beta_T*2*2)


# Calculation of direct HT cell killing and radiosensitization effect at 41 °C (t_int = 120 min)
T=41.
HT_killing_41=np.exp(-(2.05e10)*(T+273.15)*exp( S/2 - H/(2*(273.15+T)))*3600)
alpha_T, beta_T = LQ_T(T, t_int)
SF_41 = exp(-alpha_T*2-beta_T*2*2)


# Calculation of direct HT cell killing and radiosensitization effect at 43 °C (t_int = 120 min)
T=43.
HT_killing_43=np.exp(-(2.05e10)*(T+273.15)*exp( S/2 - H/(2*(273.15+T)))*3600)
alpha_T, beta_T = LQ_T(T, t_int)
SF_43 = exp(-alpha_T*2-beta_T*2*2)


# Loop to consider diverse number of HT sessions with t_int = 120 min
for i in range(0,len(N_HT)):
        #Case 1 (All sessions at 41 degrees)
        TCP_120_1[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*(SF_41**N_HT[i]))
        TCP_kill_120_1[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*((SF_41*HT_killing_41)**N_HT[i]))
        
        #Case 2 (One session at 39 degrees and the rest at 41 degrees)
        TCP_120_2[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*SF_39*(SF_41**(N_HT[i]-1)))
        TCP_kill_120_2[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*(SF_39*HT_killing_39)*((SF_41*HT_killing_41)**(N_HT[i]-1)))
        
        #Case 3 (One session at 43 degrees and the rest at 41 degrees)
        TCP_120_3[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*SF_43*(SF_41**(N_HT[i]-1)))
        TCP_kill_120_3[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*(SF_43*HT_killing_43)*((SF_41*HT_killing_41)**(N_HT[i]-1)))
        
        #Case 4 (One session at 43 degrees and the rest at 39 degrees)
        TCP_120_4[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*SF_43*(SF_39**(N_HT[i]-1)))
        TCP_kill_120_4[i] = exp(-N*(SF2**(N_fract-N_HT[i]))*(SF_43*HT_killing_43)*((SF_39*HT_killing_39)**(N_HT[i]-1)))


# Put together the results and save
schedule = [*np.zeros(len(N_HT)), *np.ones(len(N_HT)), *np.ones(len(N_HT))*2, *np.ones(len(N_HT))*3]
N_HT =[*N_HT, *N_HT, *N_HT, *N_HT]

TCP_10 = [*TCP_10_2, *TCP_10_1, *TCP_10_4, *TCP_10_3]
TCP_kill_10 = [*TCP_kill_10_2, *TCP_kill_10_1, *TCP_kill_10_4, *TCP_kill_10_3]

TCP_120 = [*TCP_120_2, *TCP_120_1, *TCP_120_4, *TCP_120_3]
TCP_kill_120 = [*TCP_kill_120_2, *TCP_kill_120_1, *TCP_kill_120_4, *TCP_kill_120_3]

np.savetxt('results/histogram.dat',(np.c_[N_HT, schedule, TCP_10, TCP_kill_10, TCP_120, TCP_kill_120]), 
           header="N_HT schedule TCP_10 TCP_kill_10 TCP_120 TCP_kill_120", comments="")
