import numpy as np
from math import exp, log
import random


# Definition of model parameters
alpha_37 = 0.386
alpha_41 = 1.73*alpha_37

beta_37 = alpha_37/17.9
beta_41 = 0.41*beta_37

mu = 0.027/60
mu_high = 0.5/60

S = 392.08
H = 147907.8


# Definition of the extended LQ model
def LQ_T(T, t_int, alpha_37, alpha_41, beta_37, beta_41, mu):
    alpha=alpha_37*exp((T-37)*log(alpha_41/alpha_37)*exp(-mu*t_int)/(41-37))
    beta=beta_37*exp((T-37)*log(beta_41/beta_37)*exp(-mu*t_int)/(41-37))
    
    return alpha, beta


# Definition of treatment conditions
N_fract = 30     # Number of RT sessions
D = 2            # Dose per fraction
t_int_min = 10   # Minimum time interval 
t_int_max = 240  # Maximum time interval
T_min=39         # Minimum temperature
T_max=43         # Minimum temperature
N_patients=10000 # Number of studied patients

volume=180       # Tumor volume in cm^3

par_orig = [alpha_37, alpha_41, beta_37, beta_41, 0.027/60, S]                   # Vector of original parameters (mu = 0.027 h^-1) 
diff = [-.25, -.2, -.15, -.1, -.05, -.03, -.01, .01, .03, .05, .1, .15, .2, .25]    # Vector of parameter modificatins


# Definition of vectors to save the results
TCP_mean_var = np.zeros((len(par_orig),len(diff)))     # Vector with mean TCP variation after parameter modification (mu = 0.027 h^-1)
TCP_mean_var_high = np.zeros((len(par_orig),len(diff)))     # Vector with mean TCP variation after parameter modification (mu = 0.027 h^-1)


#Loop over multiple parameters
for j in range(0, len(par_orig)):
    # Loop over the diverse parameters modifications
    for k in range(0,len(diff)):
        # Redefinition of the corresponding model parameter
        par_mod = [alpha_37, alpha_41, beta_37, beta_41, 0.027/60, S] 
        par_mod[j] = par_orig[j]*(1+diff[k])


        alpha_37_diff = par_mod[0]
        alpha_41_diff = par_mod[1]
        beta_37_diff = par_mod[2]
        beta_41_diff = par_mod[3]

        if j==4:
            mu_diff = 0.027*(1+diff[k])/60
            mu_high_diff = 0.5*(1+diff[k])/60
        else:
            mu_diff = 0.027/60
            mu_high_diff = 0.5/60

        S_diff = par_mod[5]

        
        random.seed(10)    # Random generator seed to obtaion always the same results
        

        # Definition of vectors to save the results
        treatments = np.zeros(4*N_patients)             # Number of HT sessions vector
        TCP = np.zeros(4*N_patients)                    # Original TCP vector (mu = 0.027 h^-1) 
        TCP_diff = np.zeros(4*N_patients)               # Modified TCP vector (mu = 0.027 h^-1) 
        TCP_high = np.zeros(4*N_patients)               # Original TCP vector (mu = 0.5 h^-1) 
        TCP_high_diff = np.zeros(4*N_patients)          # Modified TCP vector (mu = 0.5 h^-1)
        TCP_new_error = np.zeros(4*N_patients)          # Mean relative TCP change vector (mu = 0.027 h^-1) 
        TCP_high_new_error = np.zeros(4*N_patients)     # Mean relative TCP change vector (mu = 0.5 h^-1)


        for n in range(0, N_patients):
            # Loop to simulate diverse number of HT sessions
            for n_ht in range(3,7):
                # Definition of initial number of cells
                N=volume/(2e-9)
                N_diff=volume/(2e-9)
                N_high=volume/(2e-9)
                N_high_diff=volume/(2e-9)

                for i in range(0,N_fract):
                    if i<n_ht:
                        # Random selection of time interval
                        t_int = random.uniform(t_int_min, t_int_max)      

                        # Random selection of temperature
                        temp_achi = random.uniform(T_min, T_max)

                        # Calculation of survival after HT+RT fraction (original parameters)
                        HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*3600)
                        alpha_T, beta_T = LQ_T(temp_achi, t_int, alpha_37, alpha_41, beta_37, beta_41, mu)
                        N=N*exp(-alpha_T*D-beta_T*D*D)*HT_killing

                        alpha_T, beta_T = LQ_T(temp_achi, t_int, alpha_37, alpha_41, beta_37, beta_41, mu_high)
                        N_high=N_high*exp(-alpha_T*D-beta_T*D*D)*HT_killing


                        # Calculation of survival after HT+RT fraction (modified parameters)
                        HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S_diff/2 - H/(2*(273.15+temp_achi)))*3600)
                        alpha_T, beta_T = LQ_T(temp_achi, t_int, alpha_37_diff, alpha_41_diff, beta_37_diff, beta_41_diff, mu_diff)
                        N_diff=N_diff*exp(-alpha_T*D-beta_T*D*D)*HT_killing

                        alpha_T, beta_T = LQ_T(temp_achi, t_int, alpha_37_diff, alpha_41_diff, beta_37_diff, beta_41_diff, mu_high_diff)
                        N_high_diff=N_high_diff*exp(-alpha_T*D-beta_T*D*D)*HT_killing

                    else:
                        # Calculation of survival after RT fraction
                        N=N*exp(-alpha_37*D-beta_37*D*D)
                        N_diff=N_diff*exp(-alpha_37_diff*D-beta_37_diff*D*D)
                        N_high=N_high*exp(-alpha_37*D-beta_37*D*D)
                        N_high_diff=N_high_diff*exp(-alpha_37_diff*D-beta_37_diff*D*D)

                # Calculation of TCPs after the HT+RT treatment (mu = 0.027 h^-1)
                TCP[n+(n_ht-3)*N_patients]=exp(-N) 
                TCP_diff[n+(n_ht-3)*N_patients]=exp(-N_diff) 

                # Calculation of TCPs after the HT+RT treatment (mu = 0.5 h^-1)
                TCP_high[n+(n_ht-3)*N_patients]=exp(-N_high) 
                TCP_high_diff[n+(n_ht-3)*N_patients]=exp(-N_high_diff)


               # Calculation of mean TCP relative change
                TCP_new_error[n+(n_ht-3)*N_patients]=abs(TCP_diff[n+(n_ht-3)*N_patients]-TCP[n+(n_ht-3)*N_patients])/TCP[n+(n_ht-3)*N_patients]
                TCP_high_new_error[n+(n_ht-3)*N_patients]=abs(TCP_high_diff[n+(n_ht-3)*N_patients]-TCP_high[n+(n_ht-3)*N_patients])/TCP_high[n+(n_ht-3)*N_patients]

                
        # Calculation of mean TCP relative change (%)
        TCP_mean_var[j,k] = np.mean(TCP_new_error)*100
        TCP_mean_var_high[j,k] = np.mean(TCP_high_new_error)*100
            
        

# Redifinition of the vectors to properly save the results
diff = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
diff = np.append(diff,diff)
mu =  np.append(np.ones(len(TCP_mean_var[0,:]))*0.027,np.ones(len(TCP_mean_var_high[0,:]))*0.5)
TCP_0 = np.append(TCP_mean_var[0,:],TCP_mean_var_high[0,:])
TCP_1 = np.append(TCP_mean_var[1,:],TCP_mean_var_high[1,:])
TCP_2 = np.append(TCP_mean_var[2,:],TCP_mean_var_high[2,:])
TCP_3 = np.append(TCP_mean_var[3,:],TCP_mean_var_high[3,:])
TCP_4 = np.append(TCP_mean_var[4,:],TCP_mean_var_high[4,:])
TCP_5 = np.append(TCP_mean_var[5,:],TCP_mean_var_high[5,:])

np.savetxt('results/sensitivity_analysis.dat',(np.c_[diff, mu, TCP_0, TCP_1, TCP_2, TCP_3, TCP_4, TCP_5]), 
           header="diff mu TCP_0 TCP_1 TCP_2 TCP_3 TCP_4 TCP_5", comments="")
