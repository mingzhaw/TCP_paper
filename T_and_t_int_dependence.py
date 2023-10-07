import numpy as np
from math import exp, log
import random
from scipy import stats

# Definition of model parameters
alpha_37 = 0.386
alpha_41 = 1.73*alpha_37

beta_37 = alpha_37/17.9
beta_41 = 0.41*beta_37

S = 392.08
H = 147907.8

# Definition of the extended LQ model
def LQ_T(T, t_int, mu):
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


# Definition of vectors to save the results
Temp = np.zeros(4*N_patients)                   # Mean temperature vector
interval = np.zeros(4*N_patients)               # Mean time interval vector
treatments = np.zeros(4*N_patients)             # Number of HT sessions vector
TCP = np.zeros(4*N_patients)                    # TCP vector (without direct HT cell killing and mu = 0.027 h^-1) 
TCP_kill = np.zeros(4*N_patients)               # TCP vector (with direct HT cell killing and mu = 0.027 h^-1) 
TCP_high = np.zeros(4*N_patients)               # TCP vector (without direct HT cell killing and mu = 0.5 h^-1) 
TCP_kill_high = np.zeros(4*N_patients)          # TCP vector (with direct HT cell killing and mu = 0.5 h^-1)


# Calculation of the baseline TCP (only RT)
N=volume/(2e-9)
N=N*exp(-N_fract*alpha_37*2-N_fract*beta_37*2*2)
baseline = np.exp(-N)

print("Baseline TCP (only RT): ",baseline)
print()


random.seed(10)   # Random generator seed to obtaion always the same results

for n in range(0, N_patients):
    # Loop to simulate diverse number of HT sessions
    for n_ht in range(3,7):
        # Definition of initial number of cells
        N=volume/(2e-9)
        N_kill=volume/(2e-9)
        N_high=volume/(2e-9)
        N_kill_high=volume/(2e-9)
        
        # Record number of HT sessions
        treatments[n+(n_ht-3)*N_patients] = n_ht
    
        for i in range(0,N_fract):
            if i<n_ht:
                # Random selection of time interval
                t_int = random.uniform(t_int_min, t_int_max)
                
                # Record mean time interval
                interval[n+(n_ht-3)*N_patients] += t_int/n_ht       
                
                # Random selection of temperature
                temp_achi = random.uniform(T_min, T_max)
                
                # Record mean temperature
                Temp[n+(n_ht-3)*N_patients] += temp_achi/n_ht
                
                # Calculation of survival after HT+RT fraction (mu = 0.027 h^-1)
                HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*3600)
                alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.027/60)
                
                N=N*exp(-alpha_T*D-beta_T*D*D)
                N_kill=N_kill*exp(-alpha_T*D-beta_T*D*D)*HT_killing
                
                # Calculation of survival after HT+RT fraction (mu = 0.5 h^-1)
                HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*3600)
                alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.5/60)

                N_high=N_high*exp(-alpha_T*D-beta_T*D*D)
                N_kill_high=N_kill_high*exp(-alpha_T*D-beta_T*D*D)*HT_killing

            else:
                # Calculation of survival after RT fraction
                N=N*exp(-alpha_37*D-beta_37*D*D)
                N_kill=N_kill*exp(-alpha_37*D-beta_37*D*D)
                N_high=N_high*exp(-alpha_37*D-beta_37*D*D)
                N_kill_high=N_kill_high*exp(-alpha_37*D-beta_37*D*D)

        # Calculation of TCPs after the HT+RT treatment (mu = 0.027 h^-1)
        TCP[n+(n_ht-3)*N_patients]=exp(-N) 
        TCP_kill[n+(n_ht-3)*N_patients]=exp(-N_kill) 

        # Calculation of TCPs after the HT+RT treatment (mu = 0.5 h^-1)
        TCP_high[n+(n_ht-3)*N_patients]=exp(-N_high) 
        TCP_kill_high[n+(n_ht-3)*N_patients]=exp(-N_kill_high)

   
    
# Save results
np.savetxt('results/TCP_random.dat',(np.c_[treatments, interval, Temp, TCP, TCP_kill, TCP_high, TCP_kill_high]), 
           header="N_HT interval Temp TCP TCP_kill TCP_high TCP_kill_high", comments="")

# Calculation TCP improvement
print("TCP (Low mu)")
print("Treatment improvement (no kill): ",np.min(TCP)-baseline," - ",np.max(TCP)-baseline)
print("Treatment improvement (kill): ",np.min(TCP_kill)-baseline," - ",np.max(TCP_kill)-baseline)
print()
print("TCP (High mu)")
print("Treatment improvement (no kill): ",np.min(TCP_high)-baseline," - ",np.max(TCP_high)-baseline)
print("Treatment improvement (kill): ",np.min(TCP_kill_high)-baseline," - ",np.max(TCP_kill_high)-baseline)
print()


# Calculation of correlations for diverse number of HT sessions
spearman_T = np.zeros(4)
spearman_tint = np.zeros(4)
spearman_p_T = np.zeros(4)
spearman_p_tint = np.zeros(4)
kendall_T = np.zeros(4)
kendall_tint = np.zeros(4)
kendall_p_T = np.zeros(4)
kendall_p_tint = np.zeros(4)

spearman_T_high = np.zeros(4)
spearman_tint_high = np.zeros(4)
spearman_p_T_high = np.zeros(4)
spearman_p_tint_high = np.zeros(4)
kendall_T_high = np.zeros(4)
kendall_tint_high = np.zeros(4)
kendall_p_T_high = np.zeros(4)
kendall_p_tint_high = np.zeros(4)

spearman_T_kill = np.zeros(4)
spearman_tint_kill = np.zeros(4)
spearman_p_T_kill = np.zeros(4)
spearman_p_tint_kill = np.zeros(4)
kendall_T_kill = np.zeros(4)
kendall_tint_kill = np.zeros(4)
kendall_p_T_kill = np.zeros(4)
kendall_p_tint_kill = np.zeros(4)

spearman_T_high_kill = np.zeros(4)
spearman_tint_high_kill = np.zeros(4)
spearman_p_T_high_kill = np.zeros(4)
spearman_p_tint_high_kill = np.zeros(4)
kendall_T_high_kill = np.zeros(4)
kendall_tint_high_kill = np.zeros(4)
kendall_p_T_high_kill = np.zeros(4)
kendall_p_tint_high_kill = np.zeros(4)

for i in range(0,4):
    #Temperature
    res = stats.spearmanr(Temp[treatments==(i+3)], TCP[treatments==(i+3)])
    spearman_T[i] = res.correlation
    spearman_p_T[i] = res.pvalue
    kendall_T[i], kendall_p_T[i] = stats.kendalltau(Temp[treatments==(i+3)], TCP[treatments==(i+3)])
    
    res = stats.spearmanr(Temp[treatments==(i+3)], TCP_high[treatments==(i+3)])
    spearman_T_high[i] = res.correlation
    spearman_p_T_high[i] = res.pvalue
    kendall_T_high[i], kendall_p_T_high[i] = stats.kendalltau(Temp[treatments==(i+3)], TCP_high[treatments==(i+3)])
    
    res = stats.spearmanr(Temp[treatments==(i+3)], TCP_kill[treatments==(i+3)])
    spearman_T_kill[i] = res.correlation
    spearman_p_T_kill[i] = res.pvalue
    kendall_T_kill[i], kendall_p_T_kill[i] = stats.kendalltau(Temp[treatments==(i+3)], TCP_kill[treatments==(i+3)])
    
    res = stats.spearmanr(Temp[treatments==(i+3)], TCP_kill_high[treatments==(i+3)])
    spearman_T_high_kill[i] = res.correlation
    spearman_p_T_high_kill[i] = res.pvalue
    kendall_T_high_kill[i], kendall_p_T_high_kill[i] = stats.kendalltau(Temp[treatments==(i+3)], TCP_kill_high[treatments==(i+3)])
    
    #Time interval
    res = stats.spearmanr(interval[treatments==(i+3)], TCP[treatments==(i+3)])
    spearman_tint[i] = res.correlation
    spearman_p_tint[i] = res.pvalue
    kendall_tint[i], kendall_p_tint[i] = stats.kendalltau(interval[treatments==(i+3)], TCP[treatments==(i+3)])
    
    res = stats.spearmanr(interval[treatments==(i+3)], TCP_high[treatments==(i+3)])
    spearman_tint_high[i] = res.correlation
    spearman_p_tint_high[i] = res.pvalue
    kendall_tint_high[i], kendall_p_tint_high[i] = stats.kendalltau(interval[treatments==(i+3)], TCP_high[treatments==(i+3)])
    
    res = stats.spearmanr(interval[treatments==(i+3)], TCP_kill[treatments==(i+3)])
    spearman_tint_kill[i] = res.correlation
    spearman_p_tint_kill[i] = res.pvalue
    kendall_tint_kill[i], kendall_p_tint_kill[i] = stats.kendalltau(interval[treatments==(i+3)], TCP_kill[treatments==(i+3)])
    
    res = stats.spearmanr(interval[treatments==(i+3)], TCP_kill_high[treatments==(i+3)])
    spearman_tint_high_kill[i] = res.correlation
    spearman_p_tint_high_kill[i] = res.pvalue
    kendall_tint_high_kill[i], kendall_p_tint_high_kill[i] = stats.kendalltau(interval[treatments==(i+3)], TCP_kill_high[treatments==(i+3)])
    


    
print()
print("Correlation (Low mu, no killing)")    
print("Mean T Spearman's rho: ",  np.min(spearman_T),"-",np.max(spearman_T))
print("Mean T Kendall's tau: ",  np.min(kendall_T),"-",np.max(kendall_T))
print()
print("Mean t_int Spearman's rho: ",  np.min(spearman_tint),"-",np.max(spearman_tint))
print("Mean t_int Kendall's tau: ",  np.min(kendall_tint),"-",np.max(kendall_tint))
print()

print()
print("Correlation (Low mu, killing)")    
print("Mean T Spearman's rho: ",  np.min(spearman_T_kill),"-",np.max(spearman_T_kill))
print("Mean T Kendall's tau: ",  np.min(kendall_T_kill),"-",np.max(kendall_T_kill))
print()
print("Mean t_int Spearman's rho: ",  np.min(spearman_tint_kill),"-",np.max(spearman_tint_kill))
print("Mean t_int Kendall's tau: ",  np.min(kendall_tint_kill),"-",np.max(kendall_tint_kill))
print()

print()
print("Correlation (High mu, no killing)")    
print("Mean T Spearman's rho: ",  np.min(spearman_T_high),"-",np.max(spearman_T_high))
print("Mean T Kendall's tau: ",  np.min(kendall_T_high),"-",np.max(kendall_T_high))
print()
print("Mean t_int Spearman's rho: ",  np.min(spearman_tint_high),"-",np.max(spearman_tint_high))
print("Mean t_int Kendall's tau: ",  np.min(kendall_tint_high),"-",np.max(kendall_tint_high))
print()

print()
print("Correlation (High mu, killing)")    
print("Mean T Spearman's rho: ", np.min(spearman_T_high_kill),"-",np.max(spearman_T_high_kill))
print("Mean T Kendall's tau: ",  np.min(kendall_T_high_kill),"-",np.max(kendall_T_high_kill))
print()
print("Mean t_int Spearman's rho: ",  np.min(spearman_tint_high_kill),"-",np.max(spearman_tint_high_kill))
print("Mean t_int Kendall's tau: ",  np.min(kendall_tint_high_kill),"-",np.max(kendall_tint_high_kill))
print()


print("p-values: ",spearman_p_T," ",kendall_p_T," ",spearman_p_T_kill," ",kendall_p_T_kill)
print("p-values: ",spearman_p_T_high," ",kendall_p_T_high," ",spearman_p_T_high_kill," ",kendall_p_T_high_kill)

print("p-values: ",spearman_p_tint," ",kendall_p_tint," ",spearman_p_tint_kill," ",kendall_p_tint_kill)
print("p-values: ",spearman_p_tint_high," ",kendall_p_tint_high," ",spearman_p_tint_high_kill," ",kendall_p_tint_high_kill)
print()
