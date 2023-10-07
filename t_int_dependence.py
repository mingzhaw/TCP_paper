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
N_HT = 4         # Number of HT sessions
t_int_min = 10   # Minimum time interval 
t_int_max = 240  # Maximum time interval
temp_achi=39     #Temperature case 1
temp_achi_1=41   #Temperature case 2
temp_achi_2=43   #Temperature case 3
N_patients=10000 # Number of studied patients

volume=180       # Tumor volume in cm^3

# Definition of vectors to save the results
T = np.zeros(3*N_patients)                      # Temperature vector
interval = np.zeros(3*N_patients)               # Mean time interval vector
interval_min = np.ones(3*N_patients)*t_int_max  # Minimum time interval vector
interval_max = np.zeros(3*N_patients)           # Maximum time interval vector
TCP = np.zeros(3*N_patients)                    # TCP vector (without direct HT cell killing and mu = 0.027 h^-1) 
TCP_kill = np.zeros(3*N_patients)               # TCP vector (with direct HT cell killing and mu = 0.027 h^-1) 
TCP_high = np.zeros(3*N_patients)               # TCP vector (without direct HT cell killing and mu = 0.5 h^-1) 
TCP_high_kill = np.zeros(3*N_patients)          # TCP vector (with direct HT cell killing and mu = 0.5 h^-1) 


# Calculation of the baseline TCP (only RT)
N=volume/(2e-9)
N=N*exp(-N_fract*alpha_37*2-N_fract*beta_37*2*2)
baseline = np.exp(-N)

print("Baseline TCP (only RT): ",baseline)
print()


random.seed(10)   # Random generator seed to obtaion always the same results


# Loop to calculate the TCP values
for n in range(0, N_patients):
    #Definition of initial number of cells
    N=volume/(2e-9)
    N1=volume/(2e-9)
    N2=volume/(2e-9)
    
    N_high=volume/(2e-9)
    N1_high=volume/(2e-9)
    N2_high=volume/(2e-9)
    
    N_kill=volume/(2e-9)
    N1_kill=volume/(2e-9)
    N2_kill=volume/(2e-9)
    
    N_high_kill=volume/(2e-9)
    N1_high_kill=volume/(2e-9)
    N2_high_kill=volume/(2e-9)
        
    for i in range(0,N_fract):
        if i<N_HT:
            # Random selection of time interval
            t_int = random.uniform(t_int_min, t_int_max)
            
            # Record mean time interval
            interval[3*n] += t_int/4
            interval[3*n+1] += t_int/4
            interval[3*n+2] += t_int/4
            
            # Record maximum time interval
            if t_int>interval_max[3*n]:
                interval_max[3*n]=t_int
                interval_max[3*n+1]=t_int
                interval_max[3*n+2]=t_int
            
            # Record minimum time interval
            if t_int<interval_min[3*n]:
                interval_min[3*n]=t_int
                interval_min[3*n+1]=t_int
                interval_min[3*n+2]=t_int
            
            # Calculation of survival after HT+RT fraction case 1 (mu = 0.027 h^-1)
            HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*3600)
            alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.027/60)
            N=N*exp(-alpha_T*D-beta_T*D*D)
            N_kill=N_kill*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
            # Calculation of survival after HT+RT fraction case 2 (mu = 0.027 h^-1)
            HT_killing=np.exp(-(2.05e10)*(temp_achi_1+273.15)*exp( S/2 - H/(2*(273.15+temp_achi_1)))*3600)
            alpha_T, beta_T = LQ_T(temp_achi_1, t_int, 0.027/60)
            N1=N1*exp(-alpha_T*D-beta_T*D*D)
            N1_kill=N1_kill*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
            # Calculation of survival after HT+RT fraction case 3 (mu = 0.027 h^-1)
            HT_killing=np.exp(-(2.05e10)*(temp_achi_2+273.15)*exp( S/2 - H/(2*(273.15+temp_achi_2)))*3600)
            alpha_T, beta_T = LQ_T(temp_achi_2, t_int, 0.027/60)
            N2=N2*exp(-alpha_T*D-beta_T*D*D)
            N2_kill=N2_kill*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
            # Calculation of survival after HT+RT fraction case 1 (mu = 0.5 h^-1)
            HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*3600)
            alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.5/60)
            N_high=N_high*exp(-alpha_T*D-beta_T*D*D)
            N_high_kill=N_high_kill*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
            # Calculation of survival after HT+RT fraction case 2 (mu = 0.5 h^-1)
            HT_killing=np.exp(-(2.05e10)*(temp_achi_1+273.15)*exp( S/2 - H/(2*(273.15+temp_achi_1)))*3600)
            alpha_T, beta_T = LQ_T(temp_achi_1, t_int, 0.5/60)
            N1_high=N1_high*exp(-alpha_T*D-beta_T*D*D)
            N1_high_kill=N1_high_kill*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
            # Calculation of survival after HT+RT fraction case 3 (mu = 0.5 h^-1)
            HT_killing=np.exp(-(2.05e10)*(temp_achi_2+273.15)*exp( S/2 - H/(2*(273.15+temp_achi_2)))*3600)
            alpha_T, beta_T = LQ_T(temp_achi_2, t_int, 0.5/60)
            N2_high=N2_high*exp(-alpha_T*D-beta_T*D*D)
            N2_high_kill=N2_high_kill*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
        else:
            # Calculation of survival after RT fraction 
            N=N*exp(-alpha_37*D-beta_37*D*D)
            N1=N1*exp(-alpha_37*D-beta_37*D*D)
            N2=N2*exp(-alpha_37*D-beta_37*D*D)
            
            N_kill=N_kill*exp(-alpha_37*D-beta_37*D*D)
            N1_kill=N1_kill*exp(-alpha_37*D-beta_37*D*D)
            N2_kill=N2_kill*exp(-alpha_37*D-beta_37*D*D)
            
            N_high=N_high*exp(-alpha_37*D-beta_37*D*D)
            N1_high=N1_high*exp(-alpha_37*D-beta_37*D*D)
            N2_high=N2_high*exp(-alpha_37*D-beta_37*D*D)
            
            N_high_kill=N_high_kill*exp(-alpha_37*D-beta_37*D*D)
            N1_high_kill=N1_high_kill*exp(-alpha_37*D-beta_37*D*D)
            N2_high_kill=N2_high_kill*exp(-alpha_37*D-beta_37*D*D)
            
        
        
    # Record temperature    
    T[3*n] = 39
    T[3*n+1] = 41
    T[3*n+2] = 43
    
    # Calculation of TCPs after the HT+RT treatment (mu = 0.027 h^-1)
    TCP[3*n]=exp(-N)  
    TCP[3*n+1]=exp(-N1) 
    TCP[3*n+2]=exp(-N2)
    
    TCP_kill[3*n]=exp(-N_kill)  
    TCP_kill[3*n+1]=exp(-N1_kill) 
    TCP_kill[3*n+2]=exp(-N2_kill)
    
    # Calculation of TCPs after the HT+RT treatment (mu = 0.5 h^-1)
    TCP_high[3*n]=exp(-N_high) 
    TCP_high[3*n+1]=exp(-N1_high) 
    TCP_high[3*n+2]=exp(-N2_high)
    
    TCP_high_kill[3*n]=exp(-N_high_kill) 
    TCP_high_kill[3*n+1]=exp(-N1_high_kill) 
    TCP_high_kill[3*n+2]=exp(-N2_high_kill)
    

print("Low mu")
print("TCP (39 °C): ",np.min(TCP[T==39])," - ",np.max(TCP[T==39]),"  Difference=",np.max(TCP[T==39])-np.min(TCP[T==39]))
print("TCP (41 °C): ",np.min(TCP[T==41])," - ",np.max(TCP[T==41]),"  Difference=",np.max(TCP[T==41])-np.min(TCP[T==41]))
print("TCP (43 °C): ",np.min(TCP[T==43])," - ",np.max(TCP[T==43]),"  Difference=",np.max(TCP[T==43])-np.min(TCP[T==43]))
print()
print()
print("High mu")
print("TCP (39 °C): ",np.min(TCP_high[T==39])," - ",np.max(TCP_high[T==39]),"  Difference=",np.max(TCP_high[T==39])-np.min(TCP_high[T==39]))
print("TCP (41 °C): ",np.min(TCP_high[T==41])," - ",np.max(TCP_high[T==41]),"  Difference=",np.max(TCP_high[T==41])-np.min(TCP_high[T==41]))
print("TCP (43 °C): ",np.min(TCP_high[T==43])," - ",np.max(TCP_high[T==43]),"  Difference=",np.max(TCP_high[T==43])-np.min(TCP_high[T==43]))
print()

# Save results
np.savetxt('results/TCP_t_int.dat',(np.c_[interval, T, TCP, TCP_high, TCP_kill, TCP_high_kill]), 
           header="interval T TCP TCP_high TCP_kill TCP_high_kill", comments="")


# Calculation of TCP ranges
SF2 = exp(-alpha_37*2-beta_37*2*2)
N=volume/(2e-9)

print("mu = 0.027 h^-1")
temp_achi=39.
t_int=10
HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*60*60)
alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.027/60)
SF = exp(-alpha_T*2-beta_T*2*2)
TCP_max = exp(-N*(SF2**(N_fract-4))*(SF**4))
TCP_max_kill = exp(-N*(SF2**(N_fract-4))*((SF*HT_killing)**4))

t_int=240
alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.027/60)
SF = exp(-alpha_T*2-beta_T*2*2)
TCP_min = exp(-N*(SF2**(N_fract-4))*(SF**4))
TCP_min_kill = exp(-N*(SF2**(N_fract-4))*((SF*HT_killing)**4))

print("TCP (39 °C without direct HT cell killing): ",TCP_min," - ",TCP_max,"  Difference=",TCP_max-TCP_min)
print("TCP (39 °C with direct HT cell killing): ",TCP_min_kill," - ",TCP_max_kill,"  Difference=",TCP_max_kill-TCP_min_kill)

temp_achi=41.
t_int=10
HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*60*60)
alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.027/60)
SF = exp(-alpha_T*2-beta_T*2*2)
TCP_max = exp(-N*(SF2**(N_fract-4))*(SF**4))
TCP_max_kill = exp(-N*(SF2**(N_fract-4))*((SF*HT_killing)**4))

t_int=240
alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.027/60)
SF = exp(-alpha_T*2-beta_T*2*2)
TCP_min = exp(-N*(SF2**(N_fract-4))*(SF**4))
TCP_min_kill = exp(-N*(SF2**(N_fract-4))*((SF*HT_killing)**4))

print("TCP (41 °C without direct HT cell killing): ",TCP_min," - ",TCP_max,"  Difference=",TCP_max-TCP_min)
print("TCP (41 °C with direct HT cell killing): ",TCP_min_kill," - ",TCP_max_kill,"  Difference=",TCP_max_kill-TCP_min_kill)

temp_achi=43.
t_int=10
HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*60*60)
alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.027/60)
SF = exp(-alpha_T*2-beta_T*2*2)
TCP_max = exp(-N*(SF2**(N_fract-4))*(SF**4))
TCP_max_kill = exp(-N*(SF2**(N_fract-4))*((SF*HT_killing)**4))

t_int=240
alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.027/60)
SF = exp(-alpha_T*2-beta_T*2*2)
TCP_min = exp(-N*(SF2**(N_fract-4))*(SF**4))
TCP_min_kill = exp(-N*(SF2**(N_fract-4))*((SF*HT_killing)**4))

print("TCP (43 °C without direct HT cell killing): ",TCP_min," - ",TCP_max,"  Difference=",TCP_max-TCP_min)
print("TCP (43 °C with direct HT cell killing): ",TCP_min_kill," - ",TCP_max_kill,"  Difference=",TCP_max_kill-TCP_min_kill)
print()

print("mu = 0.5 h^-1")
temp_achi=39.
t_int=10
HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*60*60)
alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.5/60)
SF = exp(-alpha_T*2-beta_T*2*2)
TCP_max = exp(-N*(SF2**(N_fract-4))*(SF**4))
TCP_max_kill = exp(-N*(SF2**(N_fract-4))*((SF*HT_killing)**4))

t_int=240
alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.5/60)
SF = exp(-alpha_T*2-beta_T*2*2)
TCP_min = exp(-N*(SF2**(N_fract-4))*(SF**4))
TCP_min_kill = exp(-N*(SF2**(N_fract-4))*((SF*HT_killing)**4))

print("TCP (39 °C without direct HT cell killing): ",TCP_min," - ",TCP_max,"  Difference=",TCP_max-TCP_min)
print("TCP (39 °C with direct HT cell killing): ",TCP_min_kill," - ",TCP_max_kill,"  Difference=",TCP_max_kill-TCP_min_kill)
print("TCP (39 °C only direct HT cell killing): ",exp(-N*(SF2**(N_fract))*((HT_killing)**4)))

temp_achi=41.
t_int=10
HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*60*60)
alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.5/60)
SF = exp(-alpha_T*2-beta_T*2*2)
TCP_max = exp(-N*(SF2**(N_fract-4))*(SF**4))
TCP_max_kill = exp(-N*(SF2**(N_fract-4))*((SF*HT_killing)**4))

t_int=240
alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.5/60)
SF = exp(-alpha_T*2-beta_T*2*2)
TCP_min = exp(-N*(SF2**(N_fract-4))*(SF**4))
TCP_min_kill = exp(-N*(SF2**(N_fract-4))*((SF*HT_killing)**4))

print("TCP (41 °C without direct HT cell killing): ",TCP_min," - ",TCP_max,"  Difference=",TCP_max-TCP_min)
print("TCP (41 °C with direct HT cell killing): ",TCP_min_kill," - ",TCP_max_kill,"  Difference=",TCP_max_kill-TCP_min_kill)
print("TCP (41 °C only direct HT cell killing): ",exp(-N*(SF2**(N_fract))*((HT_killing)**4)))

temp_achi=43.
t_int=10
HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*60*60)
alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.5/60)
SF = exp(-alpha_T*2-beta_T*2*2)
TCP_max = exp(-N*(SF2**(N_fract-4))*(SF**4))
TCP_max_kill = exp(-N*(SF2**(N_fract-4))*((SF*HT_killing)**4))


t_int=240
alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.5/60)
SF = exp(-alpha_T*2-beta_T*2*2)
TCP_min = exp(-N*(SF2**(N_fract-4))*(SF**4))
TCP_min_kill = exp(-N*(SF2**(N_fract-4))*((SF*HT_killing)**4))

print("TCP (43 °C without direct HT cell killing): ",TCP_min," - ",TCP_max,"  Difference=",TCP_max-TCP_min)
print("TCP (43 °C with direct HT cell killing): ",TCP_min_kill," - ",TCP_max_kill,"  Difference=",TCP_max_kill-TCP_min_kill)
print("TCP (43 °C only direct HT cell killing): ",exp(-N*(SF2**(N_fract))*((HT_killing)**4)))
print()
print()


# Calculation of correlarions
print("Correlations (mu = 0.027 h^-1 without direct HT cell killing)")
res = stats.spearmanr(interval_min[T==39], TCP[T==39])
spearman_t_int_0 = res.correlation
spearman_p_t_int_0 = res.pvalue
kendall_t_int_0, kendall_p_t_int_0 = stats.kendalltau(interval_min[T==39], TCP[T==39])

res = stats.spearmanr(interval_min[T==41], TCP[T==41])
spearman_t_int_1 = res.correlation
spearman_p_t_int_1 = res.pvalue
kendall_t_int_1, kendall_p_t_int_1 = stats.kendalltau(interval_min[T==41], TCP[T==41])

res = stats.spearmanr(interval_min[T==43], TCP[T==43])
spearman_t_int_2 = res.correlation
spearman_p_t_int_2 = res.pvalue
kendall_t_int_2, kendall_p_t_int_2 = stats.kendalltau(interval_min[T==43], TCP[T==43])

print("Min t_int Spearman's rho (39 °C, 41 °C, 43 °C): ", spearman_t_int_0,", ",spearman_t_int_1,", ",spearman_t_int_2)
print("Min t_int Kendall's tau (39 °C, 41 °C, 43 °C): ", kendall_t_int_0,", ",kendall_t_int_1,", ",kendall_t_int_2)
print("p-values: ",spearman_p_t_int_0," ",spearman_p_t_int_1," ",spearman_p_t_int_2," ",kendall_p_t_int_0," ",kendall_p_t_int_1," ",kendall_p_t_int_2)
print()

res = stats.spearmanr(interval[T==39], TCP[T==39])
spearman_t_int_0 = res.correlation
spearman_p_t_int_0 = res.pvalue
kendall_t_int_0, kendall_p_t_int_0 = stats.kendalltau(interval[T==39], TCP[T==39])

res = stats.spearmanr(interval[T==41], TCP[T==41])
spearman_t_int_1 = res.correlation
spearman_p_t_int_1 = res.pvalue
kendall_t_int_1, kendall_p_t_int_1 = stats.kendalltau(interval[T==41], TCP[T==41])

res = stats.spearmanr(interval[T==43], TCP[T==43])
spearman_t_int_2 = res.correlation
spearman_p_t_int_2 = res.pvalue
kendall_t_int_2, kendall_p_t_int_2 = stats.kendalltau(interval[T==43], TCP[T==43])

print("Mean t_int Spearman's rho (39 °C, 41 °C, 43 °C): ", spearman_t_int_0,", ",spearman_t_int_1,", ",spearman_t_int_2)
print("Mean t_int Kendall's tau (39 °C, 41 °C, 43 °C): ", kendall_t_int_0,", ",kendall_t_int_1,", ",kendall_t_int_2)
print("p-values: ",spearman_p_t_int_0," ",spearman_p_t_int_1," ",spearman_p_t_int_2," ",kendall_p_t_int_0," ",kendall_p_t_int_1," ",kendall_p_t_int_2)
print()

res = stats.spearmanr(interval_max[T==39], TCP[T==39])
spearman_t_int_0 = res.correlation
spearman_p_t_int_0 = res.pvalue
kendall_t_int_0, kendall_p_t_int_0 = stats.kendalltau(interval_max[T==39], TCP[T==39])

res = stats.spearmanr(interval_max[T==41], TCP[T==41])
spearman_t_int_1 = res.correlation
spearman_p_t_int_1 = res.pvalue
kendall_t_int_1, kendall_p_t_int_1 = stats.kendalltau(interval_max[T==41], TCP[T==41])

res = stats.spearmanr(interval_max[T==43], TCP[T==43])
spearman_t_int_2 = res.correlation
spearman_p_t_int_2 = res.pvalue
kendall_t_int_2, kendall_p_t_int_2 = stats.kendalltau(interval_max[T==43], TCP[T==43])

print("Max t_int Spearman's rho (39 °C, 41 °C, 43 °C): ", spearman_t_int_0,", ",spearman_t_int_1,", ",spearman_t_int_2)
print("Max t_int Kendall's tau (39 °C, 41 °C, 43 °C): ", kendall_t_int_0,", ",kendall_t_int_1,", ",kendall_t_int_2)
print("p-values: ",spearman_p_t_int_0," ",spearman_p_t_int_1," ",spearman_p_t_int_2," ",kendall_p_t_int_0," ",kendall_p_t_int_1," ",kendall_p_t_int_2)
print()
print()


# Calculation of correlarions
print("Correlations (mu = 0.027 h^-1 with direct HT cell killing)")
res = stats.spearmanr(interval_min[T==39], TCP_kill[T==39])
spearman_t_int_0 = res.correlation
spearman_p_t_int_0 = res.pvalue
kendall_t_int_0, kendall_p_t_int_0 = stats.kendalltau(interval_min[T==39], TCP_kill[T==39])

res = stats.spearmanr(interval_min[T==41], TCP_kill[T==41])
spearman_t_int_1 = res.correlation
spearman_p_t_int_1 = res.pvalue
kendall_t_int_1, kendall_p_t_int_1 = stats.kendalltau(interval_min[T==41], TCP_kill[T==41])

res = stats.spearmanr(interval_min[T==43], TCP_kill[T==43])
spearman_t_int_2 = res.correlation
spearman_p_t_int_2 = res.pvalue
kendall_t_int_2, kendall_p_t_int_2 = stats.kendalltau(interval_min[T==43], TCP_kill[T==43])

print("Min t_int Spearman's rho (39 °C, 41 °C, 43 °C): ", spearman_t_int_0,", ",spearman_t_int_1,", ",spearman_t_int_2)
print("Min t_int Kendall's tau (39 °C, 41 °C, 43 °C): ", kendall_t_int_0,", ",kendall_t_int_1,", ",kendall_t_int_2)
print("p-values: ",spearman_p_t_int_0," ",spearman_p_t_int_1," ",spearman_p_t_int_2," ",kendall_p_t_int_0," ",kendall_p_t_int_1," ",kendall_p_t_int_2)
print()

res = stats.spearmanr(interval[T==39], TCP_kill[T==39])
spearman_t_int_0 = res.correlation
spearman_p_t_int_0 = res.pvalue
kendall_t_int_0, kendall_p_t_int_0 = stats.kendalltau(interval[T==39], TCP_kill[T==39])

res = stats.spearmanr(interval[T==41], TCP_kill[T==41])
spearman_t_int_1 = res.correlation
spearman_p_t_int_1 = res.pvalue
kendall_t_int_1, kendall_p_t_int_1 = stats.kendalltau(interval[T==41], TCP_kill[T==41])

res = stats.spearmanr(interval[T==43], TCP_kill[T==43])
spearman_t_int_2 = res.correlation
spearman_p_t_int_2 = res.pvalue
kendall_t_int_2, kendall_p_t_int_2 = stats.kendalltau(interval[T==43], TCP_kill[T==43])

print("Mean t_int Spearman's rho (39 °C, 41 °C, 43 °C): ", spearman_t_int_0,", ",spearman_t_int_1,", ",spearman_t_int_2)
print("Mean t_int Kendall's tau (39 °C, 41 °C, 43 °C): ", kendall_t_int_0,", ",kendall_t_int_1,", ",kendall_t_int_2)
print("p-values: ",spearman_p_t_int_0," ",spearman_p_t_int_1," ",spearman_p_t_int_2," ",kendall_p_t_int_0," ",kendall_p_t_int_1," ",kendall_p_t_int_2)
print()

res = stats.spearmanr(interval_max[T==39], TCP_kill[T==39])
spearman_t_int_0 = res.correlation
spearman_p_t_int_0 = res.pvalue
kendall_t_int_0, kendall_p_t_int_0 = stats.kendalltau(interval_max[T==39], TCP_kill[T==39])

res = stats.spearmanr(interval_max[T==41], TCP_kill[T==41])
spearman_t_int_1 = res.correlation
spearman_p_t_int_1 = res.pvalue
kendall_t_int_1, kendall_p_t_int_1 = stats.kendalltau(interval_max[T==41], TCP_kill[T==41])

res = stats.spearmanr(interval_max[T==43], TCP_kill[T==43])
spearman_t_int_2 = res.correlation
spearman_p_t_int_2 = res.pvalue
kendall_t_int_2, kendall_p_t_int_2 = stats.kendalltau(interval_max[T==43], TCP_kill[T==43])

print("Max t_int Spearman's rho (39 °C, 41 °C, 43 °C): ", spearman_t_int_0,", ",spearman_t_int_1,", ",spearman_t_int_2)
print("Max t_int Kendall's tau (39 °C, 41 °C, 43 °C): ", kendall_t_int_0,", ",kendall_t_int_1,", ",kendall_t_int_2)
print("p-values: ",spearman_p_t_int_0," ",spearman_p_t_int_1," ",spearman_p_t_int_2," ",kendall_p_t_int_0," ",kendall_p_t_int_1," ",kendall_p_t_int_2)
print()
print()


print("Correlations (mu = 0.5 h^-1 without direct HT cell killing)")
res = stats.spearmanr(interval_min[T==39], TCP_high[T==39])
spearman_t_int_0 = res.correlation
spearman_p_t_int_0 = res.pvalue
kendall_t_int_0, kendall_p_t_int_0 = stats.kendalltau(interval_min[T==39], TCP_high[T==39])

res = stats.spearmanr(interval_min[T==41], TCP_high[T==41])
spearman_t_int_1 = res.correlation
spearman_p_t_int_1 = res.pvalue
kendall_t_int_1, kendall_p_t_int_1 = stats.kendalltau(interval_min[T==41], TCP_high[T==41])

res = stats.spearmanr(interval_min[T==43], TCP_high[T==43])
spearman_t_int_2 = res.correlation
spearman_p_t_int_2 = res.pvalue
kendall_t_int_2, kendall_p_t_int_2 = stats.kendalltau(interval_min[T==43], TCP_high[T==43])

print("Min t_int Spearman's rho (39 °C, 41 °C, 43 °C): ", spearman_t_int_0,", ",spearman_t_int_1,", ",spearman_t_int_2)
print("Min t_int Kendall's tau (39 °C, 41 °C, 43 °C): ", kendall_t_int_0,", ",kendall_t_int_1,", ",kendall_t_int_2)
print("p-values: ",spearman_p_t_int_0," ",spearman_p_t_int_1," ",spearman_p_t_int_2," ",kendall_p_t_int_0," ",kendall_p_t_int_1," ",kendall_p_t_int_2)
print()

res = stats.spearmanr(interval[T==39], TCP_high[T==39])
spearman_t_int_0 = res.correlation
spearman_p_t_int_0 = res.pvalue
kendall_t_int_0, kendall_p_t_int_0 = stats.kendalltau(interval[T==39], TCP_high[T==39])

res = stats.spearmanr(interval[T==41], TCP_high[T==41])
spearman_t_int_1 = res.correlation
spearman_p_t_int_1 = res.pvalue
kendall_t_int_1, kendall_p_t_int_1 = stats.kendalltau(interval[T==41], TCP_high[T==41])

res = stats.spearmanr(interval[T==43], TCP_high[T==43])
spearman_t_int_2 = res.correlation
spearman_p_t_int_2 = res.pvalue
kendall_t_int_2, kendall_p_t_int_2 = stats.kendalltau(interval[T==43], TCP_high[T==43])

print("Mean t_int Spearman's rho (39 °C, 41 °C, 43 °C): ", spearman_t_int_0,", ",spearman_t_int_1,", ",spearman_t_int_2)
print("Mean t_int Kendall's tau (39 °C, 41 °C, 43 °C): ", kendall_t_int_0,", ",kendall_t_int_1,", ",kendall_t_int_2)
print("p-values: ",spearman_p_t_int_0," ",spearman_p_t_int_1," ",spearman_p_t_int_2," ",kendall_p_t_int_0," ",kendall_p_t_int_1," ",kendall_p_t_int_2)
print()

res = stats.spearmanr(interval_max[T==39], TCP_high[T==39])
spearman_t_int_0 = res.correlation
spearman_p_t_int_0 = res.pvalue
kendall_t_int_0, kendall_p_t_int_0 = stats.kendalltau(interval_max[T==39], TCP_high[T==39])

res = stats.spearmanr(interval_max[T==41], TCP_high[T==41])
spearman_t_int_1 = res.correlation
spearman_p_t_int_1 = res.pvalue
kendall_t_int_1, kendall_p_t_int_1 = stats.kendalltau(interval_max[T==41], TCP_high[T==41])

res = stats.spearmanr(interval_max[T==43], TCP_high[T==43])
spearman_t_int_2 = res.correlation
spearman_p_t_int_2 = res.pvalue
kendall_t_int_2, kendall_p_t_int_2 = stats.kendalltau(interval_max[T==43], TCP_high[T==43])

print("Max t_int Spearman's rho (39 °C, 41 °C, 43 °C): ", spearman_t_int_0,", ",spearman_t_int_1,", ",spearman_t_int_2)
print("Max t_int Kendall's tau (39 °C, 41 °C, 43 °C): ", kendall_t_int_0,", ",kendall_t_int_1,", ",kendall_t_int_2)
print("p-values: ",spearman_p_t_int_0," ",spearman_p_t_int_1," ",spearman_p_t_int_2," ",kendall_p_t_int_0," ",kendall_p_t_int_1," ",kendall_p_t_int_2)
print()
print()



print("Correlations (mu = 0.5 h^-1 without direct HT cell killing)")
res = stats.spearmanr(interval_min[T==39], TCP_high_kill[T==39])
spearman_t_int_0 = res.correlation
spearman_p_t_int_0 = res.pvalue
kendall_t_int_0, kendall_p_t_int_0 = stats.kendalltau(interval_min[T==39], TCP_high_kill[T==39])

res = stats.spearmanr(interval_min[T==41], TCP_high_kill[T==41])
spearman_t_int_1 = res.correlation
spearman_p_t_int_1 = res.pvalue
kendall_t_int_1, kendall_p_t_int_1 = stats.kendalltau(interval_min[T==41], TCP_high_kill[T==41])

res = stats.spearmanr(interval_min[T==43], TCP_high_kill[T==43])
spearman_t_int_2 = res.correlation
spearman_p_t_int_2 = res.pvalue
kendall_t_int_2, kendall_p_t_int_2 = stats.kendalltau(interval_min[T==43], TCP_high_kill[T==43])

print("Min t_int Spearman's rho (39 °C, 41 °C, 43 °C): ", spearman_t_int_0,", ",spearman_t_int_1,", ",spearman_t_int_2)
print("Min t_int Kendall's tau (39 °C, 41 °C, 43 °C): ", kendall_t_int_0,", ",kendall_t_int_1,", ",kendall_t_int_2)
print("p-values: ",spearman_p_t_int_0," ",spearman_p_t_int_1," ",spearman_p_t_int_2," ",kendall_p_t_int_0," ",kendall_p_t_int_1," ",kendall_p_t_int_2)
print()

res = stats.spearmanr(interval[T==39], TCP_high_kill[T==39])
spearman_t_int_0 = res.correlation
spearman_p_t_int_0 = res.pvalue
kendall_t_int_0, kendall_p_t_int_0 = stats.kendalltau(interval[T==39], TCP_high_kill[T==39])

res = stats.spearmanr(interval[T==41], TCP_high_kill[T==41])
spearman_t_int_1 = res.correlation
spearman_p_t_int_1 = res.pvalue
kendall_t_int_1, kendall_p_t_int_1 = stats.kendalltau(interval[T==41], TCP_high_kill[T==41])

res = stats.spearmanr(interval[T==43], TCP_high_kill[T==43])
spearman_t_int_2 = res.correlation
spearman_p_t_int_2 = res.pvalue
kendall_t_int_2, kendall_p_t_int_2 = stats.kendalltau(interval[T==43], TCP_high_kill[T==43])

print("Mean t_int Spearman's rho (39 °C, 41 °C, 43 °C): ", spearman_t_int_0,", ",spearman_t_int_1,", ",spearman_t_int_2)
print("Mean t_int Kendall's tau (39 °C, 41 °C, 43 °C): ", kendall_t_int_0,", ",kendall_t_int_1,", ",kendall_t_int_2)
print("p-values: ",spearman_p_t_int_0," ",spearman_p_t_int_1," ",spearman_p_t_int_2," ",kendall_p_t_int_0," ",kendall_p_t_int_1," ",kendall_p_t_int_2)
print()

res = stats.spearmanr(interval_max[T==39], TCP_high_kill[T==39])
spearman_t_int_0 = res.correlation
spearman_p_t_int_0 = res.pvalue
kendall_t_int_0, kendall_p_t_int_0 = stats.kendalltau(interval_max[T==39], TCP_high_kill[T==39])

res = stats.spearmanr(interval_max[T==41], TCP_high_kill[T==41])
spearman_t_int_1 = res.correlation
spearman_p_t_int_1 = res.pvalue
kendall_t_int_1, kendall_p_t_int_1 = stats.kendalltau(interval_max[T==41], TCP_high_kill[T==41])

res = stats.spearmanr(interval_max[T==43], TCP_high_kill[T==43])
spearman_t_int_2 = res.correlation
spearman_p_t_int_2 = res.pvalue
kendall_t_int_2, kendall_p_t_int_2 = stats.kendalltau(interval_max[T==43], TCP_high_kill[T==43])

print("Max t_int Spearman's rho (39 °C, 41 °C, 43 °C): ", spearman_t_int_0,", ",spearman_t_int_1,", ",spearman_t_int_2)
print("Max t_int Kendall's tau (39 °C, 41 °C, 43 °C): ", kendall_t_int_0,", ",kendall_t_int_1,", ",kendall_t_int_2)
print("p-values: ",spearman_p_t_int_0," ",spearman_p_t_int_1," ",spearman_p_t_int_2," ",kendall_p_t_int_0," ",kendall_p_t_int_1," ",kendall_p_t_int_2)
print()
