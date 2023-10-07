import numpy as np
from math import exp, log
import random
from scipy import stats


################################
#                              #
#   Temperature dependencies   #
#                              #
################################


# Read results with uniform distributions
data = np.loadtxt('results/TCP_T.dat',skiprows=1)
TCP_uniform = data[:,3]
TCP_kill_uniform = data[:,4]

range_uniform = np.max(TCP_uniform) - np.min(TCP_uniform)
range_kill_uniform = np.max(TCP_kill_uniform) - np.min(TCP_kill_uniform)


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


random.seed(10)   # Random generator seed to obtaion always the same results

# Loop to calculate the TCP values
for n in range(0, N_patients):
    # Definition of initial number of cells
    N=volume/(2e-9)
    N_kill=volume/(2e-9)
    
    for i in range(0, N_fract):
        if i<N_HT:
            # Random selection of temperature
            temp_achi = random.gauss(40.46, 0.68)
            while temp_achi<0:
                temp_achi = random.gauss(40.46, 0.68)
            
            
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

    
range_gauss = np.max(TCP) - np.min(TCP)
range_kill_gauss = np.max(TCP_kill) - np.min(TCP_kill)
print('Temperature dependence')
print('TCP range variation (without direct HT cell killing): ',(range_gauss - range_uniform)*100,' %')
print('TCP range variation (with direct HT cell killing): ',(range_kill_gauss - range_kill_uniform)*100,' %')
print()


##################################
#                                #
#   Time interval dependencies   #
#                                #
##################################


# Read results with uniform distributions
data = np.loadtxt('results/TCP_t_int.dat',skiprows=1)
T = data[:,3]
TCP_uniform = data[:,4]
TCP_high_uniform = data[:,5]

range_uniform_39 = np.max(TCP_uniform[T==39]) - np.min(TCP_uniform[T==39])
range_high_uniform_39 = np.max(TCP_high_uniform[T==39]) - np.min(TCP_high_uniform[T==39])

range_uniform_41 = np.max(TCP_uniform[T==41]) - np.min(TCP_uniform[T==41])
range_high_uniform_41 = np.max(TCP_high_uniform[T==41]) - np.min(TCP_high_uniform[T==41])

range_uniform_43 = np.max(TCP_uniform[T==43]) - np.min(TCP_uniform[T==43])
range_high_uniform_43 = np.max(TCP_high_uniform[T==43]) - np.min(TCP_high_uniform[T==43])

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
TCP = np.zeros(3*N_patients)                    # TCP vector (mu = 0.027 h^-1) 
TCP_high = np.zeros(3*N_patients)               # TCP vector (mu = 0.5 h^-1) 


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
        
    for i in range(0,N_fract):
        if i<N_HT:
            # Random selection of time interval
            t_int = random.gauss(61., 45.)
            while t_int<0:
                t_int = random.gauss(61., 45.)
            
            # Calculation of survival after HT+RT fraction case 1 (mu = 0.027 h^-1)
            HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*3600)
            alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.027/60)
            N=N*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
            # Calculation of survival after HT+RT fraction case 2 (mu = 0.027 h^-1)
            HT_killing=np.exp(-(2.05e10)*(temp_achi_1+273.15)*exp( S/2 - H/(2*(273.15+temp_achi_1)))*3600)
            alpha_T, beta_T = LQ_T(temp_achi_1, t_int, 0.027/60)
            N1=N1*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
            # Calculation of survival after HT+RT fraction case 3 (mu = 0.027 h^-1)
            HT_killing=np.exp(-(2.05e10)*(temp_achi_2+273.15)*exp( S/2 - H/(2*(273.15+temp_achi_2)))*3600)
            alpha_T, beta_T = LQ_T(temp_achi_2, t_int, 0.027/60)
            N2=N2*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
            # Calculation of survival after HT+RT fraction case 1 (mu = 0.5 h^-1)
            HT_killing=np.exp(-(2.05e10)*(temp_achi+273.15)*exp( S/2 - H/(2*(273.15+temp_achi)))*3600)
            alpha_T, beta_T = LQ_T(temp_achi, t_int, 0.5/60)
            N_high=N_high*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
            # Calculation of survival after HT+RT fraction case 2 (mu = 0.5 h^-1)
            HT_killing=np.exp(-(2.05e10)*(temp_achi_1+273.15)*exp( S/2 - H/(2*(273.15+temp_achi_1)))*3600)
            alpha_T, beta_T = LQ_T(temp_achi_1, t_int, 0.5/60)
            N1_high=N1_high*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
            # Calculation of survival after HT+RT fraction case 3 (mu = 0.5 h^-1)
            HT_killing=np.exp(-(2.05e10)*(temp_achi_2+273.15)*exp( S/2 - H/(2*(273.15+temp_achi_2)))*3600)
            alpha_T, beta_T = LQ_T(temp_achi_2, t_int, 0.5/60)
            N2_high=N2_high*exp(-alpha_T*D-beta_T*D*D)*HT_killing
            
        else:
            # Calculation of survival after RT fraction 
            N=N*exp(-alpha_37*D-beta_37*D*D)
            N1=N1*exp(-alpha_37*D-beta_37*D*D)
            N2=N2*exp(-alpha_37*D-beta_37*D*D)
            N_high=N_high*exp(-alpha_37*D-beta_37*D*D)
            N1_high=N1_high*exp(-alpha_37*D-beta_37*D*D)
            N2_high=N2_high*exp(-alpha_37*D-beta_37*D*D)
            
        
        
    # Record temperature    
    T[3*n] = 39
    T[3*n+1] = 41
    T[3*n+2] = 43
    
    # Calculation of TCPs after the HT+RT treatment (mu = 0.027 h^-1)
    TCP[3*n]=exp(-N)  
    TCP[3*n+1]=exp(-N1) 
    TCP[3*n+2]=exp(-N2)
    
    # Calculation of TCPs after the HT+RT treatment (mu = 0.5 h^-1)
    TCP_high[3*n]=exp(-N_high) 
    TCP_high[3*n+1]=exp(-N1_high) 
    TCP_high[3*n+2]=exp(-N2_high)


range_gauss_39 = np.max(TCP[T==39]) - np.min(TCP[T==39])
range_high_gauss_39 = np.max(TCP_high[T==39]) - np.min(TCP_high[T==39])

range_gauss_41 = np.max(TCP[T==41]) - np.min(TCP[T==41])
range_high_gauss_41 = np.max(TCP_high[T==41]) - np.min(TCP_high[T==41])

range_gauss_43 = np.max(TCP[T==43]) - np.min(TCP[T==43])
range_high_gauss_43 = np.max(TCP_high[T==43]) - np.min(TCP_high[T==43])

print('Time interval dependence')
print('TCP range variation (low mu and 39 °C): ',(range_gauss_39 - range_uniform_39)*100,' %')
print('TCP range variation (low mu and 41 °C): ',(range_gauss_41 - range_uniform_41)*100,' %')
print('TCP range variation (low mu and 43 °C): ',(range_gauss_43 - range_uniform_43)*100,' %')

print()
print('TCP range variation (high mu and 39 °C): ',(range_high_gauss_39 - range_high_uniform_39)*100,' %')
print('TCP range variation (high mu and 41 °C): ',(range_high_gauss_41 - range_high_uniform_41)*100,' %')
print('TCP range variation (high mu and 43 °C): ',(range_high_gauss_43 - range_high_uniform_43)*100,' %')
print()


#####################################
#                                   #
#   Random generation of patients   #
#                                   #
#####################################


# Read results with uniform distributions
data = np.loadtxt('results/TCP_random.dat',skiprows=1)
treatments = data[:,0]
TCP_uniform = data[:,3]
TCP_kill_uniform = data[:,4]
TCP_high_uniform = data[:,5]
TCP_kill_high_uniform = data[:,6]

range_uniform_3 = np.max(TCP_uniform[treatments==3]) - np.min(TCP_uniform[treatments==3])
range_uniform_4 = np.max(TCP_uniform[treatments==4]) - np.min(TCP_uniform[treatments==4])
range_uniform_5 = np.max(TCP_uniform[treatments==5]) - np.min(TCP_uniform[treatments==5])
range_uniform_6 = np.max(TCP_uniform[treatments==6]) - np.min(TCP_uniform[treatments==6])

range_kill_uniform_3 = np.max(TCP_kill_uniform[treatments==3]) - np.min(TCP_kill_uniform[treatments==3])
range_kill_uniform_4 = np.max(TCP_kill_uniform[treatments==4]) - np.min(TCP_kill_uniform[treatments==4])
range_kill_uniform_5 = np.max(TCP_kill_uniform[treatments==5]) - np.min(TCP_kill_uniform[treatments==5])
range_kill_uniform_6 = np.max(TCP_kill_uniform[treatments==6]) - np.min(TCP_kill_uniform[treatments==6])

range_high_uniform_3 = np.max(TCP_high_uniform[treatments==3]) - np.min(TCP_high_uniform[treatments==3])
range_high_uniform_4 = np.max(TCP_high_uniform[treatments==4]) - np.min(TCP_high_uniform[treatments==4])
range_high_uniform_5 = np.max(TCP_high_uniform[treatments==5]) - np.min(TCP_high_uniform[treatments==5])
range_high_uniform_6 = np.max(TCP_high_uniform[treatments==6]) - np.min(TCP_high_uniform[treatments==6])

range_kill_high_uniform_3 = np.max(TCP_kill_high_uniform[treatments==3]) - np.min(TCP_kill_high_uniform[treatments==3])
range_kill_high_uniform_4 = np.max(TCP_kill_high_uniform[treatments==4]) - np.min(TCP_kill_high_uniform[treatments==4])
range_kill_high_uniform_5 = np.max(TCP_kill_high_uniform[treatments==5]) - np.min(TCP_kill_high_uniform[treatments==5])
range_kill_high_uniform_6 = np.max(TCP_kill_high_uniform[treatments==6]) - np.min(TCP_kill_high_uniform[treatments==6])


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
                t_int = random.gauss(61., 45.)
                while t_int<0:
                    t_int = random.gauss(61., 45.)    
                
                # Random selection of temperature
                temp_achi = random.gauss(40.46, 0.68)
                while temp_achi<0:
                    temp_achi = random.gauss(40.46, 0.68)
                
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


range_gauss_3 = np.max(TCP[treatments==3]) - np.min(TCP[treatments==3])
range_gauss_4 = np.max(TCP[treatments==4]) - np.min(TCP[treatments==4])
range_gauss_5 = np.max(TCP[treatments==5]) - np.min(TCP[treatments==5])
range_gauss_6 = np.max(TCP[treatments==6]) - np.min(TCP[treatments==6])

range_kill_gauss_3 = np.max(TCP_kill[treatments==3]) - np.min(TCP_kill[treatments==3])
range_kill_gauss_4 = np.max(TCP_kill[treatments==4]) - np.min(TCP_kill[treatments==4])
range_kill_gauss_5 = np.max(TCP_kill[treatments==5]) - np.min(TCP_kill[treatments==5])
range_kill_gauss_6 = np.max(TCP_kill[treatments==6]) - np.min(TCP_kill[treatments==6])

range_high_gauss_3 = np.max(TCP_high[treatments==3]) - np.min(TCP_high[treatments==3])
range_high_gauss_4 = np.max(TCP_high[treatments==4]) - np.min(TCP_high[treatments==4])
range_high_gauss_5 = np.max(TCP_high[treatments==5]) - np.min(TCP_high[treatments==5])
range_high_gauss_6 = np.max(TCP_high[treatments==6]) - np.min(TCP_high[treatments==6])

range_kill_high_gauss_3 = np.max(TCP_kill_high[treatments==3]) - np.min(TCP_kill_high[treatments==3])
range_kill_high_gauss_4 = np.max(TCP_kill_high[treatments==4]) - np.min(TCP_kill_high[treatments==4])
range_kill_high_gauss_5 = np.max(TCP_kill_high[treatments==5]) - np.min(TCP_kill_high[treatments==5])
range_kill_high_gauss_6 = np.max(TCP_kill_high[treatments==6]) - np.min(TCP_kill_high[treatments==6])
    
print('Random patient generation')
print('TCP range variation (low mu, without direct HT cell killing and 3 HT sessions): ',(range_gauss_3 - range_uniform_3)*100,' %')
print('TCP range variation (low mu, without direct HT cell killing and 4 HT sessions): ',(range_gauss_4 - range_uniform_4)*100,' %')
print('TCP range variation (low mu, without direct HT cell killing and 5 HT sessions): ',(range_gauss_5 - range_uniform_5)*100,' %')
print('TCP range variation (low mu, without direct HT cell killing and 6 HT sessions): ',(range_gauss_6 - range_uniform_6)*100,' %')
print()

print('TCP range variation (high mu, without direct HT cell killing and 3 HT sessions): ',(range_high_gauss_3 - range_high_uniform_3)*100,' %')
print('TCP range variation (high mu, without direct HT cell killing and 4 HT sessions): ',(range_high_gauss_4 - range_high_uniform_4)*100,' %')
print('TCP range variation (high mu, without direct HT cell killing and 5 HT sessions): ',(range_high_gauss_5 - range_high_uniform_5)*100,' %')
print('TCP range variation (high mu, without direct HT cell killing and 6 HT sessions): ',(range_high_gauss_6 - range_high_uniform_6)*100,' %')
print()

print('TCP range variation (low mu, with direct HT cell killing and 3 HT sessions): ',(range_kill_gauss_3 - range_kill_uniform_3)*100,' %')
print('TCP range variation (low mu, with direct HT cell killing and 4 HT sessions): ',(range_kill_gauss_4 - range_kill_uniform_4)*100,' %')
print('TCP range variation (low mu, with direct HT cell killing and 5 HT sessions): ',(range_kill_gauss_5 - range_kill_uniform_5)*100,' %')
print('TCP range variation (low mu, with direct HT cell killing and 6 HT sessions): ',(range_kill_gauss_6 - range_kill_uniform_6)*100,' %')
print()

print('TCP range variation (high mu, with direct HT cell killing and 3 HT sessions): ',(range_kill_high_gauss_3 - range_kill_high_uniform_3)*100,' %')
print('TCP range variation (high mu, with direct HT cell killing and 4 HT sessions): ',(range_kill_high_gauss_4 - range_kill_high_uniform_4)*100,' %')
print('TCP range variation (high mu, with direct HT cell killing and 5 HT sessions): ',(range_kill_high_gauss_5 - range_kill_high_uniform_5)*100,' %')
print('TCP range variation (high mu, with direct HT cell killing and 6 HT sessions): ',(range_kill_high_gauss_6 - range_kill_high_uniform_6)*100,' %')
print()
