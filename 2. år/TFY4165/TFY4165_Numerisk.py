import numpy as np
import matplotlib.pyplot as plt
plt.plot(0,0)

T = 1000
dt = 1

loppeantall = [6,150,2000]

# Analytisk

def analytisk(N,T, dt):
    t = np.linspace(T,dt)
    N_a = 1/2 * (1 + np.e**(-2*t/N))
    return t, N_a


# Numerisk simulering
def numerisk(N,T,dt):
    steg = int(T/dt)
    print(steg)
    N_A = np.zeros(steg)
    N_A[0] = N
    N_B = np.zeros(steg)

    hopper = np.random.randint(0,N,T)
    
    
    loppested = np.zeros(N)
    
    for i in range(1,T):
        if loppested[hopper[i]] == 0 and N_A[i-1] > 0:
            loppested[hopper[i]] = 1
            N_B[i]=N_B[i-1]+1
            N_A[i]=N_A[i-1]-1
        elif loppested[hopper[i]] == 1 and N_B[i-1] > 0:
            loppested[hopper[i]] = 0
            N_B[i]=N_B[i-1]-1
            N_A[i]=N_A[i-1]+1
        else:
            N_B[i]=N_B[i-1]
            N_A[i]=N_A[i-1]
    t = np.arange(0,T,dt)

    return t, N_A/N

for N in loppeantall:
    T = N*20
    t_analytisk, N_a = analytisk(N, T, dt)
    t_numerisk, N_A = numerisk(N, T, dt)
    plt.plot(0,0)
    plt.plot(t_numerisk, N_A, label = 'Numerisk')
    plt.plot(t_analytisk, N_a, label = 'Analytisk')
    plt.title(f'Loppesimulering N={N}, T={T}')
    plt.ylabel('Na(t)/N')
    plt.xlabel('Tid(s)')
    plt.legend()
    plt.figure()







