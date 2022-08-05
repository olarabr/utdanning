import numpy as np
import matplotlib.pyplot as plt
import timeit

time_start = timeit.default_timer() #starter tiden


m = 0.169 #definerer massen til kula i kg
R = 0.50 #definerer radien R til kvartsirkelen i m
r = 0.057/2 #definerer radien r til kula i m

start_omega = 0 #definerer vinkelfarten ved t=0
start_alpha = 0 #definerer vinkelakselerasjonen ved t=0
start_theta = 2.9*np.pi/180   #definerer startvinkelen ved t=0
delta_t = 0.00001 #kan ikke være 0 for å kunne bruke Eulers metode


def angle(m, R, r, omega_prior, alpha, theta_prior, delta_t):
    """
    Funksjon som tar inn masse, en gitt Radius på kvartsirkelen R, radiusen til et objekt r, 
    en vinkelfart, vinkelakselerasjon, startvinkel og tidsintervall delta_t.
    
    Deretter returnerer funksjonen en streng "func_values" med informasjon om hvilken vinkel
    objektet forlater banen i. Dette kommer i både grader og radianer.
    Vi får også vite antall iterasjoner, samt tiden ballen brukte på å forlate banen.
    Vi får også et plot av energien til systemet for vår forståelse.
    
    """
    g = 9.81 #definerer tyngdeakselerasjonen i m/s^2
    
    pot_e_arr = np.zeros(2*int(1/delta_t)) #lager et tomt array for potensiell energi
    kin_e_arr = np.zeros(2*int(1/delta_t)) #lager et tomt array for kinetisk energi
    tot_e_arr = np.zeros(2*int(1/delta_t)) #lager et tomt array for total energi
  
    N = 1 #definerer en vilkårlig verdi for N større enn null slik at løkken vil kjøre
    i = 0 #indeks    
    while N > 0: #while løkke som itererer frem til N blir mindre enn 0. 
        alpha = g * np.sin(theta_prior) / (r + R) #regner ut vinkelakselerasjonen
        omega_post = omega_prior + alpha * delta_t #euler på vinkelfarten
        theta_post = theta_prior + omega_prior * delta_t #bruker euler for å finne ny vinkel

        
        v = (r + R) * omega_prior #farten langs overflaten
        N = m * g * np.cos(theta_prior) - (m * (v ** 2)) / (r + R) #regner ut normalkraften
        
        pot_e_arr[i] = m*g*(R+r)*np.cos(theta_post) #regner ut potensiell energi og legger den i sitt array
        kin_e_arr[i] = 0.5 * m * v ** 2 #regner ut kinetisk energi og legger den i sitt array
        tot_e_arr[i] = pot_e_arr[i] + kin_e_arr[i] #regner ut total energi og legger den i sitt array
        
        omega_prior = omega_post #oppdaterer omega før neste iterasjon
        theta_prior = theta_post #oppdater theta før neste iterasjoner
        
        i += 1 #legger til en i iterasjonstellingen

    plot_pot_e = np.delete(pot_e_arr, [*range(i,2*int(1/delta_t))]) #lager nytt array uten nullene fra arrayet som ikke tas i bruk
    plot_kin_e = np.delete(kin_e_arr, [*range(i,2*int(1/delta_t))]) #lager nytt array uten nullene fra arrayet som ikke tas i bruk
    plot_tot_e = np.delete(tot_e_arr, [*range(i,2*int(1/delta_t))]) #lager nytt array uten nullene fra arrayet som ikke tas i bruk
    
    time = i*delta_t #definerer sluttiden
    plot_time = np.arange(0,time, delta_t) #lager et tidsarray
    
    plt.title("Energien til kulen") #plotter tittel
    plt.xlabel("Tid") #plotter navn på x-aksen
    plt.ylabel("Energi i joule") #plotter navn på y-aksen
    plt.plot(plot_time, plot_tot_e, "r-") #plotter total energi
    plt.plot(plot_time, plot_pot_e, "b-") #plotter potensiell energi
    plt.plot(plot_time, plot_kin_e, "k-") #plotter kinetisk energi
    
    plt.legend(['Total energi', 'Potensiell energi', 'Translatorisk kinetisk energi'])    
    
    plt.grid() #legger til grid
    plt.show() #viser plottet
        
    func_values = f"""
Objektet forlot banen etter {theta_post*180/np.pi} grader,
det er {theta_post} radianer. 
Løkken kjørte {i} ganger med delta_t: {delta_t}
Objektet forlot banen etter {round(delta_t*i,6)} sekunder.
"""#lagrer en stor streng så svaret blir printet fint
    return func_values #returnerer strengen med verdiene

    
print(angle(m,R,r,start_omega,start_alpha, start_theta, delta_t)) #kjører og printer funksjonen

time_stop = timeit.default_timer() #slutter tiden

print(f"Progammet brukte {time_stop-time_start} sekunder") #printer tiden programmet brukte for å kjøre