import numpy as np
import matplotlib.pyplot as plt
import timeit

time_start = timeit.default_timer() #starter tiden

g = 9.81 #definerer tyngdeakselerasjonen i m/s^2
R = 0.50 #definerer radien R til kvartsirkelen i m

omega_prior = 0 #definerer vinkelhastigheten ved t=0 
alpha = 0 #definerer vinkelakselerasjonen ved t=0
theta_prior = 6.6*np.pi/180  #definerer startvinkelen ved t=0
delta_t = 0.00001 #kan ikke være 0 for å kunne bruke Eulers metode

objects = {'kule': [2/5, 0.169, 0.057/2], 'kompakt sylinder': [1/2, 1.099, 0.044/2]} #dictionary for objektene vi bruker i denne oppgaven


def angle(objects, object, R, omega_prior, alpha, theta_prior):
        
    """
    Funksjon som tar inn en dictionary med objekter, et gitt objekt, en gitt Radius R, 
    en vinkelfarten, vinkelakselerasjon, startvinkel og en delta t. 
    
    Deretter returnerer funksjonen en streng "func_values" med informasjon om hvilken vinkel
    objektet forlater banen i. Dette kommer i både grader og radianer.
    Vi får også vite antall iterasjoner, samt tiden ballen brukte på å forlate banen.
    Vi printer også translatorisk kinetisk energi, samt rotasjonsenergi ved slutt.
    Vi får også et plot av energien til systemet for vår forståelse.
    """
        
    values = objects.get(object) #lagrer informasjon om et gitt objekt i en variabel
    c = values[0] #henter ut c fra dict
    m = values[1] #henter ut m fra dict
    r = values[2] #henter ut r fra dict
            
    N = 1 #definerer en vilkårlig verdi for N større enn null slik at løkken vil kjøre
    i = 0 #indeks  

    pot_e_arr = np.zeros(2*int(1/delta_t)) #lager et tomt array for potensiell energi
    kin_e_arr = np.zeros(2*int(1/delta_t)) #lager et tomt array for kinetisk energi
    rot_e_arr = np.zeros(2*int(1/delta_t)) #lager et tomt array for rotasjonsenergi 
    tot_e_arr = np.zeros(2*int(1/delta_t)) #lager et tomt array for total energi
    
    I = c * m * r ** 2
    
    while N > 0: #løkke som går til normalkraften er mindre enn 0
        alpha = ((g * np.sin(theta_prior)) / (1 + c)) / (r + R) #vinkelakselerasjon
        omega_post = omega_prior + alpha * delta_t #Euler på vinkelfart
        theta_post = theta_prior + omega_prior * delta_t #vinkelen

        v = (r + R) * omega_prior #farten langs overflaten
        N = m * g * np.cos(theta_prior) - (m * (v ** 2)) / (r + R) #normalkraften

        pot_e_arr[i] = m*g*(R+r)*np.cos(theta_post)
        kin_e_arr[i] = 0.5 * m * v ** 2
        rot_e_arr[i] = 0.5 * I * (v/r) ** 2
        tot_e_arr[i] = pot_e_arr[i] + kin_e_arr[i] + rot_e_arr[i]
        
        i += 1 #legger til en i antall iterasjoner
        
        omega_prior = omega_post #oppdaterer omega for euler
        theta_prior = theta_post #oppdater theta for euler

    time = i*delta_t #finner sluttiden
    
    plot_pot_e = np.delete(pot_e_arr, [*range(i,2*int(1/delta_t))]) #lager nytt array uten nullene fra arrayet som ikke tas i bruk
    plot_kin_e = np.delete(kin_e_arr, [*range(i,2*int(1/delta_t))]) #lager nytt array uten nullene fra arrayet som ikke tas i bruk
    plot_rot_e = np.delete(rot_e_arr, [*range(i,2*int(1/delta_t))])
    plot_tot_e = np.delete(tot_e_arr, [*range(i,2*int(1/delta_t))]) 
    
    plot_time = np.arange(0,time, delta_t) #lager et array med tiden
    
    plt.title(f"Energien til {object}") #plotter tittel
    plt.xlabel("Tid") #plotter navn på x-aksen
    plt.ylabel("Energi i joule") #plotter navn på y-aksen
    plt.plot(plot_time, plot_tot_e, "r-") #plotter x og y arrayene
    plt.plot(plot_time, plot_pot_e, "b-") #plotter x og y arrayene
    plt.plot(plot_time, plot_kin_e, "k-") #plotter x og y arrayene
    plt.plot(plot_time, plot_rot_e, "g-") #plotter x og y arrayene
    
    plt.legend(['Total energi', 'Potensiell energi', 'Translatorisk kinetisk energi', 'Rotasjonsenergi'])    
    
    plt.grid() #legger til grid
    plt.show() #viser plottet
    
    func_values = f"""
Ballen forlot banen etter {theta_post*180/np.pi} grader. 
Løkken kjørte {i} ganger med delta_t: {delta_t}
Ballen forlot banen etter {round(delta_t*i,6)} sekunder.

Objektet mistet {tot_e_arr[1]-plot_tot_e[-1]} joule energi.
Den translatorisk kinetiske energien ved slutten ble {plot_kin_e[i-1]} joule.
Rotasjonsenergien ble {plot_rot_e[i-1]} joule.
"""#lagrer en stor streng så svaret blir printet fint
    return func_values
    #returnerer vinkelen i grader, samt antall kjøringer, normalkraften
    #der objektet mister kontakt med underlaget, og den siste vinkelen registrert

print(angle(objects, "kompakt sylinder", R, omega_prior, alpha, theta_prior)) #kjører og printer funksjonen

time_stop = timeit.default_timer() #slutter tiden

print(f"Progammet brukte {time_stop-time_start} sekunder") #printer tiden programmet brukte for å kjøre