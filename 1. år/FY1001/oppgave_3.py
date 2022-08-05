import numpy as np
import matplotlib.pyplot as plt
import timeit


time_start = timeit.default_timer()


objects = {'hul sylinder': [1.0, 0.254, [0.021, 0.018]], 'kompakt sylinder': [1/2, 1.099, [0.044/2]]} 
#dictionary for objektene vi bruker i denne oppgaven

R = 0.5 #radien til kvartsirkelen
omega_prior = 0 #definerer vinkelfarten  
alpha = 0 #definerer vinkelakselerasjonen 
theta_prior = 2.9 * np.pi / 180 #definerer startvinkelen
delta_t = 0.00001 #kan ikke være 0 for å kunne bruke Eulers metode

def angle(objects, object, R, omega_prior, alpha, theta_prior, delta_t):
        
    """
    Funksjon som tar inn en dictionary med objekter, et gitt objekt, en gitt Radius R, 
    en vinkelfarten, vinkelakselerasjon, startvinkel og en delta t. 
    
    Deretter returnerer funksjonen en streng "func_values" med informasjon om hvilken vinkel
    objektet forlater banen i. Dette kommer i både grader og radianer.
    Vi får også vite antall iterasjoner, samt tiden ballen brukte på å forlate banen.
    Vi printer også translatorisk kinetisk energi, samt rotasjonsenergi ved slutt.
    Vi får også et plot av energien til systemet og punktet der sluring intreffer.
    Det er ikke plot av rotasjonsenergien gjennom bevegelsen.
    """
    
    g = 9.81 #definerer gravitasjonskraften i m/s^2

    values = objects.get(object) #lagrer informasjon om et gitt objekt i en variabel
    c = values[0] #henter ut c fra dict
    m = values[1] #henter ut m fra dict
    r = sum(values[2]) # henter ut r fra dict
    
    mu_static = 0.61 #definerer den statiske friksjonskoeffisienten
    mu_kinetic = 0.47 #definerer den kinetiske friksjonskoeffisienten
    
    
    if "hul sylinder" == object: #tar hensyn til at treghetsmomentet er ulikt for hul sylinder
        r_I = 0 
        for r_n in values[2]:
            r_I += r_n ** 2
        I = c * m * r_I #treghetsmomentet hvis vi har hul sylinder  

    else: #treghetsmoment på alle andre måter
        r_I = sum(values[2]) #henter ut radius
        I = c * m * r_I ** 2



    pot_e_arr = np.zeros(2*int(1/delta_t)) #lager et tomt array for potensiell energi
    kin_e_arr = np.zeros(2*int(1/delta_t)) #lager et tomt array for kinetisk energi
    # rot_e_arr = np.zeros(2*int(1/delta_t)) #lager et tomt array for rotasjonsenergi 
    tot_e_arr = np.zeros(2*int(1/delta_t)) #lager et tomt array for total energi
    
    i = 0 #setter i for å telle antall iterasjoner
    
    f_s = 1 #definerer en f_s som er mindre enn f_s_max så løkken kjører
    f_s_max = 2 #definierer en f_s_max som er større enn f_s så løkken kjører
    
    N = 1 #definierer en vilkårlig N som er større enn 0 så løkken starter
    #løkke for ren rulling
    while N > 0 and f_s <= f_s_max: #løkke som kjører mens normalkraft er større enn 0, og friksjonen er mindre enn maksfriksjon
      
        alpha = ((g * np.sin(theta_prior)) / (1 + c)) / (r + R) #regner ut vinkelakselerasjon
        omega_post = omega_prior + alpha * delta_t #bruker euler på vinkelfart
        theta_post = theta_prior + omega_prior * delta_t #bruker euler på vinkelen

        
        v = (r + R) * omega_prior #regner ut farten
        N = m * g * np.cos(theta_prior) - (m * (v ** 2)) / (r + R) #normalkraften
        
        pot_e_arr[i] = m*g*(R+r)*np.cos(theta_post) #regner ut potensiell energi og legger den i sitt array
        kin_e_arr[i] = 0.5 * m * v ** 2 #regner ut kinetisk energi og legger den i sitt array
        # rot_e_arr[i] = 0.5 * I * (v/r) ** 2
        tot_e_arr[i] = pot_e_arr[i] + kin_e_arr[i] #+ rot_e_arr[i] #regner ut total energi og legger den i sitt array
        
        
        
        i += 1 #legger til en iterasjon i tellingen
        omega_prior = omega_post #oppdaterer omega
        theta_prior = theta_post #oppdaterer theta
        f_s = c*m*alpha*(r+R) #regner ut ny friksjonskraft
        f_s_max = mu_static * N #bruker den nye normalkraften til å finne friksjonskraft
        
    break_index = i #definerer en break index når den første løkken bryter for å finne ut hvor sluring inntreffer
             
    f_k = f_s
            
    while N > 0: #løkke som kjører mens normalkraften er større enn 0
        alpha = ((g * np.sin(theta_prior)) - (f_k / m)) / (r + R) #regner ut vinkelakselerasjon
        omega_post = omega_prior + alpha * delta_t #bruker euler på vinkelfart
        theta_post = theta_prior + omega_prior * delta_t #bruker euler på vinkelen
        
        v = (r + R) * omega_prior #regner ut hastighet
        N = m * g * np.cos(theta_prior) - (m * (v ** 2)) / (r + R) #normalkraften
        
        pot_e_arr[i] = m*g*(R+r)*np.cos(theta_post) #potensiell energi
        kin_e_arr[i] = 0.5 * m * v ** 2 #kinetisk energi
        # rot_e_arr[i] = rot_e_end
        tot_e_arr[i] = pot_e_arr[i] + kin_e_arr[i] #total energi
        
        i += 1 #legger til en iterasjon i tellingen
        omega_prior = omega_post #oppdaterer omega
        theta_prior = theta_post #oppdaterer theta
        f_k = N*mu_kinetic #regner ut ny friksjonskraft

        

        
    time = i*delta_t #finner sluttiden
    
    plot_pot_e = np.delete(pot_e_arr, [*range(i,2*int(1/delta_t))]) #lager nytt array uten nullene fra arrayet som ikke tas i bruk
    plot_kin_e = np.delete(kin_e_arr, [*range(i,2*int(1/delta_t))]) #lager nytt array uten nullene fra arrayet som ikke tas i bruk
    plot_tot_e = np.delete(tot_e_arr, [*range(i,2*int(1/delta_t))]) 
    
    plot_time = np.arange(0,time, delta_t) #lager et array med tiden
        
    plt.title(f"Energien til {object}") #plotter tittel
    plt.xlabel("Tid") #plotter navn på x-aksen
    plt.ylabel("Energi i joule") #plotter navn på y-aksen
    plt.plot(plot_time, plot_tot_e, "r-") #plotter x og y arrayene
    plt.plot(plot_time, plot_pot_e, "b-") #plotter x og y arrayene
    plt.plot(plot_time, plot_kin_e, "k-") #plotter x og y arrayene
    # plt.plot(plot_time_rot, plot_rot_e, "g-") #plotter x og y arrayene
    
    plt.legend(['Total energi', 'Potensiell energi', 'Translatorisk kinetisk energi', 'Rotasjonsenergi'])    

    plt.plot(break_index*delta_t, tot_e_arr[break_index-1], 'bo') #viser punktet rett før sluring intreffer
    plt.grid() #legger til grid
    plt.show() #viser plottet
    
    rot_e_end = 0.5 * I * (f_k * N * r * i * delta_t / I) ** 2 #regner ut rotasjonsenergien til ballen ved slutten
    
    tot_e_end = plot_kin_e[-1] + plot_pot_e[-1] + rot_e_end #regner ut den totale energien ved slutten
    
    
    
    func_values = f"""
{object} forlot banen etter {theta_post*180/np.pi} grader,
det er {theta_post} radianer. 

Løkken kjørte {i} ganger med delta_t: {delta_t}
{object} forlot banen etter {round(delta_t*i,6)} sekunder.

Objektet mistet {tot_e_arr[1]-tot_e_end} joule energi.
Den translatorisk kinetiske energien ved slutten ble {plot_kin_e[-1]} joule.
Rotasjonsenergien ved slutten ble {rot_e_end} joule.

""" #lagrer en stor streng så svaret blir printet fint
    return func_values #returnerer strengen

print(angle(objects, "kompakt sylinder", R, omega_prior, alpha, theta_prior, delta_t)) #kjører og printer funksjonen
time_stop = timeit.default_timer() #stopper tiden

print(f"Programmet brukte {round(time_stop-time_start, 5)} sekunder") #printer tiden programmet brukte

