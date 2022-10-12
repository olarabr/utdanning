"""
Analyse av pengeforbruk i kantinen
"""
import csv
from matplotlib import pyplot as plt
import numpy as np

oversikt = []

with open("transaksjonsliste.csv", 'r') as file:
    csvreader = csv.reader(file, delimiter=';') 
    for row in csvreader:
        dato = row[0]
        sted = row[1]
        beløp = row[2] 
        totalt = [dato,sted,beløp]
        oversikt.append(totalt)
        


oversikt_sit = []
for element in oversikt:
    if 'sit' in element[1].lower():
        element[2] = element[2].replace(',','.')
        oversikt_sit.append(element)

dict = {}

for kjøp in oversikt_sit:
    if kjøp[0] in dict:
        beløp = dict.get(kjøp[0])
        dict[kjøp[0]] = beløp + float(kjøp[2])
    else:
        dict[kjøp[0]] = float(kjøp[2])
    

dato = my_list = [i for i in dict.keys()]        
beløp = [i for i in dict.values()]
   
plt.bar(range(len(dict)),beløp, tick_label = dato)
plt.xticks(rotation=270)

for i in range(len(dict)):
    plt.text(i,beløp[i],round(beløp[i],1))

plt.title('Forbruk per dato')
plt.figure()


sumutvikling = np.zeros(len(beløp))
for i in range(len(beløp)):
    sumutvikling[i] = sum(beløp[:i])

plt.bar(range(len(dict)),sumutvikling, tick_label = dato)
plt.xticks(rotation=270)

for i in range(len(dict)):
    text = f'{round((sumutvikling[i]/1000),2)}k'
    plt.text(i,sumutvikling[i],text)
plt.title('Kumulativt forbruk')
