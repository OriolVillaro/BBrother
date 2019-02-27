# coding=utf-8

import pybel
import chemfp

es_actiu = 1
#allmols = [mol for mol in pybel.readfile("sdf", "/home/ori/TISA/test.sdf")]

#actius = [mol for mol in pybel.readfile("sdf", "/home/ori/TISA/test.sdf")]
#problemes = [mol for mol in pybel.readfile("sdf", "/home/ori/TISA/test.sdf")]

for mol in pybel.readfile("sdf", "/home/ori/TISA/test.sdf"):
	
	actius=(mol,es_actiu)


es_actiu = 0

for mol in pybel.readfile("sdf", "/home/ori/TISA/test.sdf"):
	
	problemes=(mol,es_actiu)

print(actius, problemes)	
"""i, j = 0, 0

matriz = [[actius[i].calcfp()|problemes[j].calcfp() for i in range(len(actius))] for j in range(len(problemes))] #Càlcul de la matriu amb tots els Tanimotos

maxims = [0] * len(problemes) #Inicialització del string que inclourà els màxims

l, m = 0, 0,

while l < len(problemes): #Càlcul del Tanimoto més gran de cada molècula, excepte si és ell mateix (quan és igual a 1) i si passa del límit establert (és actiu)
    m = 0
    while m < len(actius):

        if(maxims[l] < matriz[l][m]) and matriz[l][m] != 1.0:
            maxims[l] = matriz[l][m]
    
        m += 1
    l += 1

n, o = 0, 0

while o < len(maxims): #Ordenació dels Tanimotos de més gran a més petit mantenint l'ordre de la molècula corresponent
    n = o 
    while n < len(maxims):
       
        if(maxims[o] < maxims[n]):

            aux = maxims[o]
            maxims[o] = maxims[n]
            maxims[n] = aux

            aux1 = problemes[o]
            problemes[o] = problemes[n]
            problemes[n] = aux1

        n +=1      

    o += 1"""



