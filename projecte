# Com detectar si una molècula és activa? --> Establir un límit de Tanimoto
# Comparar actius amb actius i actius amb problemes, fer una llista total amb tots els Tanimotos i ordenar-la juntament amb una de booleans es_actiu
# Després 

""" from openbabel import PyBEL """
from rdkit import Chem
import pybel
from cinfony import rdk

print ("\nalfa\n")

""" print "*NUM. ATOMS*\n===============" """
""" for mol in suppl:
    print(mol.GetNumAtoms()) """

allmols = [mol for mol in pybel.readfile("sdf", "/Users/oriol/py-workspace/projecte/test.sdf")]
actius = [mol for mol in rdk.readfile("sdf", "/Users/oriol/py-workspace/projecte/test.sdf")]
problemes = [mol for mol in rdk.readfile("sdf", "/Users/oriol/py-workspace/projecte/test.sdf")]

""" for mymol in allmols:
    print(mymol.data) """


i, j = 0, 0

""" print ("\n\n*TANIMOTOS*\n===============")

k = 0

while k < len(allmols):
    fp = actius[k].calcfp()
    print fp
    print actius[k]
    k += 1 """

matriz = [[actius[i].calcfp()|problemes[j].calcfp() for i in range(len(actius))] for j in range(len(problemes))] #Càlcul de la matriu amb tots els Tanimotos

maxims = [0] * len(problemes) #Inicialització del string que inclourà els màxims
results = []
tipus = []

l, m = 0, 0,

while l < len(problemes): #Càlcul del Tanimoto més gran de cada molècula, excepte si és ell mateix (quan és igual a 1) i si passa del límit establert (és actiu)
    m = 0
    while m < len(actius):

        if(maxims[l] < matriz[l][m]) and matriz[l][m] != 1.0:
            maxims[l] = matriz[l][m]
    
        m += 1
    l += 1

print ("\nNO ORDENATS\n", maxims, "\n", problemes)

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

    o += 1

print ("\nORDENATS:\n", maxims, "\n", problemes)

#relacions = dict(zip(maxims, actius))

print ("\nomega\n")