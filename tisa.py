# coding=utf-8

import os, sys
import collections


import pandas as pd
#import numpy as np

from openpyxl import Workbook

from openpyxl.utils.dataframe import dataframe_to_rows


#os.environ["JPYPE_JVM"] = "/usr/lib/jvm/java-11-openjdk-amd64/lib/server/libjvm.so"
#os.environ["CLASSPATH"]="/home/ori/TISA"

from cinfony import pybel, rdk
#from cinfony import cdk
import chemfp
from rdkit.ML.Scoring import Scoring

#es_actiu = 1
#allmols = [mol for mol in pybel.readfile("sdf", "/home/ori/TISA/test.sdf")]

#rdkmol = rdk.Molecule(mol)
#pybelmol = pybel.Molecule(rdkmol)

def convertirEnInchi(list_of_molecules):
	
	x=0
	for mol in list_of_molecules:
		list_of_molecules[x] = list_of_molecules[x].write("inchi")
		x+=1
	
	return list_of_molecules
	
def convertirEnSDF(list_of_molecules):
	
	x=0
	
	for mol in list_of_molecules:
		list_of_molecules[x] = pybel.readstring("inchi", str(list_of_molecules[x])).write("sdf")
		x+=1
	
	return list_of_molecules
	
def unique_mols(mols):
	return {mol.write("inchikey"):mol for mol in mols}.values()

def eliminar_repetits(list_of_molecules):
	
	#seen = set()
	#uniq = []
	
	#list_of_molecules = convertirEnInchi(list_of_molecules)
	
	list_of_molecules = list(dict.fromkeys(list_of_molecules))
	
	#list_of_molecules = convertirEnSDF(list_of_molecules)

	return list_of_molecules

def crearLlistaTuple(list_of_molecules, esAct):
	
	#x=0
	#llista = [0] * len(list_of_molecules)
	
	#for mol in list_of_molecules:
	#	llista[x]=(mol, esAct)
	#	x+=1
	return [(mol, esAct) for mol in list_of_molecules]
	#return llista

def trobarMaxims(llista_actius, llista_totals, fp): 
	
	maxims = [[0 for x in range(4)] for y in range(len(llista_totals))]	
	for i in range(len(llista_totals)):
		
		maxims[i][1] = llista_totals[i][1]
		#print maxims[1][1]
		for j in range(len(llista_actius)):
			if (((llista_totals[i][0].calcfp(fptype = fp) | llista_actius[j][0].calcfp(fptype = fp)) > maxims[i][0]) and llista_totals[i][0] is not llista_actius[j][0]):
				maxims[i][0] = llista_totals[i][0].calcfp(fptype = fp) | llista_actius[j][0].calcfp(fptype = fp)
				maxims[i][2] = llista_actius[j][0]	
	return maxims
		
def mateixaMolecula(molecula1, molecula2):
	
	molecula1 = molecula1.write("inchi")
	molecula2 = molecula2.write("inchi")
	
	return molecula1 == molecula2
				
def ordenarTanimotos(llista_maxims):
	
	llista_maxims.sort(reverse = True)
	
	return llista_maxims

def calcularEF(factor, llistaTuplesOrdenada, num_actius):
	
	actius_trobats = 0
	ef_llargada = len(llistaTuplesOrdenada)*factor/100
	
	for i in range(ef_llargada):
		if llistaTuplesOrdenada[i][1]:
			 actius_trobats+=1
	
	ef = 100 * actius_trobats / num_actius
	
	return ef
	
def calcularBEDROC(llistaTuplesOrdenada):
	
	llista_scores = [(1-el[0], el[1]) for el in llistaTuplesOrdenada]
	bedroc = Scoring.CalcBEDROC(llista_scores, 1, 20) 

	return bedroc

	
#actius = [mol for mol in pybel.readfile("sdf", "/home/ori/TISA/Compound_014550001_014575000.sdf")]
#problemes = [mol for mol in pybel.readfile("sdf", "/home/ori/TISA/Compound_027575001_027600000.sdf")]
actius = [mol for mol in pybel.readfile("sdf", "/home/ori/TISA/test.sdf")]#problemes = [mol for mol in pybel.readfile("sdf", "/home/ori/TISA/test.sdf")]

problemes = [mol for mol in pybel.readfile("sdf", "/home/ori/TISA/Compound_027575001_027600000.sdf")]


actius=eliminar_repetits(actius) #NO FUNCIONA
problemes=eliminar_repetits(problemes)

"""
print len(actius)
actius = eliminar_repetits(actius)
problemes = eliminar_repetits(problemes)
print len(actius)
"""


llistaActius = crearLlistaTuple(actius,1)
llistaTotal = crearLlistaTuple(problemes,0)

llistaTotal = llistaActius + llistaTotal


for fp in pybel.fps: #calcular fps al principi de cada iteració, adaptar el codi (trobarMaxims). També funcionarà per rdk, però s'han de transformar les molècules (mirar a dalt els comentaris)ss
	
	
	#fp_llista_total = chain(fp_actius, fp_problemes)
	
	maxims = trobarMaxims(llistaActius, llistaTotal, fp) #Càlcul de la matriu amb tots els Tanimotos

	maxims = ordenarTanimotos(maxims)

	df = pd.DataFrame(maxims, columns =['Molècula','Tanimoto', 'És Actiu', 'Actiu més semblant'])
	
	df.to_csv(r'/home/ori/TISA/Resultats/'+fp+'.csv') #Afegir el nom del fp
	print(fp+" COMPLETAT") #fer fp/len(fps)*100 completat
	
	#FER UN ARXIU RESULTAT AMB EFs I BEDROC PER CADA fp. La única pestanya de moment seria el Tanimoto
	
	data = [('EF1%', calcularEF(1,maxims,len(llistaActius))), ('EF10%', calcularEF(10,maxims,len(llistaActius))), ('BEDROC', calcularBEDROC(maxims))]
	
	csvfile = open(fp+'.csv', 'wb')
	
	"""
	for row in data:
		line = ','.join(str(row))
		csvfile.write(line + '\n')
	

	
	print "EF1%:\t"+str(calcularEF(1,maxims,len(llistaActius)))
	print "EF10%:\t"+str(calcularEF(10,maxims,len(llistaActius)))
	print "BEDROC:\t"+str(calcularBEDROC(maxims))
	

	
	fp_metriques = Workbook()
	dest_filename1 = "fp_metriques.xlsx"
	
	
	ws1 = fp_metriques.active
	ws1.title = fp_list[fp_index]
	
	for r in dataframe_to_rows(df, index = False, header = True):
		ws1.append(r)	#ValueError: Cannot convert <cinfony.pybel.Molecule object at 0x7f72f0b08210> to Excel
		
	for cell in ws1['A'] + ws1[1]:
		cell.style = 'Pandas'
	""" 
