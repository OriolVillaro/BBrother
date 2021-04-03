# -*- coding: utf-8 -*-

import sys
from  rdkit.ML.Scoring import Scoring
import chemfp
from chemfp import bitops
import pandas as pd
import sklearn
from sklearn import metrics
from openbabel import pybel
import multiprocessing
from multiprocessing import Pool


pybel.ob.obErrorLog.SetOutputLevel(-1)

fptypes = (
    'RDKit-Pattern', 'OpenBabel-MACCS', 'RDKit-Avalon',
'RDKit-AtomPair', 'RDKit-Fingerprint', 
 'RDKit-Torsion',
 'ChemFP-Substruct-RDKit', 
'RDMACCS-OpenBabel', 
'RDMACCS-RDKit', 
'RDKit-Morgan', 
'OpenBabel-FP3', 'OpenBabel-FP2',
'ChemFP-Substruct-OpenBabel', 
'RDKit-MACCS166'
    )
       
"""fptypes = (
    'RDKit-Pattern', 'OpenEye-Path', 'OpenBabel-MACCS', 'RDKit-Avalon',
'RDKit-AtomPair', 'RDKit-Fingerprint', 'OpenEye-SMARTSScreen',
'OpenBabel-ECFP2', 'RDKit-SECFP', 'RDKit-Torsion',
'OpenBabel-ECFP8', 'ChemFP-Substruct-RDKit', 'RDMACCS-OpenEye',
'OpenBabel-ECFP6', 'RDMACCS-OpenBabel', 'OpenEye-MDLScreen',
'OpenEye-MACCS166', 'RDMACCS-RDKit', 'OpenBabel-FP4',
'OpenEye-Tree', 'RDKit-Morgan', 'ChemFP-Substruct-OpenEye',
'OpenBabel-FP3', 'OpenBabel-FP2', 'OpenBabel-ECFP0',
'ChemFP-Substruct-OpenBabel', 'OpenEye-Circular',
'OpenBabel-ECFP10', 'OpenBabel-ECFP4', 'OpenEye-MoleculeScreen',
'RDKit-MACCS166'
    )"""


"""actius = [mol for mol in pybel.readfile("sdf", sys.argv[1])]
problemes = [mol for mol in pybel.readfile("sdf", sys.argv[2])]



print("\nFITXER ACT:\t "+sys.argv[1]+"\t"+str(len(actius))+" molècules")
print("FITXER DEC:\t "+sys.argv[2]+"\t"+str(len(problemes))+" molècules")"""

def crearLlistaTuple(list_of_molecules, esAct):
	
	
	return [(mol, esAct) for mol in list_of_molecules]
	
def trobarMaxims(llista_actius, llista_totals):  #Mateixa molècula, repetits, ig
		
	
	maxims = [[0 for x in range(4)] for y in range(len(llista_totals))]	
	
	llistaFPA = [T.compute_fingerprint(llista_actius[a][0]) for a in range(len(llista_actius))]
	llistaFPT = [T.compute_fingerprint(llista_totals[b][0]) for b in range(len(llista_totals))]
	
	T=fptype.toolkit
	
	for i in range(len(llista_totals)):
		
		maxims[i][0] = T.create_string(llista_totals[i][0],"smistring")
		maxims[i][1] = 0
		
		for j in range(len(llista_actius)):
			tan=chemfp.bitops.byte_tanimoto(llistaFPT[i],llistaFPA[j])
			if ((tan > maxims[i][1]) and (llista_totals[i][0] is not llista_actius[j][0])): #Comprovar el is not / mateixa molecula
				
				maxims[i][1] = tan
				maxims[i][2] = llista_totals[i][1]
				maxims[i][3] = T.create_string(llista_actius[j][0],"smistring")	
		
	return maxims	

def ordenarTanimotos(llista_maxims):
	
	llista_maxims.sort(key=lambda x: x[2],reverse = True)
	
	return llista_maxims

def calcularEF(factor, llistaTuplesOrdenada, num_actius):
	
	actius_trobats = 0
	ef_llargada = len(llistaTuplesOrdenada)*factor/100
	
	for i in range(ef_llargada):
		if llistaTuplesOrdenada[i][2]:
			 actius_trobats+=1
	
	ef = 100 * actius_trobats / num_actius
	
	return ef
	
def calcularBEDROC(llistaTuplesOrdenada):
	
	llista_scores = [(1-el[2], el[3]) for el in llistaTuplesOrdenada]
	bedroc = Scoring.CalcBEDROC(llista_scores, 1, 20) 

	return bedroc

def eliminar_repetits(mfile):
	
	if(mfile.endswith(".sdf")):
		#mols=[mol for mol in pybel.readfile("sdf", mfile)]
		unique_mols = {mol.write("inchi") : mol for mol in pybel.readfile("sdf", mfile)}
		
	elif(mfile.endswith(".smi")):
			unique_mols = {mol.write("inchi") : mol for mol in pybel.readfile("smi", mfile)}
			
	outputsdf = pybel.Outputfile("sdf", str(mfile[:-4])+"_uniques.sdf", overwrite=True) 
	for mol in unique_mols.itervalues(): 
		outputsdf.write(mol) 
	
	outputsdf.close() 
		
def funcio_general(fingerprint):

	fptype=chemfp.get_fingerprint_type(fingerprint)
	T=fptype.toolkit
	
	mfile=str(sys.argv[1])
	mfile=mfile[:-4]+"_uniques.sdf"
			
	with T.read_molecules(mfile) as reader:
		actives=[T.copy_molecule(mol) for mol in reader]
		
	mfile=str(sys.argv[2])
	mfile=mfile[:-4]+"_uniques.sdf"
	
	with T.read_molecules(mfile) as reader:
		decoys=[T.copy_molecule(mol) for mol in reader]

	llistaActius = crearLlistaTuple(actives,1)
	llistaTotal = crearLlistaTuple(decoys,0)

	llistaTotal = llistaActius + llistaTotal
	
	df_max = pd.DataFrame()
	#df_max["Molècula"] = llistaTotal	
	
	maxims = [[0 for x in range(6)] for y in range(len(llistaTotal))]	
	
	llistaFPA = [fptype.compute_fingerprint(llistaActius[a][0]) for a in range(len(llistaActius))]
	llistaFPT = [fptype.compute_fingerprint(llistaTotal[b][0]) for b in range(len(llistaTotal))]
	
	
	for i in range(len(llistaTotal)):
		
		maxims[i][0] = T.get_id(llistaTotal[i][0])
		maxims[i][1] = T.create_string(llistaTotal[i][0],"smistring")
		maxims[i][2] = 0
		
		for j in range(len(llistaActius)):
			tan=chemfp.bitops.byte_tanimoto(llistaFPT[i],llistaFPA[j])
			id_d = T.get_id(llistaTotal[i][0])
			id_a = T.get_id(llistaActius[j][0])
			if ((tan > maxims[i][2]) and (id_d is not id_a)): #Descartar si s'està comparant una molècula amb ella mateixa
				maxims[i][2] = tan
				maxims[i][3] = llistaTotal[i][1]
				maxims[i][4] = T.get_id(llistaActius[j][0])
				maxims[i][5] = T.create_string(llistaActius[j][0],"smistring")
				
	maxims = ordenarTanimotos(maxims)
	
	df_max = pd.DataFrame(maxims, columns =['Molecule ID','Molecule SMILES','Tanimoto', 'Is Active', 'Closer Active ID', 'Closer Active SMILES'])
		
	df_max.to_csv(r'/home/ori/Ophidian/Results/'+str(fingerprint)+'.csv')
	print(str(fingerprint)+" COMPLETED")
	
	metriques[0] = fingerprint
	metriques[1] = (calcularEF(1,maxims,len(llistaActius)))
	metriques[2] = (calcularEF(10,maxims,len(llistaActius)))
	metriques[3] = sklearn.metrics.roc_auc_score(df_max['Is Active'],df_max['Tanimoto'])
	metriques[4] = calcularBEDROC(maxims)
	
	return(metriques)
	



if __name__ == '__main__':
	
	eliminar_repetits(sys.argv[1])
	eliminar_repetits(sys.argv[2])
	
	metriques = [0 for x in range(5)]
	
	pool = multiprocessing.Pool()
	resultats=pool.map(funcio_general,fptypes)
	
	df_met=pd.DataFrame(resultats, columns =['Fingerprint','EF1%','EF10%','AUC','BEDROC'])
	df_met.to_csv(r'/home/ori/Ophidian/Resultats/metrics.csv')

	print("\n\nMètriques generades:\n")
	print("\nBest fp according to EF1%:\n")
	print(df_met.iloc[df_met['EF1%'].idxmax()])
	print("\n\nBest fp according to EF10%:\n")
	print(df_met.iloc[df_met['EF10%'].idxmax()])
	print("\n\nBest fp according to AUC:\n")
	print(df_met.iloc[df_met['AUC'].idxmax()])
	print("\n\nBest fp according to BEDROC:\n")
	print(df_met.iloc[df_met['BEDROC'].idxmax()])

