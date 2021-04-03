# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------
# This is a sample controller
# this file is released under public domain and you can use without limitations
# -------------------------------------------------------------------------

from  rdkit.ML.Scoring import Scoring
import chemfp
from chemfp import bitops
import pandas as pd
import sklearn
from sklearn import metrics
from openbabel import pybel
import multiprocessing
from multiprocessing import Pool

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

def index():
	record = db.arxius(request.args(1))
	form = SQLFORM(db.arxius, record, deletable=True, upload=URL('download'))
	if form.process().accepted:
		response.flash = T('Files uploaded succesfully')
		redirect(URL('processing'))
	elif form.errors:
		response.flash = T('The form has errors')
	
	return dict(form=form)
	


# ---- action to server uploaded static content (required) ---
@cache.action()
def download():
    """
    allows downloading of uploaded files
    http://..../[app]/default/download/[filename]
    """
    return response.download(request, db)
    
def processing():
		
	row = db(db.arxius).select().first()
	
	(act, astream) = db.arxius.actius.retrieve(row.actius)
	#act = db(db.arxius).select("actius")

	#(dec, dstream) = db.arxius.decoys.retrieve(row.decoys)
	
	actlist=astream.read()
	#declist=dstream.read()
	
	amols=list(pybel.readstring("sdf", actlist))
	#dmols=[mol for mol in pybel.readfile("sdf", declist)]
	
	print len(amols)
	
	return dict(amols=amols)
	
