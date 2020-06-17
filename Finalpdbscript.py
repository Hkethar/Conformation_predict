  #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 06:48:17 2020

@author: harini
"""

from Bio.PDB import *
import Bio
import os
import pandas as pd
import math
import csv
from glob import glob
import numpy as np
pdb_data=glob('1a0s.pdb')
#pdbl = PDBList()
#pdbl.retrieve_pdb_file('/PDB data','*.pdb')
parser=PDBParser(PERMISSIVE=True)
data_out=pd.DataFrame()
def degrees(radians):
    if radians:
        angle = radians*180 / math.pi
        while angle > 180:
            angle = 360 - angle
        while angle < -180:
            angle = 360 + angle
        return angle
for filename in pdb_data:
    base=os.path.basename(filename)
    structure_id=os.path.splitext(base)[0]
    structure = parser.get_structure(structure_id, filename)
    phi_psi_degrees = []
    for model in structure:
        for chain in model:
            for residue in chain:
                polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
                for poly_index, poly in enumerate(polypeptides) :
                    #print("Model %s Chain %s" % (str(model.id), str(chain.id)),)
                    #print("(part %i of %i)" % (poly_index+1, len(polypeptides)),)
                    #print("length %i" % (len(poly)),)
                    #print("from %s%i" % (poly[0].resname, poly[0].id[1]),)
                    #print("to %s%i" % (poly[-1].resname, poly[-1].id[1]))
                    degs=poly.get_phi_psi_list()
                    phi_psi_degrees.append(degs)
                    #print(phi_psi_degrees) 
                    print(degs)
                    for res_index, residue in enumerate(poly) :
                        res_name = "%s%i" % (residue.resname, residue.id[1])
                        phi, psi = phi_psi_degrees[res_index]
                        #print(res_name,)
                        print((degrees(phi), degrees(psi)))
df_res=pd.DataFrame(data=phi_psi_degrees,columns=[structure_id])
data_out=pd.concat([data_out,df_res],axis=1)
    #print(phi_psi_degrees)                         
                        #print(data_out)
#outut as csv file, eliminate indexing
data_out.to_csv('res.csv',index=False) 
