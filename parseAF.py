#!/usr/bin/env python3
"""Iterate AlphaFold PDB files in folder and run DSSP, PropKa and S-S distance calculations."""

__author__  = "Patrick Willems"
__email__   = "willems533@gmail.com"           

import json
from urllib.request import urlopen
import urllib
import os
import sys
import numpy
import glob
import re
from Bio.PDB import *

f = open("output.txt", "w")
f.write("Accession\tPosition\tSecStructure\tRSA\tSS_partner\tSS_Distance(A)\tSS_PAE\tpKa\tBuried_perc\n")
for PDB_file in glob.glob("AF/AF*pdb"):
    print('Running now ' + PDB_file)
    accession = re.findall('AF-(.+)-F', PDB_file)[0]
    cys = {"pos": [],"ss": [],"rsa": [], "ss_cys" : [], "ss_dist" : [], "ss_pae": [], "pka": [], "buried": []} #Initialize
	
    #RUN DSSP
    p = PDBParser()
    structure = p.get_structure(accession, PDB_file)
    residues = [r for r in structure.get_residues()]
    model = structure[0]    
    dssp = DSSP(model, PDB_file)
    for i in dssp:
        pos, res, ss, rsa = i[0:4] 
        if res == 'C':
            cys["pos"].append(pos)
            cys["ss"].append(ss)
            cys["rsa"].append(format(rsa,'.4f'))
    
    #Read in AlphaFold PAE JSON
    url = 'https://alphafold.ebi.ac.uk/files/AF-' + accession + '-F1-predicted_aligned_error_v2.json'
    response = urllib.request.urlopen(url)
    PAE = json.loads(response.read())

    #Cys S-S distances
    for cys1 in cys["pos"]:
        cys1 = cys1 - 1
        cys1_coord = residues[cys1]["SG"].get_coord()
        ssFound = False
        for cys2 in cys["pos"]:
            cys2 = cys2 - 1
            if cys1 == cys2:
                 continue
            else:
                cys2_coord = residues[cys2]["SG"].get_coord()
                distance = numpy.linalg.norm(cys1_coord-cys2_coord)
                if distance < 3:
                    ssFound = True
                    cys["ss_cys"].append(cys2 + 1)
                    cys["ss_dist"].append(format(distance,'.4f'))
                    PAE_index = PAE[0]["residue1"].index(cys1 + 1) + cys2
                    cys["ss_pae"].append(PAE[0]["distance"][PAE_index])
                    break
        if ssFound == False:
            cys["ss_cys"].append('NA')
            cys["ss_dist"].append('NA')
            cys["ss_pae"].append('NA')
    
    #Run Propka3
    os.system("python3 -m propka " + PDB_file)
    pka_f = re.findall('(AF-.+.)pdb', PDB_file)[0] + 'pka'
    with open(pka_f) as file:
        for line in file:
            if re.match('^CYS\s*\d+\s+A.+%', line):
                pka = re.findall('^CYS\s*\d+\s+A\s+(\S+)\s+', line)[0]
                buried = re.findall('^CYS\s*\d+\s+A\s+\S+\s+(\d+) %', line)[0]
                cys["pka"].append(str(pka))
                cys["buried"].append(str(buried))
    os.system("rm *pka")

    #Print all output
    for i in range(0, len(cys["pos"])):
            f.write(accession+'\t'+str(cys["pos"][i])+'\t'+cys["ss"][i]+'\t'+str(cys["rsa"][i])+'\t'+str(cys["ss_cys"][i])+'\t'+str(cys["ss_dist"][i])+'\t'+str(cys["ss_pae"][i])+'\t'+cys["pka"][i]+'\t'+cys["buried"][i]+'\n')
f.close()
