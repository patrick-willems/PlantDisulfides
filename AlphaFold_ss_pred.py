#!/usr/bin/env python3
"""Iterate AlphaFold PDB files in folder and run DSSP, PropKa3, perform S-S distance calculations and parse corresponding PAE from the AF predictions."""

__author__  = "Patrick Willems"
__email__   = "willems533@gmail.com"           

import math
import os
import sys
import numpy
import re
import os.path
from Bio.PDB import *

###Pass UniProtKB proteome FASTA from argument to extract IDs
fasta = sys.argv[1]

###PREPARE AND READ CURRENT RESULT FILE###
#All results will be appended to single result file 'results.txt', write header if non-existant.
outfile = 'results.txt'
f = open(outfile, "a")
if os.stat('results.txt').st_size == 0:
    f.write("TaxID\tAccession\tAF-model\tPosition\tpLDDT\tClass\tSecStruct\tRSA\tSS_Cys\tSS_avg_pLDDT\tSS_dist\tCA_dist\tX1\tX2\tX3\tX1\'\tX2\'\tLigand\tProtein_#Cys\tProtein_Length\n")

#Check which accessions have already been performed so the script does not re-run these.
os.system("cat results.txt | cut -f2 | uniq > done_accessions")
done = []
done_acc = open("done_accessions",'r')
for line in done_acc:
    line = line.rstrip()
    done.append(line)
    
###READ FASTA###
#This will iterate over the provided FASTA file and parse the number of cysteines and protein lengths for overall statistic purposes.
#From the UniProt header lines it will store the UniProtKB accessions along with their species taxid from the OX element.
nrCys = {}
lengths = {}
taxids = {}
fasta_f = open(fasta,'r')
acc = ''
for line in fasta_f:
    line = line.rstrip()
    if re.search(r'>',line):
        acc = line.split('|')[1]
        taxid = re.findall(r'OX=(\d+) ',line)[0]
        if acc in done:
          print("Already done",acc)
          continue
        else:
          taxids.update({acc: taxid})
          nrCys.update({acc: 0})
          lengths.update({acc: 0})
    else:
        if acc in done:
          continue
         else:
          nrCys[acc] += line.count('C')
          lengths[acc] += len(line)

###LOOP OVER ACCESSIONS IN PROVIDED FASTA###
for acc in nrCys:
    taxid = taxids[acc]
    
    #There are no cysteines in this protein so skip it.
    if nrCys[acc] == 0:
        f.write(taxid + '\t' + acc + '\tNA' * 16 + '\t' + str(nrCys[acc])+'\t'+str(lengths[acc])+'\n')
    else:
        ###GET ALPHAFOLD2###
        #First of all, attempt to download the AlphaFold2 model
        PDB_file = 'AF-' + acc + '-F1-model_v3.pdb'
        path = 'https://alphafold.ebi.ac.uk/files/AF-' + acc + '-F1-model_v3.pdb'
        os.system('wget -q ' + path + ' -O ' + PDB_file)
        
        #Dealing with non-existing AlphaFold model, so skip
        if os.stat(PDB_file).st_size == 0: 
            f.write(taxid + '\t' + acc + '\tNoModel' * 16+ '\t' + str(nrCys[acc])+'\t'+str(lengths[acc])+'\n')
            os.system("rm -f " + PDB_file)
            continue
        
        ###READ PDB STRUCTURE###
        #There is an AlphaFold2 model available, start reading in with Biopython PDB module
        print('Running now ' + PDB_file)
        p = PDBParser()
        structure = p.get_structure(acc, PDB_file)
        results = {}
        for r in structure.get_residues():
            if re.search('CYS ',str(r)):
                pos = int(re.findall(r'resseq=(\d+) ',str(r))[0])
                results.update({pos: {'pLDDT': '','ss_cys': 'NA', 'ss_dist': 'NA','ss_pLDDT': 'NA', 'ca_dist': 'NA', 'X1': 'NA', 'X2': 'NA', 'X3': 'NA', 'X2_prime': 'NA', 'X1_prime' : 'NA', 'ligand': 'NA', 'class': 'none', 'secStruct': 'NA', 'rsa': 'NA', 'pka': ''}})
        residues = [r for r in structure.get_residues()]
        
        ###pLDDT VALUES###
        #Read in the pLDDT values for all cysteines within the structure. These values are specified within the B-factor fields of AlphaFold2 PDB files.
        PDB_read = open(PDB_file,'r')
        for line in PDB_read:
            if re.search(r'^ATOM.+CYS A',line):
                pos = int(re.findall(r'CYS A\s*(\d+) ',line)[0])
                pLDDT = float(line[61:66])
                fields = line.split()
                results[pos]['pLDDT'] = pLDDT

        ###RUN DSSP FOR RSA AND SEC STRUCT###
        model = structure[0]
        dssp = DSSP(model, PDB_file)
        for i in dssp:
            pos, res, ss, rsa = i[0:4]
            if res == 'C' and pos in results:
                results[pos]['rsa'] = format(rsa,'.4f')
                results[pos]['secStruct'] = ss

        ###EXTRACT DISULFIDES###
        #Calculate all pairwise cysteine S-S distances and dihedral angles in case of disulfides
        for cys1 in results:
            pos = cys1
            cys1 = cys1 - 1
            cys1_SG = residues[cys1]["SG"].get_coord()
            for cys2 in results:
                cys2 = cys2 - 1
                if cys1 == cys2:
                    continue
                else:
                    cys2_SG = residues[cys2]["SG"].get_coord()
                    dist_SG = numpy.linalg.norm(cys1_SG-cys2_SG)
                    if dist_SG < 2.5:
                        results[pos]["class"] = 'disulfide'
                        results[pos]["ss_cys"] = str(cys2 + 1) 
                        results[pos]["ss_dist"] = str(format(dist_SG,'.4f')) 
                        results[pos]["ca_dist"] = str(format(numpy.linalg.norm(residues[cys1]["CA"].get_coord()-residues[cys2]["CA"].get_coord()),'.4f'))
                        results[pos]["ss_pLDDT"] = (results[pos]["pLDDT"] + results[(cys2+1)]["pLDDT"])/2

                        ##ANGLE CALCULATION - DEFINE ALL VECTORS FIRST
                        cys1_S = residues[cys1]["SG"].get_vector()
                        cys1_N = residues[cys1]["N"].get_vector()
                        cys1_CA = residues[cys1]["CA"].get_vector()
                        cys1_CB = residues[cys1]["CB"].get_vector()
                        cys2_S = residues[cys2]["SG"].get_vector()
                        cys2_CB = residues[cys2]["CB"].get_vector()
                        cys2_CA = residues[cys2]["CA"].get_vector()
                        cys2_N = residues[cys2]["N"].get_vector()
                        
                        #CALCULATE THE FIVE DIHEDRAL ANGLES
                        X1 = calc_dihedral(cys1_N, cys1_CA, cys1_CB, cys1_S)*180/math.pi
                        X2 = calc_dihedral(cys1_CA, cys1_CB, cys1_S,cys2_S)*180/math.pi
                        X3 = calc_dihedral(cys1_CB, cys1_S,cys2_S,cys2_CB)*180/math.pi
                        X2p = calc_dihedral(cys1_S,cys2_S,cys2_CB,cys2_CA)*180/math.pi
                        X1p = calc_dihedral(cys2_S,cys2_CB,cys2_CA,cys2_N)*180/math.pi
                        results[pos]["X1"] = format(X1,'.4f')
                        results[pos]["X2"] = format(X2,'.4f')
                        results[pos]["X3"] = format(X3,'.4f')
                        results[pos]["X2_prime"] = format(X2p,'.4f')
                        results[pos]["X1_prime"] = format(X1p,'.4f')
                             
        ###RUN METALLOPROTEOME###
        #This runs the metal search algorithm described by Wehrspan et al. (2022) - see https://github.com/Elcock-Lab/Metalloproteome 
        os.system("./bin/analyze_alphafold2_metal_clusters_for_release.exe " + PDB_file + " ligand_list_FES_ZINC_RMSD_0.5_12_LIGANDS 0.0 8.0 2.0 2.5 0.0 2.5 998 > /dev/null 2>&1")
        output_metal = open("ligand_summary_info_000998.txt",'r')
        for line in output_metal:
            line = line.rstrip()
            if re.search('CYS_00', line):
                metal = re.findall(r'ligand_name (.+).pdb',line)[0]
                cys_lig = re.findall(r'CYS_0+(\d+)',line)
                for cys in cys_lig:
                    results[int(cys)]['ligand'] = metal
                    results[int(cys)]['class'] = 'metal'

        #Print all output and clean-up pdb
        for pos in results:
            f.write(taxid+'\t'+acc+'\t'+PDB_file+'\t'+str(pos)+'\t'+str(results[pos]["pLDDT"])+'\t'+results[pos]['class']+'\t'+results[pos]['secStruct']+'\t'+results[pos]['rsa']+'\t'+str(results[pos]['ss_cys'])+'\t'+str(results[pos]["ss_pLDDT"])+'\t'+str(results[pos]['ss_dist'])+'\t'+str(results[pos]['ca_dist'])+'\t'+str(results[pos]['X1'])+'\t'+str(results[pos]['X2'])+'\t'+str(results[pos]['X3'])+'\t'+str(results[pos]['X2_prime'])+'\t'+str(results[pos]['X1_prime'])+'\t'+results[pos]['ligand']+'\t'+str(nrCys[acc])+'\t'+str(lengths[acc])+'\n')
        os.system("rm -f " + PDB_file + " *0998.txt *.pka")
f.close()
