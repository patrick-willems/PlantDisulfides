#!/usr/bin/env python3
"""Iterate AlphaFold PDB files of given UniProtKB FASTA and run DSSP, PropKa3, perform S-S distance calculations and parse corresponding PAE from the AF predictions."""

__author__  = "Patrick Willems"
__email__   = "willems533@gmail.com"           

import math
import os
import sys
import numpy
import glob
import re
import os.path
from Bio.PDB import *
from Bio.PDB.ResidueDepth import min_dist
from Bio.PDB.ResidueDepth import residue_depth

#Will append to single result file, write header if non-existant.
outfile = 'results_ext.txt'
f = open(outfile, "a")
if os.stat('results_ext.txt').st_size == 0:
    f.write("TaxID\tAccession\tAF-model\tPosition\tpLDDT\tClass\tSecStruct\tRSA\tSS_Cys\tSS_avg_pLDDT\tSS_dist\tCA_dist\tX1\tX2\tX3\tX1\'\tX2\'\tLigand\tDistance_Surface_Cys\tDistance_Surface_SG\tpKa\tProtein_#Cys\tProtein_Length\n")

#Check which accessions have already been checked.
os.system("cat results_ext.txt | cut -f2 | uniq > done_accessions_ext")
done = []
done_acc = open("done_accessions_ext",'r')
for line in done_acc:
    line = line.rstrip()
    done.append(line)

#Pass UniProtKB proteome FASTA from argument to extract IDs
fasta = sys.argv[1]
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

#Perform the real work for all the accessions within the FASTA
for acc in nrCys:
    taxid = taxids[acc]
    if nrCys[acc] == 0:
        f.write(taxid + '\t' + acc + '\tNA' * 19 + '\t' + str(nrCys[acc])+'\t'+str(lengths[acc])+'\n')
    else:
        PDB_file = 'AF-' + acc + '-F1-model_v3.pdb'
        path = 'https://alphafold.ebi.ac.uk/files/AF-' + acc + '-F1-model_v3.pdb'
        os.system('wget -q ' + path + ' -O ' + PDB_file)
        if os.stat(PDB_file).st_size == 0:  #Dealing with non-existing AlphaFold model
            f.write(taxid + '\t' + acc + '\tNoModel' * 19 + '\t' + str(nrCys[acc])+'\t'+str(lengths[acc])+'\n')
            os.system("rm -f " + PDB_file)
            continue
        print('Running now ' + PDB_file)
        
        #Read structure with Biopython
        p = PDBParser()
        structure = p.get_structure(acc, PDB_file)
        results = {}
        for r in structure.get_residues():
            if re.search('CYS ',str(r)):
                pos = int(re.findall(r'resseq=(\d+) ',str(r))[0])
                results.update({pos: {'pLDDT': '','ss_cys': 'NA','distSurface': '999', 'distSurface_SG': '999', 'ss_dist': 'NA','ss_pLDDT': 'NA', 'ca_dist': 'NA', 'X1': 'NA', 'X2': 'NA', 'X3': 'NA', 'X2_prime': 'NA', 'X1_prime' : 'NA', 'ligand': 'NA', 'class': 'none', 'secStruct': 'NA', 'rsa': 'NA', 'pka': ''}})
        residues = [r for r in structure.get_residues()]
        
        ###GET THE pLDDT CONFIDENCE SCORES OF THE PREDICTION###
        PDB_read = open(PDB_file,'r')
        for line in PDB_read:
            if re.search(r'^ATOM.+CYS A',line):
                pos = int(re.findall(r'CYS A\s*(\d+) ',line)[0])
                pLDDT = float(line[61:66])
                fields = line.split()
                results[pos]['pLDDT'] = pLDDT

        ###RUN DSSP FOR RSA AND SEC STRUCT###
        model = structure[0]
        chain = model['A']
        dssp = DSSP(model, PDB_file)
        for i in dssp:
            pos, res, ss, rsa = i[0:4]
            if res == 'C' and pos in results:
                results[pos]['rsa'] = format(rsa,'.4f')
                results[pos]['secStruct'] = ss
        
        ###Get surface with ResidueDepth
        surface = get_surface(model,MSMS='/home/ubuntu/data/ss_pred/msms/msms.x86_64Linux2.2.6.1')

        ###EXTRACT DISULFIDES###
        #Calculate all pairwise cysteine S-S distances, dihedral angles and determine configuration
        for cys1 in results:
            pos = cys1
            cys1 = cys1 - 1
            cys1_SG = residues[cys1]["SG"].get_coord()
            results[pos]["distSurface_SG"] = min_dist(cys1_SG, surface)
            results[pos]["distSurface"] = residue_depth(residues[cys1],surface)
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
                        
                        #CALCULATE THE FIVE ANGLES
                        conf = ''
                        X1 = calc_dihedral(cys1_N, cys1_CA, cys1_CB, cys1_S)*180/math.pi
                        conf += '+' if X1 > 0 else '-'
                        X2 = calc_dihedral(cys1_CA, cys1_CB, cys1_S,cys2_S)*180/math.pi
                        conf += '+' if X2 > 0 else '-'
                        X3 = calc_dihedral(cys1_CB, cys1_S,cys2_S,cys2_CB)*180/math.pi
                        conf += '+' if X3 > 0 else '-'
                        X2p = calc_dihedral(cys1_S,cys2_S,cys2_CB,cys2_CA)*180/math.pi
                        conf += '+' if X2p > 0 else '-'
                        X1p = calc_dihedral(cys2_S,cys2_CB,cys2_CA,cys2_N)*180/math.pi
                        conf += '+' if X1p > 0 else '-'
                        results[pos]["X1"] = format(X1,'.4f')
                        results[pos]["X2"] = format(X2,'.4f')
                        results[pos]["X3"] = format(X3,'.4f')
                        results[pos]["X2_prime"] = format(X2p,'.4f')
                        results[pos]["X1_prime"] = format(X1p,'.4f')
                        
        ###RUN METALLOPROTEOME###
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

        ###RUN PROPKA3###
        os.system("propka -s " + PDB_file + ' > /dev/null 2>&1')
        pka_f = re.findall('(AF-.+.)pdb', PDB_file)[0] + 'pka'
        with open(pka_f) as file:
            for line in file:
                if re.match('^CYS\s*\d+\s+A.+%', line):
                    pos = re.findall('^CYS\s*(\d+)\s+A\s+\S+\s+', line)[0]
                    pka = re.findall('^CYS\s*\d+\s+A\s+(\S+)\s+', line)[0]
                    results[int(pos)]["pka"] = str(pka)

        #Print all output and clean-up pdb
        for pos in results:
            f.write(taxid+'\t'+acc+'\t'+PDB_file+'\t'+str(pos)+'\t'+str(results[pos]["pLDDT"])+'\t'+results[pos]['class']+'\t'+results[pos]['secStruct']+'\t'+results[pos]['rsa']+'\t'+str(results[pos]['ss_cys'])+'\t'+str(results[pos]["ss_pLDDT"])+'\t'+str(results[pos]['ss_dist'])+'\t'+str(results[pos]['ca_dist'])+'\t'+str(results[pos]['X1'])+'\t'+str(results[pos]['X2'])+'\t'+str(results[pos]['X3'])+'\t'+str(results[pos]['X2_prime'])+'\t'+str(results[pos]['X1_prime'])+'\t'+results[pos]['conf']+'\t'+results[pos]['ligand']+'\t'+str(results[pos]['distSurface'])+'\t'+str(results[pos]['distSurface_SG'])+'\t'+str(results[pos]['pka'])+'\t'+str(nrCys[acc])+'\t'+str(lengths[acc])+'\n')
        os.system("rm -f " + PDB_file + " *0998.txt *.pka")
f.close()
