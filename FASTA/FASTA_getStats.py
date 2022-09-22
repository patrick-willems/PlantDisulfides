#!/usr/bin/env python3
"""Iterate AlphaFold PDB files in folder and run DSSP, PropKa3, perform S-S distance calculations and parse corresponding PAE from the AF predictions."""

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

#Pass UniProtKB proteome FASTA from argument to extract IDs
fasta = sys.argv[1]
fasta_f = open(fasta,'r')
taxids = {}
motifs = {}
seq = '';
for line in fasta_f:
    line = line.rstrip()
    if re.search(r'>',line):
        if re.match('\w',seq):
            motifs[taxid]['CC'] += len(re.findall('CC',seq))
            motifs[taxid]['CXC'] += len(re.findall(r'C\wC',seq))
            motifs[taxid]['CXXC'] += len(re.findall(r'C\w{2}C',seq))
            motifs[taxid]['CXXXC'] += len(re.findall(r'C\w{3}C',seq))
            motifs[taxid]['CXXXXC'] += len(re.findall(r'C\w{4}C',seq))
            motifs[taxid]['CXXXXXC'] += len(re.findall(r'C\w{5}C',seq))
            seq = ''
        acc = line.split('|')[1]
        taxid = re.findall(r'OX=(\d+) ',line)[0]
        species = re.findall(r'OS=(.+) OX=',line)[0]
        if taxid in taxids:
            taxids[taxid]['prots'] += 1
        else:
            taxids.update({taxid: {'total': 0, 'species': species, 'C': 0, 'prots': 1}})
            motifs.update({taxid: {'CC': 0, 'CXC': 0, 'CXXC': 0, 'CXXXC': 0, 'CXXXXC': 0, 'CXXXXXC': 0}})
    else:
        seq += line
        taxids[taxid]['total'] += len(line)
        taxids[taxid]['C'] += line.count('C')

f = open('stats.txt', "w")
f.write("TaxID\tSpecies\tProteins\tAA_total\tC\tprop_C\n")
for taxid in taxids:
    prop = str(format(taxids[taxid]['C']/taxids[taxid]['total'],'.6f'))
    f.write(taxid+'\t'+taxids[taxid]['species']+'\t'+str(taxids[taxid]['prots'])+'\t'+str(taxids[taxid]['total'])+'\t'+str(taxids[taxid]['C'])+'\t'+prop+'\n')
f.write("\nTaxID\tSpecies\tMotif\tCount\tNormalized\n")
for taxid in taxids:
    for motif in motifs[taxid]:
        prop = str(format((motifs[taxid][motif]*100)/taxids[taxid]['C'],'.6f'))
        f.write(taxid+'\t'+taxids[taxid]['species']+'\t'+motif+'\t'+str(motifs[taxid][motif])+'\t'+prop+'\n')
f.close()
