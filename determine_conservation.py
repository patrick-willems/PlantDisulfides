#!/usr/bin/env python3
__author__  = "Patrick Willems"
__email__   = "willems533@gmail.com"

import json
import subprocess
import os
import sys
import numpy
import glob
import re
from pycanal import Canal


#Read Arabidopsis UniProtKB to plaza identifier
uni2acc = {}
print("Reading Arabidopsis UniProtKB to PLAZA ID")
xref_f = open("plaza/id_conversion.ath.csv", 'r')
for line in xref_f:
    line = line.rstrip()
    if re.search(r'uniprot', line):
        acc, idtype, uni = re.split(r'\t+', line)
        uni2acc.update({uni: acc})

#Read accession to plaza orthology gene family
print("Reading PLAZA ID to GeneFamily")
selected = ["ath", "ptr", "zma", "bna", "bvu", "cre", "ghi", "gma", "han", "mdo", "osa", "ppa", "sly", "smo", "vvi"]
acc2gf = {}
xref_f = open("plaza/genefamily_data.ORTHOFAM.csv", 'r')
for line in xref_f:
    if re.search(r'^\#', line): continue
    line = line.rstrip()
    gf, species, acc = re.split(r'\t', line)
    if species in selected:
        if species == 'ath':
            acc2gf.update({acc: gf})
        else:
            if gf in acc2gf.keys():
                acc2gf[gf].append(acc)
            else:
                acc2gf[gf] = [acc]

#Read the protein sequences of 15 plant proteomes
print("Reading protein sequences from 15 plant species")
acc2seq = {}
acc2species = {}
accession = 'init'
for fasta_f in glob.glob("plaza/proteome*fasta"):
    species = re.findall('proteome.selected_transcript.(\w+).fasta', fasta_f)[0]
    print("  " + species + "..")
    fa = open(fasta_f, 'r')
    for line in fa:
        line = line.rstrip()
        if re.search(r'^>', line):
            protein = re.findall('\s(\S+)$', line)[0]
            acc2species.update({protein: species})
        else:
            line = line.replace('*','')
            acc2seq.update({protein: line})

f = open("output_cons.txt", "w")
f.write("Accession\tOrthofam\tNrSpecies\tNrProteins\tPosition\tScoreEntropy\tFreq_Cys\n")
for PDB_file in glob.glob("AF/AF*pdb"):
    accession = re.findall('AF-(.+)-F', PDB_file)[0]
    print('Running now ' + accession + ' - ' + PDB_file)

    #Check if a PLAZA id (here 'acc') is crossreferenced with UniProtKB and generate MSA
    species = {}
    if accession in uni2acc.keys():
        acc = uni2acc[accession]
        if acc in acc2gf.keys():
            gf = acc2gf[acc]
            if gf in acc2gf:
                nrProteins = len(acc2gf[gf])
                out = open('in.fasta', 'w')
                out.write(">" + accession + "\n" + acc2seq[acc] + "\n")
                for protein in acc2gf[gf]:
                    species.update({acc2species[protein]: True})
                    out.write(">" + protein + "\n" + acc2seq[protein] + "\n")
                out.close
                nrSpecies = str(len(species.keys()))
            else:
                f.write(accession+'\t'+gf+'\t0\t0\tNA\tNA\tNA\n')
                continue
        else:
            f.write(accession+'\t0\t0\tNA\tNA\tNA\n')
            continue
    else:
        f.write(accession+'\t0\t0\tNA\tNA\tNA\n')
        continue

    os.system("clustalo --threads=16 -i in.fasta > msa.fasta")
    if os.stat('msa.fasta').st_size == 0: os.system("clustalo --threads=16 -i in.fasta > msa.fasta")
    if os.stat('msa.fasta').st_size == 0:
        f.write(accession+'\t0\t0\tNA\tNA\tNA\n')
        continue

    canal = Canal(fastafile='msa.fasta',ref=0,startcount=1,verbose=True)
    cons_scores = canal.analysis(include=None, method='relative')
    freq = canal.calcFrequencies()
    cysLoc = re.finditer(r"C", canal.reference_sequence)
    indices = [m.start(0) for m in cysLoc]
    for pos in indices:
        pos += 1
        f.write(accession+'\t'+str(gf)+'\t'+str(nrSpecies)+'\t'+str(nrProteins)+'\t'+str(pos)+'\t'+str(cons_scores.get_value(pos, 'relative'))+'\t'+str(freq[0].get_value("C",pos))+'\n')
    #os.system("rm *fasta")
f.close
