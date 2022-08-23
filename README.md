# ssbond_parse_plants
Parsing disulfide bonds from AlphaFold2 predictions.

## Installation and requirements
Python3 was used to run this scripts, with following required packages:
- numpy
- math
- Bio.PDB module

## Script workflow
The scripts takes a UniProtKB FASTA file as input, which should be placed in the FASTA subfolder. This is the sole required argument of the script, e.g.:

<code>python3 Alphafold_ss_pred.py FASTA/test.fasta</code>

- All FASTA entries will be parsed and AlphaFold2 models will be downloaded from https://alphafold.ebi.ac.uk/.
- For each cysteine residue the relative solvent accessibility (RSA) and secondary structure will be predicted using DSSP (available Bio.PDB)
- For each cysteine the pLDDT value will be extracted from the AlphaFold2 PDB file.
- For each cysteine possible disulfides will be defined if the 'SG' atom is within 2.5 A of any other cysteine. 
- For all disulfides, the S-S atom distance and CA-CA atom distance is stored. In addition, all five dihedral angles are calculated (X1, X2, X3, X1' and X2') using Bio.PDB.
- Metal ligand binding was predicted using the Metalloproteome algorithm (xx), required files were copied within this repository.
- Optionally, a pKa value can be predicted using propka3.0 (requires local install) - however this takes relatively long making whole proteome predictions too long, hence was skipped here.

