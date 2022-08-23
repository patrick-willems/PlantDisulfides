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
- For each cysteine:
  - The residue the relative solvent accessibility (RSA) and secondary structure will be predicted using DSSP (available Bio.PDB)
  - The pLDDT value will be extracted from the AlphaFold2 PDB file.
  - Possible disulfides will be defined if the 'SG' atom is within 2.5 A of any other cysteine. 
- For all disulfides:
  - The S-S atom distance and CA-CA atom distance is calculated.
  - Five dihedral angles are calculated (X1, X2, X3, X1' and X2') using Bio.PDB calc_dihedral().
  - The average pLDDT is calculated as (pLDDT SS_Cys1 + pLDDT SS_Cys2) / 2 
- For all proteins:
  - Metal ligand binding was predicted using the Metalloproteome algorithm ([xx](https://github.com/Elcock-Lab/Metalloproteome)), required files were copied within this repository.
  - The length and number of cysteine of the protein is provided (derived from FASTA).
  - Optionally a pKa value can be predicted using propka3.0 (requires local install) - however this takes relatively long making whole proteome predictions too long, hence was skipped here.

