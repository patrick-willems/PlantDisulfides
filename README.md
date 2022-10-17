# Plant Disulfides
Parsing disulfide bonds and metal binding sites from AlphaFold2 structure predictions. More information regarding the script is found below. The functional annotation for cysteine residues can be downloaded separately for all fifteen plant species in the DATA folder.

## Installation and requirements
Python3 was used to run this scripts, with following required packages:
- numpy
- math
- Bio.PDB module

For running the metal ligand predictions, please see requirements of the metal ligand searching algorithm ([here](https://github.com/Elcock-Lab/Metalloproteome)).

## Script workflow
The scripts takes a UniProtKB FASTA file as input, which should be placed in the FASTA subfolder. This is the sole required argument of the script, e.g.:

<code>python3 Alphafold_ss_pred.py FASTA/UP000006548_3702.fasta</code>

- All FASTA entries will be parsed and AlphaFold2 models will be downloaded from https://alphafold.ebi.ac.uk/.
- For each cysteine:
  - The residue the relative solvent accessibility (RSA) and secondary structure will be predicted using DSSP (available via Bio.PDB)
  - The pLDDT value will be extracted from the AlphaFold2 PDB file (stored in the B-factor fields).
  - Disulfides will be assigned if the sulfur ('SG') atom is within 2.5 Å distance of another cysteine sulfur atom. 
- For all disulfides:
  - The SG-SG and Cα-Cα atom distance is calculated.
  - Five dihedral angles are calculated (χ1, χ2, χ3, χ1' and χ2') using the Bio.PDB calc_dihedral() function.
  - The average pLDDT is calculated as: (pLDDT Cys1 + pLDDT Cys2) / 2 
- For all proteins:
  - Metal ligand binding was predicted using a published metal ligand searching algorithm ([github link](https://github.com/Elcock-Lab/Metalloproteome)), required files were copied within this repository (see note below).
  - The TaxID, length and number of cysteines of the protein is provided (derived from FASTA).
  - The downloaded AlphaFold2 PDB filename is provided.

Important: make sure to specify correct file paths for the metal ligand searching algorithm (here /bin and /src folders) and place the PDB templates of the ligands within the main folder.

## Output

All predictions are stored in a tab-delimited file 'results.txt' with columns providing all information calculated above for each protein cysteine individually (rows).
If running multiple FASTA files after each other, novel results/rows will be appended to this output file. UniProtKB protein accessions that already have been processed, and thus are within the current 'results.txt' file, will not be re-processed.

## Extended form of the script:

Alphafold_ss_pred_extended.py requires the same FASTA input but will include additional structural predictions such as pKa value by PROPKA3 (Olsson et al., 2011), and average cysteine residue depth as well cysteine sulfur atom residue depth calculated by MSMS (Sanner et al., 1996). This requires the correct installation of these algorithms.
