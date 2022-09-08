# PDC: Protein Data Compressor #

## Introduction ##
With recent development of high accuracy protein structure predictors, more and more predicted protein structure models have been deposited to public databases such as the [AlphaFold DB](https://alphafold.ebi.ac.uk/). This leads to hugh hard disk comsumptions. For example, the full AlphaFold DB release in year 2022 has 23 TB of data, which is only expected to increase significantly in the near future. To address this issue, the PDC package aims to convert PDB and mmCIF format protein structure models to and from the highly compressed .pdc format, both losslessly and lossily.

## Approach and Implementation ##
PDC decrease the size of protein coordinate files in PDB or mmCIF format through the following three approaches:
1. Removal of repetitive information among different atoms, such as the chain ID and residue index.
2. Use ``float`` and ``char`` instead of ``string`` to store coordinates and B-factors, respectively. 
3. Under lossy mode, the coordinates of some atoms will be discarded, including:
  * backbone atoms: CB, O, OXT (recovered from N, CA, C during decompression)
  * HIS sidechain atoms: CG, ND1, NE2 (recovered from CB, CE1, CD2 during decompression)
  * PHE sidechain atoms: CG, CD1, CD2, CZ (recovered from CB, CE1, CE2 during decompression)
  * TYR sidechain atoms: CG, CD1, CD2, CZ, OH (recovered from CB, CE1, CE2 during decompression)
  * TRP sidechain atoms: CG, CD1, CD2, NE1, CE2, CZ3, CH2 (recovered from CB, CE3, CZ2 during decompression)

## Limitations ##
PDC is specifically designed for protein models in the AlphaFold database. It is not able to convert all information of a PDB or mmCIF file, especially those from the [PDB database](https://www.rcsb.org/). In particular,
1. Information for small molecule ligands and non-standard residues (ATOM and CONECT) are ignored.
2. A file cannot be parsed if there are missing atoms in some residues.
3. Only MODEL 1 of in a multi-model structure will be converted.
4. For atoms with alternative locations, only atoms with alternative locations ' ' or 'A' will be considered.
5. Hydrogens are ignored.
