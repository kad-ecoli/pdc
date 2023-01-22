# PDC: Protein Data Compressor #

## Introduction ##
Recent development of high accuracy protein structure predictors result in more and more predicted protein structure models being deposited to public databases such as the [AlphaFold DB](https://alphafold.ebi.ac.uk/). These large datasets for protein structures leads to hugh hard disk comsumptions. For example, the full AlphaFold DB release in year 2022 has 23 TB of data, which is expected to continuously increase. To address the data storage issue, the PDC package aims to convert full atomic PDB and mmCIF format protein structure models to and from the highly compressed .pdc format, which is specifically designed for AlphaFold predicted protein structures.

## Installation ##
```bash
make
```
PDC does not run natively on Windows. However, it can be run on Windows Subsystem for Linux.

## Usage ##
Lossless compression:
```bash
pdc AF-P11532-F2-model_v3.pdb.gz AF-P11532-F2-model_v3.pdc.gz
```
Lossy compression:
```bash
pdc AF-P11532-F2-model_v3.pdb.gz AF-P11532-F2-model_v3.pdc.gz -l=2
```
Lossless compression (CA atoms only)
```bash
pdc AF-P11532-F2-model_v3.pdb.gz AF-P11532-F2-model_v3.pdc.gz -l=3
```
Lossy compression (CA atoms only)
```bash
pdc AF-P11532-F2-model_v3.pdb.gz AF-P11532-F2-model_v3.pdc.gz -l=4
```
Uncompress:
```bash
pdd AF-P11532-F2-model_v3.pdc.gz AF-P11532-F2-model_v3.pdb.gz
```

## Approach and Implementation ##
PDC decrease the size of protein coordinate files in PDB or mmCIF format through the following three approaches:
1. Removal of repetitive information among different atoms, such as the chain ID and residue index.
2. Use [int](https://en.cppreference.com/w/cpp/types/integer) and ``char`` instead of ``string`` to store coordinates and B-factors.
   Specifically, since xyz and bfactor can be expressed as %8.3f and %6.2f, they are in the range of -999.999 to 9999.999 and -99.99 to 999.99, respectively. This means that they can be expressed as integers in the range of 0 to 10999998 and 0 to 109998, respectively, both of which can be stored by unint32.
3. Delta encoding: store the difference in coordinate/bfactor from the previous value rather than the actual value, which is can be stored by int16 or int8.
4. Under lossy compression mode, store the torsion angles rather than the coordinates.

## Limitations ##
PDC is specifically designed for protein models in the AlphaFold database. It is not able to convert all information of a PDB or mmCIF file, especially those from the [PDB database](https://www.rcsb.org/). In particular,
1. Information for small molecule ligands and non-standard residues (ATOM and CONECT) are ignored.
2. A file cannot be parsed if there are missing atoms in some residues.
3. Only MODEL 1 of in a multi-model structure will be converted.
4. For atoms with alternative locations, only atoms with alternative locations ' ' or 'A' will be considered.
5. Hydrogens are ignored.

## Benchmark ##
PDC, [MMTF](https://mmtf.rcsb.org/), [PIC](https://github.com/lukestaniscia/PIC) and [BinaryCIF](https://github.com/molstar/BinaryCIF) are applied to the [E coli](https://ftp.ebi.ac.uk/pub/databases/alphafold/v3/UP000000625_83333_ECOLI_v3.tar) proteome of AlphaFold DB. The file sizes after gzip compression are shown below.
| File format | File size (MB) | Lossless/Lossy |
| :--:        | :--:           | :--:           |
| CIF         | 273            | Lossless       |
| PDB         | 196            | Lossless       |
| PIC         | 163            | Lossy          |
| BinaryCIF   | 143            | Lossless       |
| MMTF        | 76             | Lossless       |
| PDC         | 67             | Lossless       |
| PDC         | 26             | Lossy          |
