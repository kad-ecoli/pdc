# PDC: Protein Data Compressor #

## Introduction ##
With recent development of high accuracy protein structure predictors, more and more predicted protein structure models have been deposited to public databases such as the [AlphaFold DB](https://alphafold.ebi.ac.uk/). This leads to hugh hard disk comsumptions. For example, the full AlphaFold DB release in year 2022 has 23 TB of data, which is only expected to increase significantly in the near future. To address this issue, the PDC package aims to convert PDB and mmCIF format protein structure models to and from the highly compressed .pdc format, both losslessly and lossily.
