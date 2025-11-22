# Create Conda Environment with Miniconda
I am creating a conda environment for DESeq2 dta analysis with R and Bioconductor in my windows systems using windows Subsystem Linux (WSL)

## Conda in WSL terminal
- Start the WSL in Windows.
- Type the following to create a conda environmennt for R4.5 names as R4.5WS

```{cmd}
conda create --name R4.5WS
```
- Activate the R4.5WS conda environment

```{cmd}
conda activate R4.5WS
```
- Install R4.5 in the R4.5WS conda environment

```{cmd}
conda install conda-forge::r-base
```
`conda-forge` is the channel from anaconda.org. [conda-forge/r-base][https://anaconda.org/conda-forge/r-base]



```{cmd}

```

[def]: https://anaconda.org/conda-forge/r-base