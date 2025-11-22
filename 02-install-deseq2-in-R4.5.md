# Start R
To start R Statistical Programming language type R in the command line of the R4.5WS conda environment

```{cmd}
R
```

# Install Packages for Bioconductor
Installing `BiocManager` package from CRAN (Comprehensive R Archival Network) using `install.packages()` function in R base. Please type the following. 

```{R}
install.packages("BiocManager")
```

Install `Biostrings` Bioconductor package using `BiocManager::install()` function.

```{R}
BiocManager::install("Biostrings")
```