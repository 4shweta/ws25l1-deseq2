# Activate the conda environment

```{cmd}
conda activate R4.5-ok
```

# Start R in R4.5-ok
```{CMD}
R
```
# Inside R

```{R}
library("DESeq2")
library("clusterProfiler")
```

# DESeq2 tutorial
```{R}
library("tximport")
library("readr")
library("tximportData")
# locating the directory or folder name extdata
dir <- system.file("extdata", package="tximportData")
list.files(path=dir)
# Metadata
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)

samples$condition <- factor(rep(c("A","B"),each=3))

rownames(samples) <- samples$run
#rownames(samples) = samples$run
#samples$run -> rownames(samples)
#a=2

samples[,c("pop","center","run","condition")]
```

# Loadinng salmon quantified files (quant)

```{R}
files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
names(files) <- samples$run

# transcript to gene-map
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))

# loding salmon quant files
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
```
# DESeqDataSet object for Txi Object
```{R}
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(
    txi, colData = samples, design = ~ condition)
```

# Create coldata

```{R}
coldata <- samples
coldata$files <- files
coldata$names <- coldata$run
```

# Create DESeq2 object 

```{R}
library("tximeta")
se <- tximeta(coldata)
ddsTxi_alt <- DESeqDataSet(se, design = ~ condition)
```

# remove lowly expressed genes, keep higly expressed genes
```{R}
dds <- ddsTxi_alt
my_rule<- 3
keep <- rowSums(counts(dds) >= 10) >= my_rule
# To find number of true and columns
table(keep)

dds_keep <- dds[keep,]
dds_de <- DESeq(dds_keep)
res_de <- results(dds_de)

# S4 object to data frame
res_de_df <- data.frames(res_de)




# https://github.com/BioSenior/ggVolcano

```                                  

