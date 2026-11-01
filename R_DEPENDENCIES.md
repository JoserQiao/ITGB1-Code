# R dependency for IMvigor210 export

The manuscript analysis used R 4.6.1, Bioconductor 3.23 and DESeq2 1.52.0.
The R step only exports normalized expression and clinical metadata from the
public `cds.RData` object; all statistical models are run by the locked Python
environment.

```r
install.packages("BiocManager")
BiocManager::install(version = "3.23")
BiocManager::install("DESeq2")
stopifnot(packageVersion("DESeq2") == "1.52.0")
```
