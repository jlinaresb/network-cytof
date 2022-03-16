# Network analysis of cell line XXXX based on several perturbations

## Current stage

![pb](https://progress-bar.dev/50/?title=Network)
![pb](https://progress-bar.dev/0/?title=CyTOF)
![pb](https://progress-bar.dev/0/?title=Training)


## Getting Started

There are several steps to run this analyisis:

1. Identify proteins of interest
2. Identify drug targets
3. Create the network
4. CyTOF data
5. Run [CellNOpt](https://saezlab.github.io/CellNOptR/)


### Prerequisites:

You must install the following packages in R:


```{r}
install.packages(c("ggplot2", "dplyr", "ggraph", "igraph"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("OmnipathR", "CellNOptR"))
```

### How to run

Clone this repository:

```{bash}
git clone https://github.com/jlinaresb/network-cytof.git
```

