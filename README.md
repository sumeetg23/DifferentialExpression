# DifferentialExpression

R script to run a differential expression analysis using DESeq2, generate some plots and do PCA analysis.

R Packages required are:
DESeq2
pheapmap
RColorBrewer
ggplot2
plyr

Input files:
1. Counts table with each row representing gene/unique indentifier
2. columnindentity.txt - This file contains information indentifying control and experimental groups with replicates. Replicates get the same name.
3. comparisons.txt - This file contains information about which comparisons to make based on the experimental names specified in column 2 for columnidentity.txt

Ouput files:
- Significantly differentially expressed genes/id's.
- boxplots of pre and post normalized data
- Correlation heatmap
- PCA analysis plot
- Histogram showing number of counts below 2 for each column in the counts table.
