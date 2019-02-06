#! /usr/bin/Rscript

suppressMessages(library(optparse))

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="dataset file name", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default=NULL, help="Path to the input dataset file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

message("\nLoading DESeq2 and other required R packages ....")
suppressMessages(library(DESeq2))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))

message("Input File is ", opt$file)
input.filename <- opt$file

message("Setting working directory to ", opt$dir)
setwd(opt$dir)

countsTable = read.delim(input.filename, row.names=1)
groups = c(colnames(countsTable))
colnames(countsTable) <- groups
countsTable <- countsTable[rowSums(countsTable) > 0,]

tempfile <- read.delim("./columnidentity.txt")
tempfile <- subset(tempfile, tempfile[,2]!="ignore")
tempfile$columnnumber <- tempfile$columnnumber - 1
tempfile <- ddply(tempfile,.(Sample),summarise, countsofsample=length(columnnumber),columnnumber = paste(unique(columnnumber),collapse = ','))
tempfile <- subset(tempfile, tempfile[,2]>1)
rownames(tempfile) <- tempfile[,1]
tempfile[,1] <- NULL

comparisons <- read.delim("./comparisons.txt", header=FALSE)
comparisons$setting <- paste(comparisons[,1],"vs",comparisons[,2], sep="")
comparisons$expcheck <- paste(comparisons[,1] %in% row.names(tempfile), sep="")
comparisons$controlcheck <- paste(comparisons[,2] %in% row.names(tempfile), sep="")
if (length(unique(c(comparisons[,1] %in% row.names(tempfile),comparisons[,2] %in% row.names(tempfile)))) == 2) {
  message("names between comparison and columnidentity files do not match.. exiting")
}
comparisons$listofcol <- paste(tempfile[as.matrix(comparisons[,1]),"columnnumber"], tempfile[as.matrix(comparisons[,2]),"columnnumber"], sep=",")
comparisons$inputstr <- paste(comparisons$setting, tempfile[as.matrix(comparisons[,1]),"columnnumber"], tempfile[as.matrix(comparisons[,2]),"columnnumber"], sep=",")
analysisop <- c(comparisons$inputstr)

genelist <- array()

mincount <- 2
padjthreshold <- 0.01
fcthreshold <- 1

comparisons

for (i in 1:length(analysisop)){

  ops <- unlist(strsplit(analysisop[i],","))
  counts <- countsTable[,as.numeric(as.matrix(unlist(strsplit(analysisop[i],","))[-1]))]

  BelowThreshold<-colSums(counts <= mincount)
  AboveThreshold<-colSums(counts > mincount)
  geneabovethreshold <- t(as.matrix(data.frame(AboveThreshold,BelowThreshold)))

  counts <- counts[rowSums(counts) > mincount,]

  condition <- factor(c(rep(strsplit(ops[1], "vs")[[1]][1],tempfile[as.matrix(strsplit(ops[1], "vs")[[1]][1]),"countsofsample"]),rep(strsplit(ops[1], "vs")[[1]][2],tempfile[as.matrix(strsplit(ops[1], "vs")[[1]][2]),"countsofsample"])), levels = c(strsplit(ops[1], "vs")[[1]][2],strsplit(ops[1], "vs")[[1]][1]))

  coldata <- data.frame(row.names=colnames(counts), condition)
  
  dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
  dds <- DESeq(dds, fitType="mean")
  results = results(dds)
  
  resultsubset <- results
  resultsubset = resultsubset[ is.na(resultsubset$padj) == 0, ]
  resultsubset = resultsubset[ resultsubset$padj < padjthreshold, ]
  resultsubset = resultsubset[ abs(resultsubset$log2FoldChange) > fcthreshold, ]
  resultoutput <- merge(as.data.frame(resultsubset), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resultoutput)[1] <- "id"
  write.table(resultoutput, file=paste(ops[1],"_Significant.xls",sep=""), sep="\t", quote=F, row.names=F)
  genelist <- c(genelist,resultoutput$id)

  jpeg(paste(ops[1], "_VolcanoPlot.jpeg", sep=""), width=600, height=600)
  plot(results[,2], -log10(results[,6]), cex=0.5, xlab="log2 fold change", ylab="log10 adjusted p-value", main = paste(ops[1], "_VolcanoPlot.jpeg", sep=""))
  points(resultsubset[,2], -log2(resultsubset[,6]), cex=0.5, col="green")
  legend("bottomright", c("genes", "Sig."), text.col = c("black", "green"))
  dev.off()

  counts.normalized = round(counts(dds, normalized=TRUE), 2)
  colnames(counts.normalized) = paste(colnames(counts.normalized), "norm", sep="_")

  jpeg(paste(ops[1], "_RawCountAboveThreshold.jpeg", sep=""), width=600, height=600)
  barplot(geneabovethreshold, main = "Count of Genes Above a Threshold (> 2 Counts)",xlab = "Samples",ylab = "Gene Count", col=c("darkblue","red"), legend = rownames(geneabovethreshold))
  dev.off()

  jpeg(paste(ops[1], "_Top20ExpressionHeatMap.jpeg", sep=""), width=600, height=600)
  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(dds)[,c("condition")])
  rownames(df) <- colnames(counts.normalized)
  colnames(df) <- "Condition"
  pheatmap(log2(counts.normalized[select,]+1), cluster_rows=FALSE, cluster_cols=FALSE, main="Top 20 Genes - Log2 Expression")
  dev.off()
  
  if(length(row.names(resultsubset))>0) {
	  jpeg(paste(ops[1], "_TopDEGenesHeatMap.jpeg", sep=""), width=600, height=600)
	  degenes <- counts(dds,normalized=TRUE)[row.names(resultsubset),]
	  if(length(row.names(resultsubset))>20) {
		degenes <- degenes[order(rowMeans(degenes),decreasing=TRUE)[1:20],]
	  } 
	  df <- as.data.frame(colData(dds)[,c("condition")])
	  rownames(df) <- colnames(counts.normalized)
	  colnames(df) <- "Condition"
	  pheatmap(log2(counts.normalized[row.names(degenes),]+1), cluster_rows=FALSE, cluster_cols=FALSE, main="Top (Upto 20) DE Genes - Log2 Expression")
	  dev.off()
  }

  jpeg(paste(ops[1], "_boxplotsofcounts.jpeg", sep=""), width=600, height=600)
  par(mfrow=c(1,2))
  boxplot(log2(counts), xlab="", ylab="Log2 raw counts",las=2,main="Raw log counts")
  boxplot(log2(counts.normalized+1), xlab="", ylab="Log2 counts",las=2,main="Normalised log counts")
  dev.off()

  jpeg(paste(ops[1], "_CorrelationHeatMap.jpeg", sep=""), width=600, height=600)
  par(mfrow=c(1,1))
  sampleDistMatrix <- as.matrix(dist(t(log2(counts.normalized+1))))
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,col=colors)
  dev.off()

  pcaData <- plotPCA(vst(dds, blind=FALSE), intgroup=c("condition"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  jpeg(paste(ops[1], "_PCAAnalysis.jpeg", sep=""), width=600, height=600)
  pcaplot <- ggplot(pcaData, aes(PC1, PC2, color=condition)) + geom_point(size=3) +  xlab(paste0("PC1: ",percentVar[1],"% variance")) +  ylab(paste0("PC2: ",percentVar[2],"% variance")) +   coord_fixed()
  print(pcaplot)
  dev.off()


  controlstring <- paste(ops[1], "_Control", sep="")
  expstring <- paste(ops[1], "_Treated", sep="")
  
  #colnames(results)[3] = sub("A", paste(".", controlstring, sep=""), colnames(results)[3])
  #colnames(results)[4] = sub("B", paste(".", expstring, sep=""), colnames(results)[4])
  #colnames(results)[5] = paste(expstring, controlstring, sep="/")
  #colnames(results)[6] = paste("log2(", paste(expstring, controlstring, sep="/"), ")", sep="")

  colnames(results) = paste(colnames(results), ops[1], sep="_")

  output.filename = paste(ops[1], ".xls", sep="")
  # Print output (including norm counts)
  
  write.table(cbind(as.data.frame(results), as.data.frame(counts), as.data.frame(counts.normalized)), file=output.filename, sep="\t", quote=F, col.names=NA, row.names = TRUE)
  
  if(i == 1){
    mergedData <- as.data.frame(results)
  }
  else {
    mergedData <- merge(mergedData, as.data.frame(results), by = "row.names", suffixes=c("",ops[1]))
    rownames(mergedData) <- mergedData[,1]
    mergedData[,1] <- NULL
  }
  
}

write.table(mergedData, file="MergedData.xls", sep="\t", quote=F, col.names=NA, row.names = TRUE)

genelist <- unique(genelist)
write.table(genelist, file="SignificantGeneList.txt", sep="\t", quote=F, row.names=F)
