

##########################
#######Index genome!######
##########################

cd /home/tcherubino/doc/genome

rm eugenioidesChromossomes.fa unplacedContigs.fa canephoraChromossomes.fa

rm *ht2

#create the unified C.a genome without unplaced contigs
cat c.c.fasta c.e.fa > C.a.fa

#create a unified GFF file

cat ~/doc/subCanephora/finalPrediction/processed.C.canephora.gff ~/doc/subEugenioides/finalPrediction/processed.C.eugenioides.gff > C.a.gff

sed -i"" 's/transcript_id /transcript_id=/g' C.a.gff

sed -i"" 's/gene_id /gene_id=/g' C.a.gff

awk '{if ($3 == "exon") print $0}' C.a.gff > compatible.gff

/usr/local/bin/STAR runThreadN 24 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ./C.a.gff --sjdbGTFtagExonParentTranscript gene_id --sjdbGTFfeatureExon exon --genomeFastaFiles ./C.a.fa --genomeSAindexNbases 13

cd /home/tcherubino/doc

mkdir G.libraries

rm *PREFILTERED*

rm *ADAPTER*

mkdir ../analysis
cd ../analysis

allPath=/home/tcherubino/doc/G.libraries/DLVR214183Lom-2/
for library in `ls ../DLVR214183Lom-2/*1.fastq.gz`
do

if [ -f "`basename -s 1.fastq.gz $library`counts.txt" ]; then
echo "`basename -s 1.fastq.gz $library`counts.txt exists."
else
echo "`basename -s 1.fastq.gz $library`counts.txt does not exist."

echo $library
pigz -p 8 -d -f -k -c $library > `basename -s .gz $library` &

#Add a 10 seconds pause just to avoid problems while decompressing the second file
sleep 10s
echo Resumming

pigz -p 8 -d -f -k -c $allPath`basename -s 1.fastq.gz $library`2.fastq.gz > `basename -s 1.fastq.gz $library`2.fastq

#aling paired end reads

/usr/local/bin/STAR --runThreadN 24 --genomeDir ../../genome/ --readFilesIn `basename -s .gz $library` `basename -s 1.fastq.gz $library`2.fastq --outFileNamePrefix  `basename -s 1.fastq.gz $library` --outSAMtype BAM Unsorted

rm `basename -s .gz $library` `basename -s 1.fastq.gz $library`2.fastq

#keep only the mapped reads in bam file
#samtools view -F 4 -@ 8 `basename -s 1.fq.gz $library`Aligned.out.bam > `basename -s 1.fq.gz $library`onlyMapped.bam
#not a good idea the last command, for some reason the resulting file was file 3 times bigger than the original

#Mark dup with picard

#The program can take either coordinate-sorted or query-sorted inputs, however the behavior is slightly different. When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates. However, when the input is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads.


java -XX:ParallelGCThreads=24 -jar ~/bin/picard.jar SortSam I=`basename -s 1.fastq.gz $library`Aligned.out.bam O=`basename -s 1.fastq.gz $library`sorted.bam SORT_ORDER=queryname

rm `basename -s 1.fastq.gz $library`Aligned.out.bam

java -XX:ParallelGCThreads=24 -jar ~/bin/picard.jar MarkDuplicates I=`basename -s 1.fastq.gz $library`sorted.bam O=`basename -s 1.fastq.gz $library`sorted.rmdup.bam REMOVE_DUPLICATES=true M=`basename -s 1.fastq.gz $library`rmdupReport.txt

rm `basename -s 1.fastq.gz $library`sorted.bam

/usr/local/bioinformatic/miniconda2/bin/htseq-count -a 10 -t exon -i gene_id -f bam --stranded=no `basename -s 1.fastq.gz $library`sorted.rmdup.bam ../../genome/compatible.gff > `basename -s 1.fastq.gz $library`counts.txt &
fi
done



############################################
##############Regulatory Networks##########
##########################################
mkdir WGCNA
cd WGCNA/
mkdir counts

export ALLOW_WGCNA_THREADS=8

R
library("edgeR")

targets <- readTargets()

names <- targets$description

#create a DGE matrix
matrix_input <- readDGE(targets, comment.char = "!")

sum(as.matrix(tail(matrix_input)[6,]))

#remove meta Tags
MetaTags <- grep("^__", rownames(matrix_input))
matrix_input <- matrix_input[-MetaTags, ]

write.table(matrix_input$counts, file ="counts.matrix.tab",sep="\t")

reads_before <- sum(matrix_input$counts)
#remove low expressed genes
rnaseqmatrix <- matrix_input$counts

#nrow(rnaseqmatrix[rowMeans(rnaseqmatrix) > 0,])
rnaseqmatrix <- rnaseqmatrix[rowMeans(rnaseqmatrix) >=20,]

rnaseqmatrix <- rnaseqmatrix[rowMeans(rnaseqmatrix) >=25,]

conditions = matrix_input$samples[,2]
colnames(rnaseqmatrix) <- names
write.table(rnaseqmatrix, file ="rnaseq.matrix.tab",sep="\t")

#/home/tcherubino/bin/R-3.5.3/bin/R

library(WGCNA)

rnaseqMatrix <- read.delim("rnaseq.matrix.tab",sep = "\t",header=T, row.names=1)

enableWGCNAThreads(8)

rnaseqMatrix <- rnaseqmatrix

analysis_matrix <- DGEList(counts = rnaseqMatrix,group = targets$group)
design <- model.matrix(~0+group, data=analysis_matrix$samples)
colnames(design) <- levels(analysis_matrix$samples$group)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

analysis_matrix <- calcNormFactors(analysis_matrix)

#To estimate common dispersion:
analysis_matrix <- estimateGLMCommonDisp(analysis_matrix, design)
#To estimate trended dispersions:
analysis_matrix <- estimateGLMTrendedDisp(analysis_matrix, design)
#To estimate tagwise dispersions:
analysis_matrix <- estimateGLMTagwiseDisp(analysis_matrix, design)

#Fit a negative binomial generalized log-linear model to the read counts for each gene.
fit <- glmFit(analysis_matrix,design)

sum(matrix_input$counts)

sum(fit$counts)/24
#The square root of dispersion is the coe cient of biological variation (BCV)

pdf(file = "edgeR_BCV.pdf")
plotBCV(analysis_matrix)
dev.off()

#An MDS plots shows distances, in terms of biological coeficient of variation (BCV) - An MDS plot shows the relative similarities of the samples.
pdf(file = "edgeR_MDS.pdf")
plotMDS(analysis_matrix)
dev.off()


#visualize the libraries sizes
pdf("libraries_sizes.pdf", wi=12,he=8)
barplot(analysis_matrix$samples$lib.size)
abline(h=mean(analysis_matrix$samples$lib.size)-2*sd(analysis_matrix$samples$lib.size),col="black",lwd = 3)
abline(h=mean(analysis_matrix$samples$lib.size),lty = 2,lwd = 3)
dev.off()


#samples_trees
cpm.matrix <- cpm(analysis_matrix,normalized.lib.sizes=TRUE)
colnames(cpm.matrix) <- names
t.cpm.matrix <- t(cpm.matrix)
sampleTree <- hclust(dist(t.cpm.matrix), method = "average");


# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
#abline(h = 15, col = "red")
dev.off()

matrix_input <- read.delim("counts.matrix.tab",sep = "\t",header=T, row.names=1)

#analyse densisties before normalization - input
pdf("Densities_input.pdf")
plotDensities( log(matrix_input), legend = "topright")
dev.off()

#analyse densisties low expression filter
pdf("Densities_low_expression_fiter.pdf")
plotDensities( log(rnaseqMatrix), legend = "topright")
dev.off()

#analyse densisties after normalization
pdf("Densities_normalization.pdf")
plotDensities( log(cpm.matrix), legend = "topright")
dev.off()

pdf("Densities_log_cpm_fitted_norm.pdf")
plotDensities(log(cpm(fit$fitted.values, normalized.lib.sizes=TRUE)), legend = "topright")
dev.off()

#histogram of densities log10
pdf("Log10_histogram_normilized.pdf", h=10,w=10)
hist(log(cpm.matrix+1,10), col=gray.colors(19, start = 0.9, end = 0.3))
dev.off()

#histogram of densities Log2
pdf("Log2_histogram_normilized.pdf", h=10,w=10)
hist(log(cpm.matrix+1,2), col=gray.colors(19, start = 0.9, end = 0.3))
dev.off()

#histogram of densities no log
pdf("histogram_normilized.pdf", h=10,w=10)
hist(cpm.matrix, col=gray.colors(19, start = 0.9, end = 0.3))
dev.off()

# Load library for pheatmap
cpm.matrix.corr <- cor(cpm.matrix, method="spearman",use="pairwise.complete.obs")
colnames(cpm.matrix.corr) <- fit$samples$group
rownames(cpm.matrix.corr) <- fit$samples$group
library(pheatmap)
library(gplots)
pdf(file = "sampleClusteringHeatmap.pdf", width = 35, height = 25)
pheatmap(cpm.matrix.corr, fontsize=30,cellwidth=70, cellheight=60, treeheight_col= 200, treeheight_row=200,angle_col=0, legend=T)
#heatmap.2(cpm.matrix.corr*10000,trace = "none",margins = c(5, 11),keysize =1 , key.title="",col ="bluered",density.info="none")
dev.off()

###########################################################


cpm.matrix <- cpm(analysis_matrix,normalized.lib.sizes=FALSE)
cpm.matrix <- log(cpm.matrix+1,2)
cpm.matrix <- t(cpm.matrix)


#We first check for genes and samples with too many missing values:

gsg = goodSamplesGenes(cpm.matrix, verbose = 3);
gsg$allOK


#If the last statement returns TRUE , all genes have passed the cuts. If not, we remove the offending genes and samplesfrom the data:
if (!gsg$allOK)
{
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(cpm.matrix)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(cpm.matrix)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
cpm.matrix0 = cpm.matrix[gsg$goodSamples, gsg$goodGenes]
}

#Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.

sampleTree = hclust(dist(cpm.matrix), method = "average");


# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 15, col = "red")
dev.off()


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(cpm.matrix, powerVector = powers, verbose = 5)

# Plot the results:
pdf("pick_soft_treshold.pdf", width = 12, height = 9)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower = 6;

adjacency = adjacency(cpm.matrix, power = softPower);

#Topological Overlap Matrix (TOM)
#To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

#Clustering using TOM

#We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes. Note that we use the function hclust that provides a much faster hierarchical clustering routine than the standard hclust function.

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
pdf("Gene_clustering_on_TOM-based_dissimilarity.pdf", width = 12, height = 9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04);
dev.off()


#The clustering dendrogram plotted by the last command is shown in Figure 2. In the clustering tree (dendrogram), each leaf, that is a short vertical line, corresponds to a gene. Branches of the dendrogram group together densely interconnected, highly co-expressed genes. Module identification amounts to the identification of individual branches (”cutting the branches off the dendrogram”). There are several methods for branch cutting; our standard method is the Dynamic Tree Cut from the package dynamicTreeCut. The next snippet of code illustrates its use.

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
#I choosed to set this parameter to 100 to reduce the total number of modules - we are reporting more than 20 here
#minModuleSize = 100;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath

pdf("Gene_clustering_on_TOM-based_dissimilarity_colors.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
dev.off()

#Merging of modules whose expression profiles are very similar
#The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent to merge such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire modules, we calculate their eigengenes and cluster them on their correlation:

# Calculate eigengenes
MEList = moduleEigengenes(cpm.matrix, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf("Clustering_of_module_eigengenes.pdf", width = 10, height = 6)

plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

#We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge = mergeCloseModules(cpm.matrix, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

#To see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged module colors underneath

pdf(file = "geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
#colorOrder = c("grey", standardColors(50));
colorOrder = c(standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

####################################################
###Visualizing the gene network#####################
####################################################

#One way to visualize a weighted network is to plot its heatmap, Fig. 1. Each row and column of the heatmap correspond to a single gene. The heatmap can depict adjacencies or topological overlaps, with light colors denoting low adjacency (overlap) and darker colors higher adjacency (overlap). In addition, the gene dendrograms and module colors are plotted along the top and left side of the heatmap. The package provides a convenient function to create such network plots; Fig. 1 was created using the following code. This code can be executed only if the network was calculated using a single-block approach (that is, using the 1-step automatic or the step-by-step tutorials). If the networks were calculated using the block-wise approach, the user will need to modify this code to perform the visualization in each block separately. The modification is simple and we leave it as an exercise for the interested reader.

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
#dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function


png(file = "network_heatmap_4000.png", width = 4000, height = 4000);
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes", margin = c(50,50))
dev.off()


#terrainColors
png(file = "network_heatmap_4000terrainCorlors.png", width = 4000, height = 4000);
TOMplot(plotTOM, geneTree, moduleColors,terrainColors=T,  main = "Network heatmap plot, all genes", margin = c(50,50))
dev.off()


#Visualizing the network of eigengenes
#It is often interesting to study the relationships among the found modules. One can use the eigengenes as representative profiles and quantify module similarity by eigengene correlation. The package contains a convenient function plotEigengeneNetworks that generates a summary plot of the eigengene network. It is usually informative to add a clinical trait (or multiple traits) to the eigengenes to see how the traits fit into the eigengene network:

# Recalculate module eigengenes
#MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
#weight = as.data.frame(datTraits$weight_g);
#names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(MEs)
# Plot the relationships among the eigengenes and the trait
#sizeGrWindow(5,7.5);
pdf(file = "network_eigengenes_heatmap.pdf", width = 5, height = 7.5);

par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
= 90)
dev.off()

pdf("TOMhistogram.pdf")
plot(hist(TOM))
dev.off()


pdf("TOMdensity.pdf")
plot(density(TOM))
dev.off()

####################################
#export group members ID and annot#
##################################
descriptions <- read.delim("Gene.description.tab",header=F)
probes = colnames(cpm.matrix)

#antiquewhite4
inModule = is.finite(match(moduleColors, "antiquewhite4"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="antiquewhiteDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#brown4
inModule = is.finite(match(moduleColors, "brown4"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="brownDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#coral2
inModule = is.finite(match(moduleColors, "coral2"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="coralDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#cyan
inModule = is.finite(match(moduleColors, "cyan"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="cyanDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#darkorange
inModule = is.finite(match(moduleColors, "darkorange"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="darkorangeDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#darkred
inModule = is.finite(match(moduleColors, "darkred"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="darkredDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#darkseagreen4
inModule = is.finite(match(moduleColors, "darkseagreen4"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="darkseagreenDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#darkslateblue
inModule = is.finite(match(moduleColors, "darkslateblue"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="darkslateblueDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#greenyellow
inModule = is.finite(match(moduleColors, "greenyellow"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="greenyellowDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#grey60
inModule = is.finite(match(moduleColors, "grey60"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="greyDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#honeydew1
inModule = is.finite(match(moduleColors, "honeydew1"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="honeydewDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#lavenderblush3

inModule = is.finite(match(moduleColors, "lavenderblush3"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="lavenderblushDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#lightpink4

inModule = is.finite(match(moduleColors, "lightpink4"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="lightpink4Description.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#lightsteelblue1

inModule = is.finite(match(moduleColors, "lightsteelblue1"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="lightsteelblue1Description.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#orange

inModule = is.finite(match(moduleColors, "orange"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="orangeDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#pink

inModule = is.finite(match(moduleColors, "pink"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="pinkDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#plum2

inModule = is.finite(match(moduleColors, "plum2"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="plum2Description.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#salmon4

inModule = is.finite(match(moduleColors, "salmon4"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="salmon4Description.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#sienna3

inModule = is.finite(match(moduleColors, "sienna3"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="sienna3Description.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#tan

inModule = is.finite(match(moduleColors, "tan"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="tanDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#thistle1

inModule = is.finite(match(moduleColors, "thistle1"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="thistle1Description.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#turquoise

inModule = is.finite(match(moduleColors, "turquoise"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="turquoiseDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#violet

inModule = is.finite(match(moduleColors, "violet"));
inModuleDescriptions <- descriptions[is.finite(match(descriptions[,1],probes[inModule])),]
write.table(inModuleDescriptions,file="violetDescription.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)


##################################
#####Export gene ontology Terms###
##################################

GO.terms <- read.delim("GO.one.per.row.tab",header=F)

#antiquewhite4
inModule = is.finite(match(moduleColors, "antiquewhite4"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="antiquewhiteGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#brown4
inModule = is.finite(match(moduleColors, "brown4"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="brownGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)


table(moduleColors)

#coral2
inModule = is.finite(match(moduleColors, "coral2"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="coralGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)


table(moduleColors)

#cyan
inModule = is.finite(match(moduleColors, "cyan"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="cyanGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#darkorange
inModule = is.finite(match(moduleColors, "darkorange"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="darkorangeGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#darkred
inModule = is.finite(match(moduleColors, "darkred"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="darkredGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#darkseagreen4
inModule = is.finite(match(moduleColors, "darkseagreen4"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="darkseagreenGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#darkslateblue
inModule = is.finite(match(moduleColors, "darkslateblue"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="darkslateblueGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#greenyellow
inModule = is.finite(match(moduleColors, "greenyellow"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="greenyellowGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#grey60
inModule = is.finite(match(moduleColors, "grey60"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="greyGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

table(moduleColors)

#honeydew1
inModule = is.finite(match(moduleColors, "honeydew1"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="honeydewGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#lavenderblush3

inModule = is.finite(match(moduleColors, "lavenderblush3"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="lavenderblushGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#lightpink4

inModule = is.finite(match(moduleColors, "lightpink4"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="lightpink4GO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#lightsteelblue1

inModule = is.finite(match(moduleColors, "lightsteelblue1"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="lightsteelblue1GO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#orange

inModule = is.finite(match(moduleColors, "orange"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="orangeGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#pink

inModule = is.finite(match(moduleColors, "pink"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="pinkGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#plum2

inModule = is.finite(match(moduleColors, "plum2"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="plum2GO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#salmon4

inModule = is.finite(match(moduleColors, "salmon4"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="salmon4GO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#sienna3

inModule = is.finite(match(moduleColors, "sienna3"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="sienna3GO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#tan

inModule = is.finite(match(moduleColors, "tan"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="tanGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#thistle1

inModule = is.finite(match(moduleColors, "thistle1"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="thistle1GO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#turquoise

inModule = is.finite(match(moduleColors, "turquoise"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="turquoiseGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)

#violet

inModule = is.finite(match(moduleColors, "violet"));
inModuleGO <- GO.terms[is.finite(match(GO.terms[,1],probes[inModule])),]
write.table(inModuleGO,file="violetGO.tab",sep="\t", row.names = T, col.names=FALSE,quote=F)


###########################################################################
########6. Exporting a gene network to external visualization software#####
###########################################################################


#########################################
#GroupA
#########################################

probes = colnames(cpm.matrix)

modules_GroupA = c("plum2", "darkslateblue","thistle1")

inModule = is.finite(match(moduleColors, modules_GroupA));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

summary(adj[adj != 1])

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,1))

#lets just consider whats above the 95th Qu.
adj[adj > 0.049] = 1
adj[adj != 1] = 0

(network <- graph.adjacency(adj))
(network <- simplify(network))  # removes self-loops
(network <- delete.vertices(network, degree(network)==0))


ColorTable <- cbind(probes,mergedColors)

V(network)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(network)))),2]

# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)

#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network. It can be thought of as how critical the vertex is to the flow of information through a network. Individuals with high betweenness are key bridges between different parts of a network. In our measles transmission network, vertices with high betweenness are those children who were central to passing on the disease to other parts of the network. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the network adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(network, directed = TRUE)

# Show histogram of vertex betweenness
pdf("vertexBetweeness_groupA.pdf")
hist(g.b, breaks = 80)
dev.off()


png(file = "nicely.network_groupA.png", width = 4000, height = 4000);
plot(network,layout=layout_nicely(network), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()


png(file = "fr.network_groupA.png", width = 4000, height = 4000);
plot(network,layout=layout_with_fr(network), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_groupA.png", width = 4000, height = 4000);
plot(network,layout=layout_as_tree(network), edge.arrow.size = 0.05)
dev.off()



#########################################
#GroupB
#########################################

probes = colnames(cpm.matrix)

modules_GroupB = c("sienna3", "tan")

inModule = is.finite(match(moduleColors, modules_GroupB));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,1))

#lets just consider whats above the 95th Qu.
adj[adj > 0.13] = 1
adj[adj != 1] = 0

(network <- graph.adjacency(adj))
(network <- simplify(network))  # removes self-loops
(network <- delete.vertices(network, degree(network)==0))

ColorTable <- cbind(probes,mergedColors)

V(network)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(network)))),2]


#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network. It can be thought of as how critical the vertex is to the flow of information through a network. Individuals with high betweenness are key bridges between different parts of a network. In our measles transmission network, vertices with high betweenness are those children who were central to passing on the disease to other parts of the network. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the network adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(network, directed = TRUE)

# Show histogram of vertex betweenness
pdf("vertexBetweenessB.pdf")
hist(g.b, breaks = 80)
dev.off()

png(file = "nicely.network_groupB.png", width = 4000, height = 4000);
plot(network,layout=layout_nicely(network), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_groupB.png", width = 4000, height = 4000);
plot(network,layout=layout_with_fr(network), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_groupB.png", width = 4000, height = 4000);
plot(network,layout=layout_as_tree(network), edge.arrow.size = 0.05)
dev.off()

#######################
#GroupC################
#######################


probes = colnames(cpm.matrix)

modules_GroupD = c("coral2", "darkred","turquoise")

inModule = is.finite(match(moduleColors, modules_GroupC));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,1))

#lets just consider whats above the 90th Qu.
adj[adj > 0.202] = 1
adj[adj != 1] = 0

(network <- graph.adjacency(adj))
(network <- simplify(network))  # removes self-loops
(network <- delete.vertices(network, degree(network)==0))

ColorTable <- cbind(probes,mergedColors)

V(network)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(network)))),2]


#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network. It can be thought of as how critical the vertex is to the flow of information through a network. Individuals with high betweenness are key bridges between different parts of a network. In our measles transmission network, vertices with high betweenness are those children who were central to passing on the disease to other parts of the network. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the network adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(network, directed = TRUE)

# Show histogram of vertex betweenness
pdf("vertexBetweenessGroupC.pdf")
hist(g.b, breaks = 80)
dev.off()


png(file = "nicely.network_groupC.png", width = 4000, height = 4000);
plot(network,layout=layout_nicely(network), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_groupC.png", width = 4000, height = 4000);
plot(network,layout=layout_with_fr(network), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_groupC.png", width = 4000, height = 4000);
plot(network,layout=layout_as_tree(network), edge.arrow.size = 0.05)
dev.off()


#######################
#GroupD################
#######################


probes = colnames(cpm.matrix)

modules_GroupC = c("greenyellow", "orange","honeydew1")

inModule = is.finite(match(moduleColors, modules_GroupD));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

summary(adj[adj != 1])

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,1))

#lets just consider whats above the 95th Qu.
adj[adj > 0.067] = 1
adj[adj != 1] = 0

(network <- graph.adjacency(adj))
(network <- simplify(network))  # removes self-loops
(network <- delete.vertices(network, degree(network)==0))


ColorTable <- cbind(probes,mergedColors)

V(network)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(network)))),2]


#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network. It can be thought of as how critical the vertex is to the flow of information through a network. Individuals with high betweenness are key bridges between different parts of a network. In our measles transmission network, vertices with high betweenness are those children who were central to passing on the disease to other parts of the network. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the network adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(network, directed = TRUE)


# Show histogram of vertex betweenness
pdf("vertexBetweenessGroupD.pdf")
hist(g.b, breaks = 80)
dev.off()

png(file = "nicely.network_groupD.png", width = 4000, height = 4000);
plot(network,layout=layout_nicely(network), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_groupD.png", width = 4000, height = 4000);
plot(network,layout=layout_with_fr(network), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_groupD.png", width = 4000, height = 4000);
plot(network,layout=layout_as_tree(network), edge.arrow.size = 0.05)
dev.off()


#######################
#GroupE################
#######################

probes = colnames(cpm.matrix)

modules_GroupE = c("cyan", "lightsteelblue1","antiquewhite4","violet")

inModule = is.finite(match(moduleColors, modules_GroupE));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,1))

#lets just consider whats above the 95th Qu.
adj[adj > 0.03] = 1
adj[adj != 1] = 0

(network <- graph.adjacency(adj))
(network <- simplify(network,  degree(network)==0))  # removes self-loops
(network <- delete.vertices(network, degree(network)==0))


ColorTable <- cbind(probes,mergedColors)

V(network)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(network)))),2]

#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network. It can be thought of as how critical the vertex is to the flow of information through a network. Individuals with high betweenness are key bridges between different parts of a network. In our measles transmission network, vertices with high betweenness are those children who were central to passing on the disease to other parts of the network. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the network adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(network, directed = TRUE)

# Show histogram of vertex betweenness
pdf("vertexBetweenessGroupE.pdf")
hist(g.b, breaks = 80)
dev.off()

png(file = "nicely.network_groupE.png", width = 4000, height = 4000);
plot(network,layout=layout_nicely(network), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()



png(file = "fr.network_groupE.png", width = 4000, height = 4000);
plot(network,layout=layout_with_fr(network), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_groupE.png", width = 4000, height = 4000);
plot(network,layout=layout_as_tree(network), edge.arrow.size = 0.05)
dev.off()



#######################
#GroupF################
#######################


probes = colnames(cpm.matrix)

modules_GroupF = c("lightpink4", "pink")

inModule = is.finite(match(moduleColors, modules_GroupF));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,.96,.97,.98,.99,1))

#lets just consider whats above the 95th Qu.
adj[adj > 0.14] = 1
adj[adj != 1] = 0

(network <- graph.adjacency(adj))
(network <- simplify(network,  degree(network)==0))  # removes self-loops

(network <- delete.vertices(network, degree(network)==0))

ColorTable <- cbind(probes,mergedColors)

V(network)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(network)))),2]

#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network. It can be thought of as how critical the vertex is to the flow of information through a network. Individuals with high betweenness are key bridges between different parts of a network. In our measles transmission network, vertices with high betweenness are those children who were central to passing on the disease to other parts of the network. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the network adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(network, directed = TRUE)

# Show histogram of vertex betweenness
pdf("vertexBetweenessGroupF.pdf")
hist(g.b, breaks = 80)
dev.off()

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_groupF.png", width = 4000, height = 4000);
plot(network,layout=layout_nicely(network), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_groupF.png", width = 4000, height = 4000);
plot(network,layout=layout_with_fr(network), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_groupF.png", width = 4000, height = 4000);
plot(network,layout=layout_as_tree(network), edge.arrow.size = 0.05)
dev.off()

#######################
#GroupG################
#######################


probes = colnames(cpm.matrix)

modules_GroupG = c("brown4", "lavenderblush3")

inModule = is.finite(match(moduleColors, modules_GroupG));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,.96,.97,.98,.99,1))

#lets just consider whats above the 95th Qu.
adj[adj > 0.048] = 1
adj[adj != 1] = 0

(network <- graph.adjacency(adj))
(network <- simplify(network,  degree(network)==0))  # removes self-loops

network <- delete.vertices(network, degree(network)==0)

ColorTable <- cbind(probes,mergedColors)

V(network)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(network)))),2]

#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network. It can be thought of as how critical the vertex is to the flow of information through a network. Individuals with high betweenness are key bridges between different parts of a network. In our measles transmission network, vertices with high betweenness are those children who were central to passing on the disease to other parts of the network. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the network adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(network, directed = TRUE)

# Show histogram of vertex betweenness
pdf("vertexBetweenessGroupG.pdf")
hist(g.b, breaks = 80)
dev.off()
png(file = "nicely.network_groupG.png", width = 4000, height = 4000);
plot(network,layout=layout_nicely(network), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_groupG.png", width = 4000, height = 4000);
plot(network,layout=layout_with_fr(network), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_groupG.png", width = 4000, height = 4000);
plot(network,layout=layout_as_tree(network), edge.arrow.size = 0.05)
dev.off()


#######################
#GroupH################
#######################


probes = colnames(cpm.matrix)

modules_GroupH = c("salmon4", "darkseagreen4", "darkorange", "grey60")

inModule = is.finite(match(moduleColors, modules_GroupH));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,1))

#lets just consider whats above the 95th Qu.
adj[adj > 0.04] = 1
adj[adj != 1] = 0

(network <- graph.adjacency(adj))
(network <- simplify(network,  degree(network)==0))  # removes self-loops

(network <- delete.vertices(network, degree(network)==0))

ColorTable <- cbind(probes,mergedColors)

V(network)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(network)))),2]

#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network. It can be thought of as how critical the vertex is to the flow of information through a network. Individuals with high betweenness are key bridges between different parts of a network. In our measles transmission network, vertices with high betweenness are those children who were central to passing on the disease to other parts of the network. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the network adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(network, directed = TRUE)

# Show histogram of vertex betweenness
pdf("vertexBetweenessGroupH.pdf")
hist(g.b, breaks = 80)
dev.off()

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_groupH.png", width = 4000, height = 4000);
plot(network,layout=layout_nicely(network), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_groupH.png", width = 4000, height = 4000);
plot(network,layout=layout_with_fr(network), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_groupH.png", width = 4000, height = 4000);
plot(network,layout=layout_as_tree(network), edge.arrow.size = 0.05)
dev.off()

#######################
# all Groups ##############
#######################


probes = colnames(cpm.matrix)

modules_Group= moduleColors

inModule = is.finite(match(moduleColors, modules_Group));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]

dimnames(adj) <- list(modProbes,modProbes)

quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1))

#lets just consider whats above the 99th Qu.
adj[adj > 0.17] = 1
adj[adj != 1] = 0

(network <- graph.adjacency(adj))
(network <- simplify(network))  # removes self-loops
network <- delete.vertices(network, degree(network)==0)


ColorTable <- cbind(probes,mergedColors)

V(network)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(network)))),2]

#Betweenness

#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network. It can be thought of as how critical the vertex is to the flow of information through a network. Individuals with high betweenness are key bridges between different parts of a network. In our measles transmission network, vertices with high betweenness are those children who were central to passing on the disease to other parts of the network. In this exercise, you will identify the betweenness score for each vertex and then make a new plot of the network adjusting the vertex size by its betweenness score to highlight these key vertices.

g.b <- betweenness(network, directed = TRUE)

# Show histogram of vertex betweenness
pdf("vertexBetweeness.pdf")
hist(g.b, breaks = 80)
dev.off()

png(file = "nicely.network_ALL.groups.png", width = 4000, height = 4000);
plot(network,layout=layout_nicely(network), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()
png(file = "fr.network_ALL.groups.png", width = 4000, height = 4000);
plot(network,layout=layout_with_fr(network), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_ALL.groups.png", width = 4000, height = 4000);
plot(network,layout=layout_as_tree(network), edge.arrow.size = 0.05)
dev.off()

# Plot threejs plot of graph setting vertex size to v

library(threejs)

#png(file = "3D.network_groupABC.png", width = 4000, height = 4000);
#rglplot(network)
#dev.off()

/Users/johnjoyce/Dropbox/Bioinformatica/pipelines/C.c_subgenomePred_alingments_and_denovo_and_network.sh


#######################
####Custon analysis####
#######################

library("igraph")
#plum2 module

probes = colnames(cpm.matrix)

plum2 = "plum2"

inModule = is.finite(match(moduleColors, plum2));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.038] = 0

(networkplum2 <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkplum2 <- simplify(networkplum2,  degree(networkplum2)==0))  # removes self-loops

((networkplum2 <- delete.vertices(networkplum2, degree(networkplum2)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkplum2)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkplum2)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_plum2.png", width = 4000, height = 4000);
plot(networkplum2,layout=layout_nicely(networkplum2), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_plum2.png", width = 4000, height = 4000);
plot(networkplum2,layout=layout_with_fr(networkplum2), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_plum2.png", width = 4000, height = 4000);
plot(networkplum2,layout=layout_as_tree(networkplum2), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_plum2.png", width = 4000, height = 4000);
plot(networkplum2,layout=layout_in_circle(networkplum2), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkplum2)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkplum2, directed = FALSE)
get_diameter(networkplum2, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkplum2)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkplum2)

# Get the average path length of the graph g
g.apl <- mean_distance(networkplum2, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an important method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph
Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkplum2)

gorder(networkplum2)

g.random <- erdos.renyi.game(n = gorder(networkplum2), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random
mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networkplum2), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthplum2.pdf")
hist(gl.apls,xlim=c(1.6,1.8),col="plum2",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.005, y= 150,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,2),col="red",line=.46,cex=.9)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkplum2, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenessplum2.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="plum2", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkplum2, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreeplum2.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="plum2", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkplum2)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkplum2 <- make_ego_graph(networkplum2, 1, "SubC.e_24520", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

pOrfans <- c("SubC.c_4068","SubC.c_25356","SubC.c_29281", "SubC.e_10567", "SubC.e_18855","SubC.e_38581")

V(egoNetworkplum2)$color <- ifelse( names(V(egoNetworkplum2)) %in% pOrfans, "black", "plum2")

V(egoNetworkplum2)[names(V(egoNetworkplum2)) %in% pOrfans]

png(file = "nicely.network_ego_plum2.png", width = 4000, height = 4000);
plot(egoNetworkplum2,layout=layout_nicely(egoNetworkplum2), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#######################

#darkslateblue module

probes = colnames(cpm.matrix)

darkslateblue = "darkslateblue"

inModule = is.finite(match(moduleColors, darkslateblue));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.039] = 0

(networkdarkslateblue <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkdarkslateblue <- simplify(networkdarkslateblue,  degree(networkdarkslateblue)==0))  # removes self-loops

((networkdarkslateblue <- delete.vertices(networkdarkslateblue, degree(networkdarkslateblue)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkdarkslateblue)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkdarkslateblue)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_darkslateblue.png", width = 4000, height = 4000);
plot(networkdarkslateblue,layout=layout_nicely(networkdarkslateblue), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_darkslateblue.png", width = 4000, height = 4000);
plot(networkdarkslateblue,layout=layout_with_fr(networkdarkslateblue), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_darkslateblue.png", width = 4000, height = 4000);
plot(networkdarkslateblue,layout=layout_as_tree(networkdarkslateblue), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_darkslateblue.png", width = 4000, height = 4000);
plot(networkdarkslateblue,layout=layout_in_circle(networkdarkslateblue), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkdarkslateblue)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkdarkslateblue, directed = FALSE)
get_diameter(networkdarkslateblue, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkdarkslateblue)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkdarkslateblue)

# Get the average path length of the graph g
g.apl <- mean_distance(networkdarkslateblue, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an important method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph
Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkdarkslateblue)

gorder(networkdarkslateblue)

g.random <- erdos.renyi.game(n = gorder(networkdarkslateblue), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random
mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networkdarkslateblue), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthDarkslateblue.pdf")
hist(gl.apls,xlim=c(1.6,1.8),col="darkslateblue",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.005, y= 150,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,2),col="red",line=.46,cex=.9)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkdarkslateblue, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenessDarkslateblue.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="darkslateblue", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkdarkslateblue, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreeDarkslateblue.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="darkslateblue", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkdarkslateblue)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkdarkslateblue <- make_ego_graph(networkdarkslateblue, 1, "SubC.e_38709", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

pOrfans <- c("SubC.c_418","SubC.c_1367","SubC.c_6143", "SubC.c_24693", "SubC.c_29304", "SubC.c_29702","SubC.e_15337","SubC.e_15435","SubC.e_17236", "SubC.e_21997", "SubC.e_34722", "SubC.e_35277")

V(egoNetworkdarkslateblue)$color <- ifelse( names(V(egoNetworkdarkslateblue)) %in% pOrfans, "black", "darkslateblue")

V(egoNetworkdarkslateblue)[names(V(egoNetworkdarkslateblue)) %in% pOrfans]

png(file = "nicely.network_ego_darkslateblue.png", width = 4000, height = 4000);
plot(egoNetworkdarkslateblue,layout=layout_nicely(egoNetworkdarkslateblue), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()


#########################################

#thistle1 module

probes = colnames(cpm.matrix)

thistle1 = "thistle1"

inModule = is.finite(match(moduleColors, thistle1));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.044] = 0

(networkthistle1 <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkthistle1 <- simplify(networkthistle1,  degree(networkthistle1)==0))  # removes self-loops

((networkthistle1 <- delete.vertices(networkthistle1, degree(networkthistle1)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkthistle1)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkthistle1)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_thistle1.png", width = 4000, height = 4000);
plot(networkthistle1,layout=layout_nicely(networkthistle1), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_thistle1.png", width = 4000, height = 4000);
plot(networkthistle1,layout=layout_with_fr(networkthistle1), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_thistle1.png", width = 4000, height = 4000);
plot(networkthistle1,layout=layout_as_tree(networkthistle1), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_thistle1.png", width = 4000, height = 4000);
plot(networkthistle1,layout=layout_in_circle(networkthistle1), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkthistle1)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkthistle1, directed = FALSE)
get_diameter(networkthistle1, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkthistle1)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkthistle1)

# Get the average path length of the graph g
g.apl <- mean_distance(networkthistle1, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an important method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph
Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkthistle1)

gorder(networkthistle1)

g.random <- erdos.renyi.game(n = gorder(networkthistle1), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random
mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networkthistle1), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengththistle1.pdf")
hist(gl.apls,xlim=c(1.6,2.6),col="thistle1",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.02, y= 100,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,2),col="red",line=.46,cex=.9)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkthistle1, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenessthistle1.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="thistle1", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkthistle1, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreethistle1.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="thistle1", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkthistle1)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkthistle1 <- make_ego_graph(networkthistle1, 1, "SubC.c_13014", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

pOrfans <- c("SubC.c_20186","SubC.c_20494","SubC.c_12173", "SubC.c_13014", "SubC.c_20726","SubC.e_16257", "SubC.e_38916")

V(egoNetworkthistle1)$color <- ifelse( names(V(egoNetworkthistle1)) %in% pOrfans, "black", "thistle1")

V(egoNetworkthistle1)[names(V(egoNetworkthistle1)) %in% pOrfans]

png(file = "nicely.network_ego_thistle1.png", width = 4000, height = 4000);
plot(egoNetworkthistle1,layout=layout_nicely(egoNetworkthistle1), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#####################################

#sienna3 module

probes = colnames(cpm.matrix)

sienna3 = "sienna3"

inModule = is.finite(match(moduleColors, sienna3));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.088] = 0

(networksienna3 <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networksienna3 <- simplify(networksienna3,  degree(networksienna3)==0))  # removes self-loops

((networksienna3 <- delete.vertices(networksienna3, degree(networksienna3)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networksienna3)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networksienna3)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_sienna3.png", width = 4000, height = 4000);
plot(networksienna3,layout=layout_nicely(networksienna3), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_sienna3.png", width = 4000, height = 4000);
plot(networksienna3,layout=layout_with_fr(networksienna3), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_sienna3.png", width = 4000, height = 4000);
plot(networksienna3,layout=layout_as_tree(networksienna3), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_sienna3.png", width = 4000, height = 4000);
plot(networksienna3,layout=layout_in_circle(networksienna3), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networksienna3)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networksienna3, directed = FALSE)
get_diameter(networksienna3, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networksienna3)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networksienna3)

# Get the average path length of the graph g
g.apl <- mean_distance(networksienna3, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an important method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph
Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networksienna3)

gorder(networksienna3)

g.random <- erdos.renyi.game(n = gorder(networksienna3), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random
mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networksienna3), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthsienna3.pdf")
hist(gl.apls,xlim=c(1.6,2.6),col="sienna3",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.02, y= 100,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,2),col="red",line=.46,cex=.9)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networksienna3, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenesssienna3.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="sienna3", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networksienna3, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreesienna3.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="sienna3", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networksienna3)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworksienna3 <- make_ego_graph(networksienna3, 1, "SubC.c_6100", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

pOrfans <- c("SubC.c_10229","SubC.c_12080","SubC.c_12583", "SubC.c_20338", "SubC.c_28103","SubC.c_28461", "SubC.c_29405", "SubC.e_841","SubC.e_2040", "SubC.e_3731","SubC.e_6951", "SubC.e_11133", "SubC.e_16078", "SubC.e_16169", "SubC.e_16795","SubC.e_16936", "SubC.e_24637")

V(egoNetworksienna3)$color <- ifelse( names(V(egoNetworksienna3)) %in% pOrfans, "black", "sienna3")

V(egoNetworksienna3)[names(V(egoNetworksienna3)) %in% pOrfans]

png(file = "nicely.network_ego_sienna3.png", width = 4000, height = 4000);
plot(egoNetworksienna3,layout=layout_nicely(egoNetworksienna3), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#####################################

#tan module

probes = colnames(cpm.matrix)

tan = "tan"

inModule = is.finite(match(moduleColors, tan));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.093] = 0

(networktan <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networktan <- simplify(networktan,  degree(networktan)==0))  # removes self-loops

((networktan <- delete.vertices(networktan, degree(networktan)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networktan)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networktan)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_tan.png", width = 4000, height = 4000);
plot(networktan,layout=layout_nicely(networktan), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_tan.png", width = 4000, height = 4000);
plot(networktan,layout=layout_with_fr(networktan), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_tan.png", width = 4000, height = 4000);
plot(networktan,layout=layout_as_tree(networktan), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_tan.png", width = 4000, height = 4000);
plot(networktan,layout=layout_in_circle(networktan), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networktan)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networktan, directed = FALSE)
get_diameter(networktan, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networktan)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networktan)

# Get the average path length of the graph g
g.apl <- mean_distance(networktan, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an important method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph
Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networktan)

gorder(networktan)

g.random <- erdos.renyi.game(n = gorder(networktan), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random
mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networktan), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthtan.pdf")
hist(gl.apls,xlim=c(1.45,1.5),col="tan",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl +0.0015, y= 150,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.46,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networktan, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenesstan.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="tan", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networktan, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreetan.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="tan", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networktan)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworktan <- make_ego_graph(networktan, 1, "SubC.c_21050", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene


temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworktan)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworktan)$color <- ifelse( names(V(egoNetworktan)) %in% pOrfans$V1, "black", "tan")


na.omit(pOrfans[match(names(V(egoNetworktan)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworktan)), pOrfans$V1),]),file = "tanHighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_tan.png", width = 4000, height = 4000);
plot(egoNetworktan,layout=layout_nicely(egoNetworktan), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#####################################

#coral2 module

probes = colnames(cpm.matrix)

coral2 = "coral2"

inModule = is.finite(match(moduleColors, coral2));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.052] = 0

(networkcoral2 <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkcoral2 <- simplify(networkcoral2,  degree(networkcoral2)==0))  # removes self-loops

((networkcoral2 <- delete.vertices(networkcoral2, degree(networkcoral2)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkcoral2)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkcoral2)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_coral2.png", width = 4000, height = 4000);
plot(networkcoral2,layout=layout_nicely(networkcoral2), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_coral2.png", width = 4000, height = 4000);
plot(networkcoral2,layout=layout_with_fr(networkcoral2), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_coral2.png", width = 4000, height = 4000);
plot(networkcoral2,layout=layout_as_tree(networkcoral2), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_coral2.png", width = 4000, height = 4000);
plot(networkcoral2,layout=layout_in_circle(networkcoral2), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkcoral2)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkcoral2, directed = FALSE)
get_diameter(networkcoral2, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkcoral2)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkcoral2)

# Get the average path length of the graph g
g.apl <- mean_distance(networkcoral2, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an imporcoral2t method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph
Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkcoral2)

gorder(networkcoral2)

g.random <- erdos.renyi.game(n = gorder(networkcoral2), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networkcoral2), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthcoral2.pdf")
hist(gl.apls,xlim=c(1.45,1.8),col="coral2",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.007, y= 100,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.46,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the imporcoral2ce of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkcoral2, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenesscoral2.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="coral2", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic discoral2ce between vertices. Vertices that are connected to each other have a geodesic discoral2ce of 1. Those that share a neighbor in common but are not connected to each other have a geodesic discoral2ce of 2 and so on.

#Perhaps the most straightforward measure of vertex imporcoral2ce is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkcoral2, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreecoral2.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="coral2", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkcoral2)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkcoral2 <- make_ego_graph(networkcoral2, 1, "SubC.e_5720", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkcoral2)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworkcoral2)$color <- ifelse( names(V(egoNetworkcoral2)) %in% pOrfans$V1, "black", "coral2")


na.omit(pOrfans[match(names(V(egoNetworkcoral2)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworkcoral2)), pOrfans$V1),]),file = "coral2HighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_coral2.png", width = 4000, height = 4000);
plot(egoNetworkcoral2,layout=layout_nicely(egoNetworkcoral2), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#####################################

#darkred module

probes = colnames(cpm.matrix)

darkred = "darkred"

inModule = is.finite(match(moduleColors, darkred));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.046] = 0

(networkdarkred <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkdarkred <- simplify(networkdarkred,  degree(networkdarkred)==0))  # removes self-loops

((networkdarkred <- delete.vertices(networkdarkred, degree(networkdarkred)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkdarkred)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkdarkred)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_darkred.png", width = 4000, height = 4000);
plot(networkdarkred,layout=layout_nicely(networkdarkred), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_darkred.png", width = 4000, height = 4000);
plot(networkdarkred,layout=layout_with_fr(networkdarkred), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_darkred.png", width = 4000, height = 4000);
plot(networkdarkred,layout=layout_as_tree(networkdarkred), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_darkred.png", width = 4000, height = 4000);
plot(networkdarkred,layout=layout_in_circle(networkdarkred), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkdarkred)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkdarkred, directed = FALSE)
get_diameter(networkdarkred, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkdarkred)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkdarkred)

# Get the average path length of the graph g
g.apl <- mean_distance(networkdarkred, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an impordarkredt method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkdarkred)

gorder(networkdarkred)

g.random <- erdos.renyi.game(n = gorder(networkdarkred), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networkdarkred), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthdarkred.pdf")
hist(gl.apls,xlim=c(1.82,1.94),col="darkred",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.004, y= 200,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkdarkred, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenessdarkred.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="darkred", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkdarkred, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreedarkred.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="darkred", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkdarkred)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkdarkred <- make_ego_graph(networkdarkred, 1, "SubC.c_15222", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkdarkred)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]

V(egoNetworkdarkred)$color <- ifelse( names(V(egoNetworkdarkred)) %in% pOrfans$V1, "black", "darkred")

na.omit(pOrfans[match(names(V(egoNetworkdarkred)), pOrfans$V1),])

write.table(na.omit(pOrfans[match(names(V(egoNetworkdarkred)), pOrfans$V1),]),file = "darkredHighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_darkred.png", width = 4000, height = 4000);
plot(egoNetworkdarkred,layout=layout_nicely(egoNetworkdarkred), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#####################################

#turquoise module

probes = colnames(cpm.matrix)

turquoise = "turquoise"

inModule = is.finite(match(moduleColors, turquoise));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.160] = 0

(networkturquoise <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkturquoise <- simplify(networkturquoise,  degree(networkturquoise)==0))  # removes self-loops

((networkturquoise <- delete.vertices(networkturquoise, degree(networkturquoise)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkturquoise)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkturquoise)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_turquoise.png", width = 4000, height = 4000);
plot(networkturquoise,layout=layout_nicely(networkturquoise), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_turquoise.png", width = 4000, height = 4000);
plot(networkturquoise,layout=layout_with_fr(networkturquoise), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_turquoise.png", width = 4000, height = 4000);
plot(networkturquoise,layout=layout_as_tree(networkturquoise), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_turquoise.png", width = 4000, height = 4000);
plot(networkturquoise,layout=layout_in_circle(networkturquoise), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkturquoise)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkturquoise, directed = FALSE)
get_diameter(networkturquoise, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkturquoise)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkturquoise)

# Get the average path length of the graph g
g.apl <- mean_distance(networkturquoise, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an importurquoiset method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkturquoise)

gorder(networkturquoise)

g.random <- erdos.renyi.game(n = gorder(networkturquoise), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networkturquoise), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthturquoise.pdf")
hist(gl.apls,xlim=c(1.635,1.645),col="turquoise",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl +0.0003, y= 120,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkturquoise, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenessturquoise.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="turquoise", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkturquoise, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreeturquoise.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="turquoise", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkturquoise)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkturquoise <- make_ego_graph(networkturquoise, 1, "SubC.c_28986", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkturquoise)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworkturquoise)$color <- ifelse( names(V(egoNetworkturquoise)) %in% pOrfans$V1, "black", "turquoise")


na.omit(pOrfans[match(names(V(egoNetworkturquoise)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworkturquoise)), pOrfans$V1),]),file = "turquoiseHighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_turquoise.png", width = 4000, height = 4000);
plot(egoNetworkturquoise,layout=layout_nicely(egoNetworkturquoise), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#####################################

#orange module

probes = colnames(cpm.matrix)

orange = "orange"

inModule = is.finite(match(moduleColors, orange));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.053] = 0

(networkorange <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkorange <- simplify(networkorange,  degree(networkorange)==0))  # removes self-loops

((networkorange <- delete.vertices(networkorange, degree(networkorange)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkorange)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkorange)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_orange.png", width = 4000, height = 4000);
plot(networkorange,layout=layout_nicely(networkorange), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_orange.png", width = 4000, height = 4000);
plot(networkorange,layout=layout_with_fr(networkorange), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_orange.png", width = 4000, height = 4000);
plot(networkorange,layout=layout_as_tree(networkorange), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_orange.png", width = 4000, height = 4000);
plot(networkorange,layout=layout_in_circle(networkorange), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkorange)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkorange, directed = FALSE)
get_diameter(networkorange, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkorange)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkorange)

# Get the average path length of the graph g
g.apl <- mean_distance(networkorange, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an impororanget method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkorange)

gorder(networkorange)

g.random <- erdos.renyi.game(n = gorder(networkorange), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networkorange), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthorange.pdf")
hist(gl.apls,xlim=c(2.4,2.65),col="orange",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.005, y= 130,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkorange, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenessorange.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="orange", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkorange, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreeorange.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="orange", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkorange)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkorange <- make_ego_graph(networkorange, 1, "SubC.e_6180", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkorange)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworkorange)$color <- ifelse( names(V(egoNetworkorange)) %in% pOrfans$V1, "black", "orange")


na.omit(pOrfans[match(names(V(egoNetworkorange)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworkorange)), pOrfans$V1),]),file = "orangeHighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_orange.png", width = 4000, height = 4000);
plot(egoNetworkorange,layout=layout_nicely(egoNetworkorange), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#####################################

#greenyellow module

probes = colnames(cpm.matrix)

greenyellow = "greenyellow"

inModule = is.finite(match(moduleColors, greenyellow));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.127] = 0

(networkgreenyellow <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkgreenyellow <- simplify(networkgreenyellow,  degree(networkgreenyellow)==0))  # removes self-loops

((networkgreenyellow <- delete.vertices(networkgreenyellow, degree(networkgreenyellow)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkgreenyellow)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkgreenyellow)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_greenyellow.png", width = 4000, height = 4000);
plot(networkgreenyellow,layout=layout_nicely(networkgreenyellow), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_greenyellow.png", width = 4000, height = 4000);
plot(networkgreenyellow,layout=layout_with_fr(networkgreenyellow), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_greenyellow.png", width = 4000, height = 4000);
plot(networkgreenyellow,layout=layout_as_tree(networkgreenyellow), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_greenyellow.png", width = 4000, height = 4000);
plot(networkgreenyellow,layout=layout_in_circle(networkgreenyellow), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkgreenyellow)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkgreenyellow, directed = FALSE)
get_diameter(networkgreenyellow, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkgreenyellow)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkgreenyellow)

# Get the average path length of the graph g
g.apl <- mean_distance(networkgreenyellow, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an imporgreenyellowt method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkgreenyellow)

gorder(networkgreenyellow)

g.random <- erdos.renyi.game(n = gorder(networkgreenyellow), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networkgreenyellow), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthgreenyellow.pdf")
hist(gl.apls,xlim=c(1.78,1.92),col="greenyellow",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.003, y= 70,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkgreenyellow, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenessgreenyellow.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="greenyellow", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkgreenyellow, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreegreenyellow.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="greenyellow", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkgreenyellow)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkgreenyellow <- make_ego_graph(networkgreenyellow, 1, "SubC.c_3326", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkgreenyellow)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworkgreenyellow)$color <- ifelse( names(V(egoNetworkgreenyellow)) %in% pOrfans$V1, "black", "greenyellow")


na.omit(pOrfans[match(names(V(egoNetworkgreenyellow)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworkgreenyellow)), pOrfans$V1),]),file = "greenyellowHighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_greenyellow.png", width = 4000, height = 4000);
plot(egoNetworkgreenyellow,layout=layout_nicely(egoNetworkgreenyellow), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#####################################

#honeydew1 module

probes = colnames(cpm.matrix)

honeydew1 = "honeydew1"

inModule = is.finite(match(moduleColors, honeydew1));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.039] = 0

(networkhoneydew1 <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkhoneydew1 <- simplify(networkhoneydew1,  degree(networkhoneydew1)==0))  # removes self-loops

((networkhoneydew1 <- delete.vertices(networkhoneydew1, degree(networkhoneydew1)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkhoneydew1)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkhoneydew1)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_honeydew1.png", width = 4000, height = 4000);
plot(networkhoneydew1,layout=layout_nicely(networkhoneydew1), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_honeydew1.png", width = 4000, height = 4000);
plot(networkhoneydew1,layout=layout_with_fr(networkhoneydew1), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_honeydew1.png", width = 4000, height = 4000);
plot(networkhoneydew1,layout=layout_as_tree(networkhoneydew1), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_honeydew1.png", width = 4000, height = 4000);
plot(networkhoneydew1,layout=layout_in_circle(networkhoneydew1), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkhoneydew1)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkhoneydew1, directed = FALSE)
get_diameter(networkhoneydew1, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkhoneydew1)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkhoneydew1)

# Get the average path length of the graph g
g.apl <- mean_distance(networkhoneydew1, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an imporhoneydew1t method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkhoneydew1)

gorder(networkhoneydew1)

g.random <- erdos.renyi.game(n = gorder(networkhoneydew1), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networkhoneydew1), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthhoneydew1.pdf")
hist(gl.apls,xlim=c(1.6,1.8),col="honeydew1",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl +0.006, y= 90,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkhoneydew1, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenesshoneydew1.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="honeydew1", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkhoneydew1, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreehoneydew1.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="honeydew1", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkhoneydew1)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkhoneydew1 <- make_ego_graph(networkhoneydew1, 1, "SubC.c_5490", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkhoneydew1)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworkhoneydew1)$color <- ifelse( names(V(egoNetworkhoneydew1)) %in% pOrfans$V1, "black", "honeydew1")


na.omit(pOrfans[match(names(V(egoNetworkhoneydew1)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworkhoneydew1)), pOrfans$V1),]),file = "honeydew1HighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_honeydew1.png", width = 4000, height = 4000);
plot(egoNetworkhoneydew1,layout=layout_nicely(egoNetworkhoneydew1), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#####################################

#cyan module

probes = colnames(cpm.matrix)

cyan = "cyan"

inModule = is.finite(match(moduleColors, cyan));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.029] = 0

(networkcyan <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkcyan <- simplify(networkcyan,  degree(networkcyan)==0))  # removes self-loops

((networkcyan <- delete.vertices(networkcyan, degree(networkcyan)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkcyan)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkcyan)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_cyan.png", width = 4000, height = 4000);
plot(networkcyan,layout=layout_nicely(networkcyan), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_cyan.png", width = 4000, height = 4000);
plot(networkcyan,layout=layout_with_fr(networkcyan), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_cyan.png", width = 4000, height = 4000);
plot(networkcyan,layout=layout_as_tree(networkcyan), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_cyan.png", width = 4000, height = 4000);
plot(networkcyan,layout=layout_in_circle(networkcyan), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkcyan)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkcyan, directed = FALSE)
get_diameter(networkcyan, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkcyan)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkcyan)

# Get the average path length of the graph g
g.apl <- mean_distance(networkcyan, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an imporcyant method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkcyan)

gorder(networkcyan)

g.random <- erdos.renyi.game(n = gorder(networkcyan), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networkcyan), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthcyan.pdf")
hist(gl.apls,xlim=c(1.69,1.75),col="cyan",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl +0.0025, y= 140,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkcyan, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenesscyan.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="cyan", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkcyan, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreecyan.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="cyan", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkcyan)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkcyan <- make_ego_graph(networkcyan, 1, "SubC.c_9375", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkcyan)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworkcyan)$color <- ifelse( names(V(egoNetworkcyan)) %in% pOrfans$V1, "black", "cyan")


na.omit(pOrfans[match(names(V(egoNetworkcyan)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworkcyan)), pOrfans$V1),]),file = "cyanHighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_cyan.png", width = 4000, height = 4000);
plot(egoNetworkcyan,layout=layout_nicely(egoNetworkcyan), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#####################################

#lightsteelblue1 module

probes = colnames(cpm.matrix)

lightsteelblue1 = "lightsteelblue1"

inModule = is.finite(match(moduleColors, lightsteelblue1));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.040] = 0

(networklightsteelblue1 <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networklightsteelblue1 <- simplify(networklightsteelblue1,  degree(networklightsteelblue1)==0))  # removes self-loops

((networklightsteelblue1 <- delete.vertices(networklightsteelblue1, degree(networklightsteelblue1)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networklightsteelblue1)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networklightsteelblue1)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_lightsteelblue1.png", width = 4000, height = 4000);
plot(networklightsteelblue1,layout=layout_nicely(networklightsteelblue1), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_lightsteelblue1.png", width = 4000, height = 4000);
plot(networklightsteelblue1,layout=layout_with_fr(networklightsteelblue1), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_lightsteelblue1.png", width = 4000, height = 4000);
plot(networklightsteelblue1,layout=layout_as_tree(networklightsteelblue1), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_lightsteelblue1.png", width = 4000, height = 4000);
plot(networklightsteelblue1,layout=layout_in_circle(networklightsteelblue1), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networklightsteelblue1)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networklightsteelblue1, directed = FALSE)
get_diameter(networklightsteelblue1, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networklightsteelblue1)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networklightsteelblue1)

# Get the average path length of the graph g
g.apl <- mean_distance(networklightsteelblue1, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an imporlightsteelblue1t method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networklightsteelblue1)

gorder(networklightsteelblue1)

g.random <- erdos.renyi.game(n = gorder(networklightsteelblue1), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networklightsteelblue1), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthlightsteelblue1.pdf")
hist(gl.apls,xlim=c(1.65,1.75),col="lightsteelblue1",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl +0.0025, y= 170,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networklightsteelblue1, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenesslightsteelblue1.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="lightsteelblue1", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networklightsteelblue1, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreelightsteelblue1.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="lightsteelblue1", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networklightsteelblue1)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworklightsteelblue1 <- make_ego_graph(networklightsteelblue1, 1, "SubC.c_21423", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworklightsteelblue1)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworklightsteelblue1)$color <- ifelse( names(V(egoNetworklightsteelblue1)) %in% pOrfans$V1, "black", "lightsteelblue1")


na.omit(pOrfans[match(names(V(egoNetworklightsteelblue1)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworklightsteelblue1)), pOrfans$V1),]),file = "lightsteelblue1HighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_lightsteelblue1.png", width = 4000, height = 4000);
plot(egoNetworklightsteelblue1,layout=layout_nicely(egoNetworklightsteelblue1), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()


#####################################

#antiquewhite4 module

probes = colnames(cpm.matrix)

antiquewhite4 = "antiquewhite4"

inModule = is.finite(match(moduleColors, antiquewhite4));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.045] = 0

(networkantiquewhite4 <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkantiquewhite4 <- simplify(networkantiquewhite4,  degree(networkantiquewhite4)==0))  # removes self-loops

((networkantiquewhite4 <- delete.vertices(networkantiquewhite4, degree(networkantiquewhite4)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkantiquewhite4)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkantiquewhite4)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_antiquewhite4.png", width = 4000, height = 4000);
plot(networkantiquewhite4,layout=layout_nicely(networkantiquewhite4), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_antiquewhite4.png", width = 4000, height = 4000);
plot(networkantiquewhite4,layout=layout_with_fr(networkantiquewhite4), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_antiquewhite4.png", width = 4000, height = 4000);
plot(networkantiquewhite4,layout=layout_as_tree(networkantiquewhite4), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_antiquewhite4.png", width = 4000, height = 4000);
plot(networkantiquewhite4,layout=layout_in_circle(networkantiquewhite4), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkantiquewhite4)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkantiquewhite4, directed = FALSE)
get_diameter(networkantiquewhite4, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkantiquewhite4)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkantiquewhite4)

# Get the average path length of the graph g
g.apl <- mean_distance(networkantiquewhite4, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an imporantiquewhite4t method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkantiquewhite4)

gorder(networkantiquewhite4)

g.random <- erdos.renyi.game(n = gorder(networkantiquewhite4), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networkantiquewhite4), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthantiquewhite4.pdf")
hist(gl.apls,xlim=c(1.6,1.8),col="antiquewhite4",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl +0.0065, y= 140,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkantiquewhite4, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenessantiquewhite4.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="antiquewhite4", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkantiquewhite4, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreeantiquewhite4.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="antiquewhite4", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkantiquewhite4)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkantiquewhite4 <- make_ego_graph(networkantiquewhite4, 1, "SubC.c_16701", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkantiquewhite4)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworkantiquewhite4)$color <- ifelse( names(V(egoNetworkantiquewhite4)) %in% pOrfans$V1, "black", "antiquewhite4")


na.omit(pOrfans[match(names(V(egoNetworkantiquewhite4)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworkantiquewhite4)), pOrfans$V1),]),file = "antiquewhite4HighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_antiquewhite4.png", width = 4000, height = 4000);
plot(egoNetworkantiquewhite4,layout=layout_nicely(egoNetworkantiquewhite4), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#####################################

#violet module

probes = colnames(cpm.matrix)

violet = "violet"

inModule = is.finite(match(moduleColors, violet));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.073] = 0

(networkviolet <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkviolet <- simplify(networkviolet,  degree(networkviolet)==0))  # removes self-loops

((networkviolet <- delete.vertices(networkviolet, degree(networkviolet)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkviolet)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkviolet)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_violet.png", width = 4000, height = 4000);
plot(networkviolet,layout=layout_nicely(networkviolet), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_violet.png", width = 4000, height = 4000);
plot(networkviolet,layout=layout_with_fr(networkviolet), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_violet.png", width = 4000, height = 4000);
plot(networkviolet,layout=layout_as_tree(networkviolet), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_violet.png", width = 4000, height = 4000);
plot(networkviolet,layout=layout_in_circle(networkviolet), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkviolet)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkviolet, directed = FALSE)
get_diameter(networkviolet, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkviolet)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkviolet)

# Get the average path length of the graph g
g.apl <- mean_distance(networkviolet, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an imporviolett method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkviolet)

gorder(networkviolet)

g.random <- erdos.renyi.game(n = gorder(networkviolet), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networkviolet), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthviolet.pdf")
hist(gl.apls,xlim=c(1.5,1.56),col="violet",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl +0.0025, y= 150,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkviolet, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenessviolet.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="violet", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkviolet, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreeviolet.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="violet", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkviolet)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkviolet <- make_ego_graph(networkviolet, 1, "SubC.e_6351", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkviolet)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworkviolet)$color <- ifelse( names(V(egoNetworkviolet)) %in% pOrfans$V1, "black", "violet")


na.omit(pOrfans[match(names(V(egoNetworkviolet)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworkviolet)), pOrfans$V1),]),file = "violetHighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_violet.png", width = 4000, height = 4000);
plot(egoNetworkviolet,layout=layout_nicely(egoNetworkviolet), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#####################################

#lightpink4 module

probes = colnames(cpm.matrix)

lightpink4 = "lightpink4"

inModule = is.finite(match(moduleColors, lightpink4));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.065] = 0

(networklightpink4 <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networklightpink4 <- simplify(networklightpink4,  degree(networklightpink4)==0))  # removes self-loops

((networklightpink4 <- delete.vertices(networklightpink4, degree(networklightpink4)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networklightpink4)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networklightpink4)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_lightpink4.png", width = 4000, height = 4000);
plot(networklightpink4,layout=layout_nicely(networklightpink4), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_lightpink4.png", width = 4000, height = 4000);
plot(networklightpink4,layout=layout_with_fr(networklightpink4), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_lightpink4.png", width = 4000, height = 4000);
plot(networklightpink4,layout=layout_as_tree(networklightpink4), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_lightpink4.png", width = 4000, height = 4000);
plot(networklightpink4,layout=layout_in_circle(networklightpink4), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networklightpink4)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networklightpink4, directed = FALSE)
get_diameter(networklightpink4, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networklightpink4)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networklightpink4)

# Get the average path length of the graph g
g.apl <- mean_distance(networklightpink4, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an imporlightpink4t method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networklightpink4)

gorder(networklightpink4)

g.random <- erdos.renyi.game(n = gorder(networklightpink4), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
    gl[[i]] <- erdos.renyi.game(n = gorder(networklightpink4), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthlightpink4.pdf")
hist(gl.apls,xlim=c(1.74,1.9),col="lightpink4",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "green", lty = 3, lwd=2)
text(x=g.apl +0.0045, y= 100,labels="Calculated Average Path Length", col="green",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networklightpink4, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenesslightpink4.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="lightpink4", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networklightpink4, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreelightpink4.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="lightpink4", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networklightpink4)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworklightpink4 <- make_ego_graph(networklightpink4, 1, "SubC.c_9200", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworklightpink4)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworklightpink4)$color <- ifelse( names(V(egoNetworklightpink4)) %in% pOrfans$V1, "black", "lightpink4")


na.omit(pOrfans[match(names(V(egoNetworklightpink4)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworklightpink4)), pOrfans$V1),]),file = "lightpink4HighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_lightpink4.png", width = 4000, height = 4000);
plot(egoNetworklightpink4,layout=layout_nicely(egoNetworklightpink4), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

#####################################

#pink module

probes = colnames(cpm.matrix)

pink = "pink"

inModule = is.finite(match(moduleColors, pink));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.062] = 0

(networkpink <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkpink <- simplify(networkpink,  degree(networkpink)==0))  # removes self-loops

((networkpink <- delete.vertices(networkpink, degree(networkpink)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkpink)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkpink)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_pink.png", width = 4000, height = 4000);
plot(networkpink,layout=layout_nicely(networkpink), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_pink.png", width = 4000, height = 4000);
plot(networkpink,layout=layout_with_fr(networkpink), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_pink.png", width = 4000, height = 4000);
plot(networkpink,layout=layout_as_tree(networkpink), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_pink.png", width = 4000, height = 4000);
plot(networkpink,layout=layout_in_circle(networkpink), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkpink)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkpink, directed = FALSE)
get_diameter(networkpink, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkpink)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkpink)

# Get the average path length of the graph g
g.apl <- mean_distance(networkpink, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an imporpinkt method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkpink)

gorder(networkpink)

g.random <- erdos.renyi.game(n = gorder(networkpink), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 100 random graphs
gl <- vector('list', 100)

for(i in 1:100){
	print(i)
    gl[[i]] <- erdos.renyi.game(n = gorder(networkpink), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

summary(gl.apls)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthpink.pdf")
hist(gl.apls,xlim=c(1.665,1.69),col="pink",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.0015, y= 15,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkpink, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenesspink.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="pink", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkpink, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreepink.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="pink", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkpink)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkpink <- make_ego_graph(networkpink, 1, "SubC.e_20190", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkpink)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworkpink)$color <- ifelse( names(V(egoNetworkpink)) %in% pOrfans$V1, "black", "pink")


na.omit(pOrfans[match(names(V(egoNetworkpink)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworkpink)), pOrfans$V1),]),file = "pinkHighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_pink.png", width = 4000, height = 4000);
plot(egoNetworkpink,layout=layout_nicely(egoNetworkpink), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

###############################################

#brown4 module

probes = colnames(cpm.matrix)

brown4 = "brown4"

inModule = is.finite(match(moduleColors, brown4));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.047] = 0

(networkbrown4 <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkbrown4 <- simplify(networkbrown4,  degree(networkbrown4)==0))  # removes self-loops

((networkbrown4 <- delete.vertices(networkbrown4, degree(networkbrown4)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkbrown4)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkbrown4)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_brown4.png", width = 4000, height = 4000);
plot(networkbrown4,layout=layout_nicely(networkbrown4), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_brown4.png", width = 4000, height = 4000);
plot(networkbrown4,layout=layout_with_fr(networkbrown4), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_brown4.png", width = 4000, height = 4000);
plot(networkbrown4,layout=layout_as_tree(networkbrown4), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_brown4.png", width = 4000, height = 4000);
plot(networkbrown4,layout=layout_in_circle(networkbrown4), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkbrown4)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkbrown4, directed = FALSE)
get_diameter(networkbrown4, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkbrown4)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkbrown4)

# Get the average path length of the graph g
g.apl <- mean_distance(networkbrown4, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an imporbrown4t method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkbrown4)

gorder(networkbrown4)

g.random <- erdos.renyi.game(n = gorder(networkbrown4), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
	print(i)
    gl[[i]] <- erdos.renyi.game(n = gorder(networkbrown4), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

summary(gl.apls)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthbrown4.pdf")
hist(gl.apls,xlim=c(2,2.3),col="brown4",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.01, y= 135,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkbrown4, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenessbrown4.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="brown4", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkbrown4, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreebrown4.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="brown4", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkbrown4)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkbrown4 <- make_ego_graph(networkbrown4, 1, "SubC.e_581", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkbrown4)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworkbrown4)$color <- ifelse( names(V(egoNetworkbrown4)) %in% pOrfans$V1, "black", "brown4")


na.omit(pOrfans[match(names(V(egoNetworkbrown4)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworkbrown4)), pOrfans$V1),]),file = "brown4HighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_brown4.png", width = 4000, height = 4000);
plot(egoNetworkbrown4,layout=layout_nicely(egoNetworkbrown4), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()


#lavenderblush3 module

probes = colnames(cpm.matrix)

lavenderblush3 = "lavenderblush3"

inModule = is.finite(match(moduleColors, lavenderblush3));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.052] = 0

(networklavenderblush3 <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networklavenderblush3 <- simplify(networklavenderblush3,  degree(networklavenderblush3)==0))  # removes self-loops

((networklavenderblush3 <- delete.vertices(networklavenderblush3, degree(networklavenderblush3)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networklavenderblush3)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networklavenderblush3)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_lavenderblush3.png", width = 4000, height = 4000);
plot(networklavenderblush3,layout=layout_nicely(networklavenderblush3), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_lavenderblush3.png", width = 4000, height = 4000);
plot(networklavenderblush3,layout=layout_with_fr(networklavenderblush3), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_lavenderblush3.png", width = 4000, height = 4000);
plot(networklavenderblush3,layout=layout_as_tree(networklavenderblush3), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_lavenderblush3.png", width = 4000, height = 4000);
plot(networklavenderblush3,layout=layout_in_circle(networklavenderblush3), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networklavenderblush3)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networklavenderblush3, directed = FALSE)
get_diameter(networklavenderblush3, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networklavenderblush3)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networklavenderblush3)

# Get the average path length of the graph g
g.apl <- mean_distance(networklavenderblush3, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an imporlavenderblush3t method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networklavenderblush3)

gorder(networklavenderblush3)

g.random <- erdos.renyi.game(n = gorder(networklavenderblush3), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
	print(i)
    gl[[i]] <- erdos.renyi.game(n = gorder(networklavenderblush3), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

summary(gl.apls)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthlavenderblush3.pdf")
hist(gl.apls,xlim=c(1.8,2.1),col="lavenderblush3",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.01, y= 150,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networklavenderblush3, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenesslavenderblush3.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="lavenderblush3", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networklavenderblush3, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreelavenderblush3.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="lavenderblush3", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networklavenderblush3)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworklavenderblush3 <- make_ego_graph(networklavenderblush3, 1, "SubC.c_10098", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworklavenderblush3)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworklavenderblush3)$color <- ifelse( names(V(egoNetworklavenderblush3)) %in% pOrfans$V1, "black", "lavenderblush3")


na.omit(pOrfans[match(names(V(egoNetworklavenderblush3)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworklavenderblush3)), pOrfans$V1),]),file = "lavenderblush3HighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_lavenderblush3.png", width = 4000, height = 4000);
plot(egoNetworklavenderblush3,layout=layout_nicely(egoNetworklavenderblush3), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()


#salmon4 module

probes = colnames(cpm.matrix)

salmon4 = "salmon4"

inModule = is.finite(match(moduleColors, salmon4));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.086] = 0

(networksalmon4 <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networksalmon4 <- simplify(networksalmon4,  degree(networksalmon4)==0))  # removes self-loops

((networksalmon4 <- delete.vertices(networksalmon4, degree(networksalmon4)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networksalmon4)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networksalmon4)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_salmon4.png", width = 4000, height = 4000);
plot(networksalmon4,layout=layout_nicely(networksalmon4), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_salmon4.png", width = 4000, height = 4000);
plot(networksalmon4,layout=layout_with_fr(networksalmon4), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_salmon4.png", width = 4000, height = 4000);
plot(networksalmon4,layout=layout_as_tree(networksalmon4), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_salmon4.png", width = 4000, height = 4000);
plot(networksalmon4,layout=layout_in_circle(networksalmon4), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networksalmon4)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networksalmon4, directed = FALSE)
get_diameter(networksalmon4, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networksalmon4)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networksalmon4)

# Get the average path length of the graph g
g.apl <- mean_distance(networksalmon4, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an imporsalmon4t method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networksalmon4)

gorder(networksalmon4)

g.random <- erdos.renyi.game(n = gorder(networksalmon4), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
	print(i)
    gl[[i]] <- erdos.renyi.game(n = gorder(networksalmon4), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

summary(gl.apls)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthsalmon4.pdf")
hist(gl.apls,xlim=c(1.50,1.6),col="salmon4",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.003, y= 90,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networksalmon4, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenesssalmon4.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="salmon4", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networksalmon4, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreesalmon4.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="salmon4", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networksalmon4)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworksalmon4 <- make_ego_graph(networksalmon4, 1, "SubC.c_8356", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworksalmon4)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworksalmon4)$color <- ifelse( names(V(egoNetworksalmon4)) %in% pOrfans$V1, "black", "salmon4")


na.omit(pOrfans[match(names(V(egoNetworksalmon4)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworksalmon4)), pOrfans$V1),]),file = "salmon4HighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_salmon4.png", width = 4000, height = 4000);
plot(egoNetworksalmon4,layout=layout_nicely(egoNetworksalmon4), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()


#darkseagreen4 module

probes = colnames(cpm.matrix)

darkseagreen4 = "darkseagreen4"

inModule = is.finite(match(moduleColors, darkseagreen4));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.042] = 0

(networkdarkseagreen4 <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkdarkseagreen4 <- simplify(networkdarkseagreen4,  degree(networkdarkseagreen4)==0))  # removes self-loops

((networkdarkseagreen4 <- delete.vertices(networkdarkseagreen4, degree(networkdarkseagreen4)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkdarkseagreen4)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkdarkseagreen4)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_darkseagreen4.png", width = 4000, height = 4000);
plot(networkdarkseagreen4,layout=layout_nicely(networkdarkseagreen4), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_darkseagreen4.png", width = 4000, height = 4000);
plot(networkdarkseagreen4,layout=layout_with_fr(networkdarkseagreen4), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_darkseagreen4.png", width = 4000, height = 4000);
plot(networkdarkseagreen4,layout=layout_as_tree(networkdarkseagreen4), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_darkseagreen4.png", width = 4000, height = 4000);
plot(networkdarkseagreen4,layout=layout_in_circle(networkdarkseagreen4), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkdarkseagreen4)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkdarkseagreen4, directed = FALSE)
get_diameter(networkdarkseagreen4, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkdarkseagreen4)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkdarkseagreen4)

# Get the average path length of the graph g
g.apl <- mean_distance(networkdarkseagreen4, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an impordarkseagreen4t method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkdarkseagreen4)

gorder(networkdarkseagreen4)

g.random <- erdos.renyi.game(n = gorder(networkdarkseagreen4), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
	print(i)
    gl[[i]] <- erdos.renyi.game(n = gorder(networkdarkseagreen4), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

summary(gl.apls)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthdarkseagreen4.pdf")
hist(gl.apls,xlim=c(1.65,1.76),col="darkseagreen4",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.003, y= 175,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkdarkseagreen4, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenessdarkseagreen4.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="darkseagreen4", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkdarkseagreen4, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreedarkseagreen4.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="darkseagreen4", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkdarkseagreen4)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkdarkseagreen4 <- make_ego_graph(networkdarkseagreen4, 1, "SubC.c_7265", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkdarkseagreen4)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworkdarkseagreen4)$color <- ifelse( names(V(egoNetworkdarkseagreen4)) %in% pOrfans$V1, "black", "darkseagreen4")


na.omit(pOrfans[match(names(V(egoNetworkdarkseagreen4)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworkdarkseagreen4)), pOrfans$V1),]),file = "darkseagreen4HighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_darkseagreen4.png", width = 4000, height = 4000);
plot(egoNetworkdarkseagreen4,layout=layout_nicely(egoNetworkdarkseagreen4), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()


#darkorange module

probes = colnames(cpm.matrix)

darkorange = "darkorange"

inModule = is.finite(match(moduleColors, darkorange));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.037] = 0

(networkdarkorange <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkdarkorange <- simplify(networkdarkorange,  degree(networkdarkorange)==0))  # removes self-loops

((networkdarkorange <- delete.vertices(networkdarkorange, degree(networkdarkorange)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkdarkorange)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkdarkorange)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_darkorange.png", width = 4000, height = 4000);
plot(networkdarkorange,layout=layout_nicely(networkdarkorange), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_darkorange.png", width = 4000, height = 4000);
plot(networkdarkorange,layout=layout_with_fr(networkdarkorange), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_darkorange.png", width = 4000, height = 4000);
plot(networkdarkorange,layout=layout_as_tree(networkdarkorange), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_darkorange.png", width = 4000, height = 4000);
plot(networkdarkorange,layout=layout_in_circle(networkdarkorange), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkdarkorange)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkdarkorange, directed = FALSE)
get_diameter(networkdarkorange, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkdarkorange)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkdarkorange)

# Get the average path length of the graph g
g.apl <- mean_distance(networkdarkorange, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an impordarkoranget method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkdarkorange)

gorder(networkdarkorange)

g.random <- erdos.renyi.game(n = gorder(networkdarkorange), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
	print(i)
    gl[[i]] <- erdos.renyi.game(n = gorder(networkdarkorange), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

summary(gl.apls)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthdarkorange.pdf")
hist(gl.apls,xlim=c(1.68,1.73),col="darkorange",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.001, y= 150,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkdarkorange, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenessdarkorange.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="darkorange", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkdarkorange, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreedarkorange.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="darkorange", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkdarkorange)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkdarkorange <- make_ego_graph(networkdarkorange, 1, "SubC.e_19195", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkdarkorange)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworkdarkorange)$color <- ifelse( names(V(egoNetworkdarkorange)) %in% pOrfans$V1, "black", "darkorange")


na.omit(pOrfans[match(names(V(egoNetworkdarkorange)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworkdarkorange)), pOrfans$V1),]),file = "darkorangeHighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_darkorange.png", width = 4000, height = 4000);
plot(egoNetworkdarkorange,layout=layout_nicely(egoNetworkdarkorange), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()


#grey60 module

probes = colnames(cpm.matrix)

grey60 = "grey60"

inModule = is.finite(match(moduleColors, grey60));
modProbes <- probes[inModule]

adj <- TOM[inModule,inModule]
dimnames(adj) <- list(modProbes,modProbes)

round(quantile(adj[adj != 1], probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.96,.97,.98,.99,.999,1)),3)

adj[adj < 0.033] = 0

(networkgrey60 <- graph.adjacency(adj, mode="undirected", weighted = TRUE))

(networkgrey60 <- simplify(networkgrey60,  degree(networkgrey60)==0))  # removes self-loops

((networkgrey60 <- delete.vertices(networkgrey60, degree(networkgrey60)==0)))

ColorTable <- cbind(probes,mergedColors)

V(networkgrey60)$color <- ColorTable[is.finite(match(ColorTable[,1], names(V(networkgrey60)))),2]

#vertex.size = sqrt(g.b)+1
png(file = "nicely.network_grey60.png", width = 4000, height = 4000);
plot(networkgrey60,layout=layout_nicely(networkgrey60), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

png(file = "fr.network_grey60.png", width = 4000, height = 4000);
plot(networkgrey60,layout=layout_with_fr(networkgrey60), edge.arrow.size = 0.05)
dev.off()

#layout_as_tree

png(file = "tree.network_grey60.png", width = 4000, height = 4000);
plot(networkgrey60,layout=layout_as_tree(networkgrey60), edge.arrow.size = 0.05)
dev.off()

png(file = "circle.network_grey60.png", width = 4000, height = 4000);
plot(networkgrey60,layout=layout_in_circle(networkgrey60), edge.arrow.size = 0.05)
dev.off()

# Get density of a graph - This is essentially the proportion of all potential edges between vertices that actually exist in the network graph.
gd <- edge_density(networkgrey60)

#Another measure of how interconnected a network is average path length. This is calculated by determining the mean of the lengths of the shortest paths between all pairs of vertices in the network. The longest path length between any pair of vertices is called the diameter of the network graph.

# Get the diameter of the graph g
diameter(networkgrey60, directed = FALSE)
get_diameter(networkgrey60, directed = FALSE)

# Which two vertices are the furthest apart in the graph ?
farthest_vertices(networkgrey60)
# Shows the path sequence between two furthest apart vertices.
get_diameter(networkgrey60)

# Get the average path length of the graph g
g.apl <- mean_distance(networkgrey60, directed = FALSE)
g.apl

#Random graphs
#Generating random graphs is an imporgrey60t method for investigating how likely or unlikely other network metrics are likely to occur given certain properties of the original graph. The simplest random graph is one that has the same number of vertices as your original graph and approximately the same density as the original graph Generate a random graph using the function erdos.renyi.game(). The first argument n should be the number of nodes of the graph g which can be calculated using gorder(), the second argument p.or.m should be the density of the graph g which you previously stored as the object gd. The final argument is set as type='gnp' to tell the function that you are using the density of the graph to generate a random graph. Store this new graph as the vector g.random

# Get density of a graph
gd <- edge_density(networkgrey60)

gorder(networkgrey60)

g.random <- erdos.renyi.game(n = gorder(networkgrey60), p.or.m = gd, type = "gnp")

# Get density of new random graph `g.random`
edge_density(g.random)

#Get the average path length of the random graph g.random

mean_distance(g.random, directed = FALSE)

#Generate 1000 random graphs of the original graph g by executing the code that creates the list object gl and the for loop.

# Generate 1000 random graphs
gl <- vector('list', 1000)

for(i in 1:1000){
	print(i)
    gl[[i]] <- erdos.renyi.game(n = gorder(networkgrey60), p.or.m = gd, type = "gnp")
}

#Calculate the average path length of the 1000 random graphs using lapply(). Create a vector gl.apls of these 1000 values by executing the code that uses unlist()

# Calculate average path length of 1000 random graphs
gl.apl <- lapply(gl, mean_distance, directed = FALSE)
gl.apls <- unlist(gl.apl)

summary(gl.apls)

#Plot a histogram of the average path lengths of the 1000 random graphs using hist() on the vector gl.apls. Add a red dashed line to the plot using abline() with the x-intercept being the value of the average path length of the original graph g.apl. You calculated this value in the previous exercise.

# Plot the distribution of average path lengths
pdf("distributionOfAveragePathLengthgrey60.pdf")
hist(gl.apls,xlim=c(1.695,1.71),col="grey60",main="Theoretical Distribution Vs Calculated Average Path Lenth",xlab="Average Path Length")
abline(v = g.apl, col = "red", lty = 3, lwd=2)
text(x=g.apl -0.0005, y= 95,labels="Calculated Average Path Length", col="red",srt=90)
axis(1, at=g.apl,col="red",labels=F)

mtext(1,at=g.apl,text=round(g.apl,3),col="red",line=.3,cex=.85)

dev.off()

#Calculate the proportion of times that the values of the average path length of random graphs gl.apls are lower than the value of the original graph g.apl. This is essentially the probability that we would expect our observed average path length by chance given the original density and number of vertices of the original graph.

# Calculate the proportion of graphs with an average path length lower than our observed
sum(gl.apls < g.apl)/1000

#Betweenness
#Another measure of the importance of a given vertex is its betweenness. This is an index of how frequently the vertex lies on shortest paths between any two vertices in the network.

g.b <- betweenness(networkgrey60, directed = FALSE)

which.max(g.b)

g.b[which.max(g.b)]

# Show histogram of vertex betweenness
pdf("histVertexBetweenessgrey60.pdf")
hist(g.b,breaks = g.b[which.max(g.b)], col ="grey60", main="Betweenness Histogram",xlab="Betweenness")
dev.off()

#Distances between vertices

#The inter-connectivity of a network can be assessed by examining the number and length of paths between vertices. A path is simply the chain of connections between vertices. The number of intervening edges between two vertices represents the geodesic distance between vertices. Vertices that are connected to each other have a geodesic distance of 1. Those that share a neighbor in common but are not connected to each other have a geodesic distance of 2 and so on.

#Perhaps the most straightforward measure of vertex importance is the degree of a vertex. The out-degree of a vertex is the number of other individuals to which a vertex has an outgoing edge directed to.

g.outd <- degree(networkgrey60, mode = c("all"))

#Make a histogram
# Find the vertex that has the maximum out-degree
which.max(g.outd)

g.outd[which.max(g.outd)]

pdf ("histDegreegrey60.pdf")
hist(g.outd, breaks = g.outd[which.max(g.outd)]
, col ="grey60", main="Degree Histogram",xlab="Degree")
dev.off()

#Eigenvector centrality scores correspond to the values of the first eigenvector of the graph adjacency matrix; these scores may, in turn, be interpreted as arising from a reciprocal process in which the centrality of each actor is proportional to the sum of  the centralities of those actors to whom he or she is connected. In general, vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).  (The perceptive may realize that this implies that the largest values will be obtained by   individuals in large cliques (or high-density substructures).   This is also intelligible from an algebraic point of view, with  the first eigenvector being closely related to the best rank-1 approximation of the adjacency matrix.

# Identify key nodes using eigenvector centrality

g.ec <- eigen_centrality(networkgrey60)
g.ec$vector[which.max(g.ec$vector)]

#Use make_ego_graph() to create a subset of our network comprised of vertices that are connected to vertex most connected vertex.

#Identify vertices that are reachable within one connection

#These functions find the vertices not farther than a given limit from another fixed vertex, these are called the neighborhood of the vertex.

egoNetworkgrey60 <- make_ego_graph(networkgrey60, 1, "SubC.e_9654", mode=c("all"))[[1]]

#selected the putative orphan genes closed to the egengene

temp <- descriptions[is.finite(match(descriptions$V1,names(V(egoNetworkgrey60)))),]

pOrfans <- temp[grep("---NA---|uncharacterized protein", temp$V2),]


V(egoNetworkgrey60)$color <- ifelse( names(V(egoNetworkgrey60)) %in% pOrfans$V1, "black", "grey60")


na.omit(pOrfans[match(names(V(egoNetworkgrey60)), pOrfans$V1),])


write.table(na.omit(pOrfans[match(names(V(egoNetworkgrey60)), pOrfans$V1),]),file = "grey60HighCentralityOrphans.tab",sep ="\t")

png(file = "nicely.network_ego_grey60.png", width = 4000, height = 4000);
plot(egoNetworkgrey60,layout=layout_nicely(egoNetworkgrey60), edge.arrow.size =0.05,vertex.size=3,vertex.label = NA)
dev.off()

###############################################
#####Expression Barplots ######################
###############################################
GeneModules <- cbind(probes,moduleColors)
GeneModules <- cbind(GeneModules, moduleColors)
colnames(GeneModules) <- c("probes","moduleColors", "displayColor")

table(GeneModules[,3])

GeneModules[,3] <- gsub("antiquewhite4","antiquewhite",GeneModules[,3])

GeneModules[,3] <- gsub("honeydew1","honeydew",GeneModules[,3])

GeneModules[,3] <- gsub("thistle1","thistle",GeneModules[,3])

GeneModules[,3] <- gsub("brown4","brown",GeneModules[,3])

GeneModules[,3] <- gsub("darkseagreen4","darkseagreen",GeneModules[,3])

GeneModules[,3] <- gsub("lavenderblush3","lavenderblush",GeneModules[,3])

GeneModules[,3] <- gsub("plum2","plum",GeneModules[,3])

GeneModules[,3] <- gsub("coral2","coral",GeneModules[,3])

GeneModules[,3] <- gsub("lightpink4","lightpink",GeneModules[,3])

GeneModules[,3] <- gsub("salmon4","salmon",GeneModules[,3])

GeneModules[,3] <- gsub("lightsteelblue1","lightsteelblue",GeneModules[,3])

GeneModules[,3] <- gsub("sienna3","sienna",GeneModules[,3])

GeneModules[,3] <- gsub("grey60","grey",GeneModules[,3])



Expression.matrix <- cpm(analysis_matrix,normalized.lib.sizes=F)

#"#7bc437"
largura = 190/25.4
altura = largura/1.618

i=10
fontColor = "black"
system("mkdir PointExpression")
for (gene in rownames(Expression.matrix)){
print(gene)
limitOfY <- max(Expression.matrix[gene,])*1.2
pdf(paste0("./PointExpression/exp_",GeneModules[match(gene, GeneModules[,1]),3],"_",gene,".pdf"),h= altura,w= largura)
plot(Expression.matrix[gene,],col="black", yaxt='n',,xaxt='n',ylim=c(0,limitOfY),main =paste("Expression of",gene),cex.main=1.5,cex =1.2,type="p", ylab="",xlab ="Leaf RNAseq libraries")

title(ylab="Counts Per Million",cex.lab=1.5,line=2.3)

if( limitOfY <= 50 ){
i=3
}
if( limitOfY > 50 & limitOfY <= 100 ){
i=10
}
if( limitOfY > 100 & limitOfY <= 1000 ){
i=50
}
if( limitOfY > 1000 & limitOfY <= 10000 ){
i=100
}
if( limitOfY > 10000 & limitOfY <= 100000 ){
i=1000
}
if( limitOfY > 100000 & limitOfY <= 1000000 ){
i=10000
}
if( limitOfY > 1000000 & limitOfY <= 10000000 ){
i=100000
}

axis(2,seq(0,limitOfY,i),labels=F)
mtext(side=2,text=seq(0,limitOfY,i),outer=F,las=2,line=.8,at=seq(0,limitOfY,i),cex=.8)

mtext(1, at= seq(1,24,1), text=colnames(Expression.matrix),cex=0.7,las=2, line = 0.5)

mtext(3, at = 12, text = descriptions[match(gene,descriptions[,1]),2], cex = 1, las = 1, line = 0.2)

rect(xleft = -2, ybottom =  (limitOfY + limitOfY*0.15), xright = 4, ytop =(limitOfY + limitOfY*0.3), density = NA, col =GeneModules[match(gene, GeneModules[,1]),2], xpd=T)

#GeneModules[match(gene, GeneModules[,1]),2] == "darkorange" |

if (GeneModules[match(gene, GeneModules[,1]),2] == "antiquewhite4" |  GeneModules[match(gene, GeneModules[,1]),2] == "lightpink" | GeneModules[match(gene, GeneModules[,1]),2] == "darkslateblue" | GeneModules[match(gene, GeneModules[,1]),2] == "darkred" | GeneModules[match(gene, GeneModules[,1]),2] == "brown4" | GeneModules[match(gene, GeneModules[,1]),2] == "darkseagreen4" ){
fontColor = "white"
}
else{
	fontColor = "black"
}
text(x=0.52, y=(limitOfY + limitOfY*0.22), labels = paste(GeneModules[match(gene, GeneModules[,1]),3]),xpd=T , cex= 1.25, col=fontColor)

dev.off()

GeneModules[,1] == rownames(Expression.matrix)

rowMeans(Expression.matrix)

GeneModules <- cbind(GeneModules, rowMeans(Expression.matrix))

df.GeneModules <- as.data.frame(GeneModules)

df.GeneModules$V4 <- as.numeric(df.GeneModules$V4)

colnames(df.GeneModules) <- c("probes","moduleColors", "displayColor", "meanCPM")

write.table(df.GeneModules, file = "GeneModules.tab",sep="\t")



######################
####Calculate coverage

#' Calculate gene density on chromosome
#'
#' @description Calculate density of genes, CDS, mRNAs, tRNAs, etc on chromosome from a given GFF file with a slide window
#'
#' @param gff_file GFF file
#' @param karyotype_info tab-separated file (or a data frame) containing karyotype information, with first column being "Chr", second being "Start", third being "End". File should NOT have header.
#' @param feature interested feature to calculate, should be in GFF file
#' @param window length of each slide window, defaut 100kb
#'
#' @importFrom tidyr separate
#' @author Yujie Liu
#' @export
#'


geneDensity <- calcDensityGFF("GeneCoordnates.tab", "chromossomes.tab",feature = "gene", window = 1e5)

calcDensityGFF <- function(gff_file,
                           karyotype_info,
                           feature = "gene",
                           window = 1e5) {
  # read in files
  gff <- read.delim(gff_file, header = FALSE, comment.char = "#")
  if (is.character(karyotype_info)) {
    karyotype <- read.delim(karyotype_info, header = FALSE)
  } else {
    karyotype <- karyotype_info
  }

  # select only interested feature
  all_chrs <- karyotype[[1]]
  gff <- gff[gff[[1]] %in% all_chrs & gff[[3]] == feature, ]

  # make output data frame
  density <- data.frame()

  # loop against each chromosome
  for (chr in all_chrs) {
    chr_end <- karyotype[karyotype[[1]] == chr, 3]

    # count, but only consider the start of each gene
    tmp <-
      data.frame(table(cut(gff[gff[[1]] == chr, 4],
                           breaks = c(
                             seq(0, chr_end, window),
                             chr_end
                           ))))

    # make outcome tidy
    tmp <-
      tidyr::separate(tmp,
                      1,
                      into = c("Start", "End"),
                      sep = ",")
    tmp$Start <-
      as.numeric(gsub("\\(", "", tmp$Start)) + 1
    tmp$End <-
      as.numeric(gsub("\\]", "", tmp$End))

    # add chr info
    tmp$Chr <- chr

    # reorder columns
    tmp <- tmp[c(4, 1:3)]
    colnames(tmp) <- c("Chr", "Start", "End", "Count")

    # modify last end_pos to chr_end
    tmp[nrow(tmp), "End"] <- chr_end

    # combine
    density <- rbind(density, tmp)
  }

  return(density)
}

############################
#####Circle plots###########
############################


#Unify old IDS with new ones - done - thank you past me

"processedC.arabicaGeneDescriptions.tab"

#Export a table with IDs, module, real color and mean CPM

options ( stringsAsFactors = TRUE) ;
library(OmicCircos)

X <- read.delim("chromossomes.tab",sep = "\t",header=F)

AllGenes <- read.delim("AllGenes.mapping.tab",sep="\t",header=F)

geneDescriptions <- read.delim("processedC.arabicaGeneDescriptions.tab", sep="\t",header=F)

AllGenes <- cbind(AllGenes,geneDescriptions[,c(2,3)])

colnames(AllGenes) <- c("chr","pos","uniID","subGenID","description")

colnames(X) <- c("seg.name", "seg.Start", "seg.End", "name", "color")

GeneModules <- read.delim("GeneModules.tab",sep="\t",header=T)

exp.leaf.GeneModules <- AllGenes[(AllGenes$subGenID %in% GeneModules$probes),]

exp.leaf.GeneModules <- cbind(exp.leaf.GeneModules, GeneModules[GeneModules$probes %in% exp.leaf.GeneModules$subGenID,c(2,4)] )

exp.leaf.GeneModules <- cbind(exp.leaf.GeneModules, log(exp.leaf.GeneModules[,7],2))

db <- segAnglePo(seg.dat=X, seg= X$seg.name) ;

colors <- as.vector(X$color)

mapX <- as.data.frame(cbind(X$seg.name,as.integer(X$seg.End),as.character(X$name)))

par(mar=c(2,2,2,2))

plot(c(1,800),c(1,800),type="n", axes=FALSE,xlab="",ylab="",main="")

#chromossome names
adjust = 20

circos(R=425-adjust, cir=db, W=0, mapping= mapX, type="label2",side="out", col="black" , cex=0.5);

circos(R=400-adjust,cir=db,type="chr",col= colors,print.chr.lab=FALSE,W=40,scale=TRUE)

#All putative genes
circos(R=370-adjust, cir=db, W=20, mapping=AllGenes[grep("---NA---", AllGenes$description, invert=T),] , type="s2", col="black",B=TRUE,cex =.004);

#NOVEL genes
circos(R=350-adjust, cir=db, W=20, mapping=AllGenes[grep("---NA---", AllGenes$description),] , type="s2", col="red",B=TRUE,cex =.04);

#Verified expression in fully expanded leaves

circos(R=325-adjust, cir=db, W=20, mapping=exp.leaf.GeneModules , type="s2", col="darkseagreen",B=TRUE,cex =.04);

circos(R=300-adjust, cir=db, W=20, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "pink",] , type="s2", col="pink",B=TRUE,cex =.04);

circos(R=275-adjust, cir=db, W=20, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "grey60",] , type="s2", col="grey60",B=TRUE,cex =.04);

circos(R=250-adjust, cir=db, W=20, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "turquoise",] , type="s2", col="turquoise",B=TRUE,cex =.04);

circos(R=225-adjust, cir=db, W=20, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "darkred",] , type="s2", col="darkred",B=TRUE,cex =.04);

circos(R=200-adjust, cir=db, W=20, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "darkorange",] , type="s2", col="darkorange",B=TRUE,cex =.04);

circos(R=175-adjust, cir=db, W=20, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "cyan",] , type="s2", col="cyan",B=TRUE,cex =.04);

circos(R=150-adjust, cir=db, W=20, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "lightsteelblue1",] , type="s2", col="lightsteelblue1",B=TRUE,cex =.04);

circos(R=125-adjust, cir=db, W=20, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "greenyellow",] , type="s2", col="greenyellow",B=TRUE,cex =.04);

circos(R=100-adjust, cir=db, W=20, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "tan",] , type="s2", col="tan",B=TRUE,cex =.04);

circos(R=75-adjust, cir=db, W=20, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "orange",] , type="s2", col="orange",B=TRUE,cex =.04);

circos(R=50-adjust, cir=db, W=20, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "violet",] , type="s2", col="violet",B=TRUE,cex =.04);

#circos(R=25-adjust, cir=db, W=20, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "sienna3",] , type="s2", col="sienna3",B=TRUE,cex =.04);

###########################################################
######Only general data of gen position

par(mar=c(2,2,2,2))

plot(c(1,800),c(1,800),type="n", axes=FALSE,xlab="",ylab="",main="")

#chromossome names
adjust = 20

circos(R=425-adjust, cir=db, W=0, mapping= mapX, type="label2",side="out", col="grey20" , cex=0.5);

circos(R=400-adjust,cir=db,type="chr",col= colors,print.chr.lab=FALSE,W=40,scale=TRUE)

#All putative genes
circos(R=340-adjust, cir=db, W=50, mapping=density,col.v=4 , type="h", col="black",B=TRUE,cex =.1);

#NOVEL genes
circos(R=290-adjust, cir=db, W=50, mapping=AllGenes[grep("---NA---", AllGenes$description),] , type="s2", col="#976464",B=TRUE,cex =.1);

#Verified expression in fully expanded leaves

circos(R=240-adjust, cir=db, W=50, mapping=exp.leaf.GeneModules , type="s2", col="#5F6D54",B=TRUE,cex =.1);

circos(R=190-adjust, cir=db, W=50, mapping=exp.leaf.GeneModules,col.v=8 , type="b", col="#C6CA97",B=TRUE,cex =.1);

######################
#Group specific plots#
######################

#Group A
par(mar=c(2,2,2,2))

plot(c(1,800),c(1,800),type="n", axes=FALSE,xlab="",ylab="",main="")

#chromossome names
adjust = 20

circos(R=425-adjust, cir=db, W=0, mapping= mapX, type="label2",side="out", col="grey20" , cex=0.5);

circos(R=400-adjust,cir=db,type="chr",col= colors,print.chr.lab=FALSE,W=40,scale=TRUE)

#All putative genes
circos(R=340-adjust, cir=db, W=50, mapping=density,col.v=4 , type="h", col="black",B=TRUE,cex =.1);

circos(R=310-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "plum2",] , type="s", col="plum2",B=TRUE,cex =.8);

circos(R=280-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "darkslateblue",] , type="s", col="darkslateblue",B=TRUE,cex =.8);

circos(R=250-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "thistle1",] , type="s", col="black",B=TRUE,cex =.8);

circos(R=190-adjust, cir=db, W=50, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "thistle1" | exp.leaf.GeneModules$moduleColors == "darkslateblue" | exp.leaf.GeneModules$moduleColors == "plum2",],col.v=8 , type="b", col="darkslateblue",B=TRUE,cex =.1);


######################

#Group B
par(mar=c(2,2,2,2))

plot(c(1,800),c(1,800),type="n", axes=FALSE,xlab="",ylab="",main="")

#chromossome names
adjust = 20

circos(R=425-adjust, cir=db, W=0, mapping= mapX, type="label2",side="out", col="grey20" , cex=0.5);

circos(R=400-adjust,cir=db,type="chr",col= colors,print.chr.lab=FALSE,W=40,scale=TRUE)

#All putative genes
circos(R=340-adjust, cir=db, W=50, mapping=density,col.v=4 , type="h", col="black",B=TRUE,cex =.1);

circos(R=310-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "tan",] , type="s", col="tan",B=TRUE,cex =.8);

circos(R=280-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "sienna3",] , type="s", col="sienna3",B=TRUE,cex =.8);

#circos(R=250-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "thistle1",] , type="s", col="black",B=TRUE,cex =.8);

circos(R=220-adjust, cir=db, W=50, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "sienna3" | exp.leaf.GeneModules$moduleColors == "tan" ,],col.v=8 , type="b", col="sienna3",B=TRUE,cex =.1);

######################

#Group C
par(mar=c(2,2,2,2))

plot(c(1,800),c(1,800),type="n", axes=FALSE,xlab="",ylab="",main="")

#chromossome names
adjust = 20

circos(R=425-adjust, cir=db, W=0, mapping= mapX, type="label2",side="out", col="grey20" , cex=0.5);

circos(R=400-adjust,cir=db,type="chr",col= colors,print.chr.lab=FALSE,W=40,scale=TRUE)

#All putative genes
circos(R=340-adjust, cir=db, W=50, mapping=density,col.v=4 , type="h", col="black",B=TRUE,cex =.1);

circos(R=310-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "turquoise",] , type="s", col="turquoise",B=TRUE,cex =.4);

circos(R=280-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "darkred",] , type="s", col="darkred",B=TRUE,cex =.4);

circos(R=250-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "coral2",] , type="s", col="coral2",B=TRUE,cex =.4);

circos(R=190-adjust, cir=db, W=50, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "turquoise" | exp.leaf.GeneModules$moduleColors == "darkred" | exp.leaf.GeneModules$moduleColors == "coral2",],col.v=8 , type="b", col="turquoise",B=TRUE,cex =.1);

######################

#Group D
par(mar=c(2,2,2,2))

plot(c(1,800),c(1,800),type="n", axes=FALSE,xlab="",ylab="",main="")

#chromossome names
adjust = 20

circos(R=425-adjust, cir=db, W=0, mapping= mapX, type="label2",side="out", col="grey20" , cex=0.5);

circos(R=400-adjust,cir=db,type="chr",col= colors,print.chr.lab=FALSE,W=40,scale=TRUE)

#All putative genes
circos(R=340-adjust, cir=db, W=50, mapping=density,col.v=4 , type="h", col="black",B=TRUE,cex =.1);

circos(R=310-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "greenyellow",] , type="s", col="greenyellow",B=TRUE,cex =.8);

circos(R=280-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "orange",] , type="s", col="orange",B=TRUE,cex =.8);

circos(R=250-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "honeydew1",] , type="s", col="black",B=TRUE,cex =.8);

circos(R=190-adjust, cir=db, W=50, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "greenyellow" | exp.leaf.GeneModules$moduleColors == "orange" | exp.leaf.GeneModules$moduleColors == "honeydew1",],col.v=8 , type="b", col="orange",B=TRUE,cex =.1);

###########################
#Group E
par(mar=c(2,2,2,2))

plot(c(1,800),c(1,800),type="n", axes=FALSE,xlab="",ylab="",main="")

#chromossome names
adjust = 20

circos(R=425-adjust, cir=db, W=0, mapping= mapX, type="label2",side="out", col="grey20" , cex=0.5);

circos(R=400-adjust,cir=db,type="chr",col= colors,print.chr.lab=FALSE,W=40,scale=TRUE)

#All putative genes
circos(R=340-adjust, cir=db, W=50, mapping=density,col.v=4 , type="h", col="black",B=TRUE,cex =.1);

circos(R=310-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "cyan",] , type="s", col="cyan",B=TRUE,cex =.8);

circos(R=280-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "lightsteelblue1",] , type="s", col="lightsteelblue1",B=TRUE,cex =.8);

circos(R=250-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "violet",] , type="s", col="violet",B=TRUE,cex =.8);

circos(R=220-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "antiquewhite4",] , type="s", col="antiquewhite4",B=TRUE,cex =.8);

circos(R=160-adjust, cir=db, W=50, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "cyan" | exp.leaf.GeneModules$moduleColors == "lightsteelblue1" | exp.leaf.GeneModules$moduleColors == "violet" |  exp.leaf.GeneModules$moduleColors == "antiquewhite4" ,],col.v=8 , type="b", col="antiquewhite4",B=TRUE,cex =.1);


######################

#Group F
par(mar=c(2,2,2,2))

plot(c(1,800),c(1,800),type="n", axes=FALSE,xlab="",ylab="",main="")

#chromossome names
adjust = 20

circos(R=425-adjust, cir=db, W=0, mapping= mapX, type="label2",side="out", col="grey20" , cex=0.5);

circos(R=400-adjust,cir=db,type="chr",col= colors,print.chr.lab=FALSE,W=40,scale=TRUE)

#All putative genes
circos(R=340-adjust, cir=db, W=50, mapping=density,col.v=4 , type="h", col="black",B=TRUE,cex =.1);

circos(R=310-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "pink",] , type="s", col="pink",B=TRUE,cex =.4);

circos(R=280-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "lightpink4",] , type="s", col="lightpink4",B=TRUE,cex =.4);

#circos(R=250-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "thistle1",] , type="s", col="black",B=TRUE,cex =.8);

circos(R=220-adjust, cir=db, W=50, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "pink" | exp.leaf.GeneModules$moduleColors == "lightpink4" ,],col.v=8 , type="b", col="lightpink4",B=TRUE,cex =.1);

######################

#Group G
par(mar=c(2,2,2,2))

plot(c(1,800),c(1,800),type="n", axes=FALSE,xlab="",ylab="",main="")

#chromossome names
adjust = 20

circos(R=425-adjust, cir=db, W=0, mapping= mapX, type="label2",side="out", col="grey20" , cex=0.5);

circos(R=400-adjust,cir=db,type="chr",col= colors,print.chr.lab=FALSE,W=40,scale=TRUE)

#All putative genes
circos(R=340-adjust, cir=db, W=50, mapping=density,col.v=4 , type="h", col="black",B=TRUE,cex =.1);

circos(R=310-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "lavenderblush3",] , type="s", col="lavenderblush3",B=TRUE,cex =.4);

circos(R=280-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "brown4",] , type="s", col="brown4",B=TRUE,cex =.4);

#circos(R=250-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "thistle1",] , type="s", col="black",B=TRUE,cex =.8);

circos(R=220-adjust, cir=db, W=50, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "lavenderblush3" | exp.leaf.GeneModules$moduleColors == "brown4" ,],col.v=8 , type="b", col="brown4",B=TRUE,cex =.1);

######################
#Group H
par(mar=c(2,2,2,2))

plot(c(1,800),c(1,800),type="n", axes=FALSE,xlab="",ylab="",main="")

#chromossome names
adjust = 20

circos(R=425-adjust, cir=db, W=0, mapping= mapX, type="label2",side="out", col="grey20" , cex=0.5);

circos(R=400-adjust,cir=db,type="chr",col= colors,print.chr.lab=FALSE,W=40,scale=TRUE)

#All putative genes
circos(R=340-adjust, cir=db, W=50, mapping=density,col.v=4 , type="h", col="black",B=TRUE,cex =.1);

circos(R=310-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "grey60",] , type="s", col="grey60",B=TRUE,cex =.4);

circos(R=280-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "darkorange",] , type="s", col="darkorange",B=TRUE,cex =.4);

circos(R=250-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "salmon4",] , type="s", col="salmon4",B=TRUE,cex =.4);

circos(R=220-adjust, cir=db, W=25, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "darkseagreen4",] , type="s", col="darkseagreen4",B=TRUE,cex =.4);

circos(R=160-adjust, cir=db, W=50, mapping=exp.leaf.GeneModules[exp.leaf.GeneModules$moduleColors == "grey60" | exp.leaf.GeneModules$moduleColors == "darkorange" | exp.leaf.GeneModules$moduleColors == "salmon4" |  exp.leaf.GeneModules$moduleColors == "darkseagreen4" ,],col.v=8 , type="b", col="grey60",B=TRUE,cex =.1);
