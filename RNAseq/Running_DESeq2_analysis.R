###Set your working directory. Create the folder structure required for the script output

setwd($YOURWORKINGDIRECTORY)

dir.create('data') # Transfer your count tables here

dir.create('plots') # Output plots will be saved here
dir.create('plots/MAplots')
dir.create('plots/PCAplots')
dir.create('plots/VolcanoPlots')

dir.create('results') #Output data tables will be saved here
dir.create('results/DESeq2OutputTables')
dir.create('results/DiffGeneBodyCoordinates')
dir.create('results/DiffRefSeqID')


########Run DESeq2 script for spike-in calibrated RNA-seq data
########This version performs differential gene expression analysis for two conditions (UNT vs TAM) and 3 biological replicates


cell_line = c('BAP1ff', 'PRC1CPM_BAP1ff', 'PRC1CPM') #Cell line names should match the ones in the read count file
rnaseq_app = 'Nuc' #For Nuclear RNA-seq name convention
cutoff=1.5 #Fold-change cut-off for calling DEGs

source('RNAseq_DESeq2.R')
