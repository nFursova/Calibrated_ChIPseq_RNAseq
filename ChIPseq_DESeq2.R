#################################################################################################################################################################
###################Running DESeq2 to compare ChIP enrichment in one cell line in two condtions (UNT and TAM) for 3 biological reps###############################
#################################################################################################################################################################

####Load required packages


library("DESeq2")
library('plyr')
library('dplyr')
library('ggplot2')
library('reshape2')
library('RColorBrewer')
library('ggExtra')
library('purrr')
library('VennDiagram')
library('grid')
library('apeglm')


####To use it with the counts produced in multiBamSummary, first get rid of "#" on the first line

####Arguments: required to correctly read in files with file names and column names specified in the correct format

print(CellLine) #Name of Cell Line
print(chip) # Type of ChIP-seq - XChIPseq or NChIPSeq
print(line_factor) # Name of the cell line and line_factor/modification to analyse, i.e. 'PRC1CPM_SUZ12'
print(intervals) #Name of intervals to use for the analysis, i.e. 'mESC_PRC1cKO_PRC1CPM_UNT_RING1B_peaks'
print(spikein_norm) #Name of the spikein intervals to use for spikein_norm in DESeq2, i.e. 'hg19.CpG_Islands.annotation'
print(cutoff) #Fold change cutoff to use for significant changes, i.e. 2




deseq=function(CellLine, chip, line_factor, intervals, spikein_norm, cutoff){
  
  #############################################################
  
  #Preparing spike-in counts table
  
  print(paste0(line_factor, '_', spikein_norm))
  
  hg19=read.table(file=paste0('data/', CellLine ,'.', chip, '.', spikein_norm, '.Normalised.Read.Counts.txt'), sep='\t', header=T, stringsAsFactors = F)
  
  hg19$ID = paste0('Peak_', 1:nrow(hg19))
  
  hg19=hg19[,c("ID", paste0(line_factor, "_TAM_rep1"), paste0(line_factor, "_TAM_rep2"), paste0(line_factor, "_TAM_rep3"), paste0(line_factor, "_UNT_rep1"), paste0(line_factor, "_UNT_rep2"), paste0(line_factor, "_UNT_rep3"))]
  
  rownames(hg19) = make.names(c(as.character(hg19$ID)), unique=TRUE)
  
  hg19=hg19[,2:ncol(hg19)]
  
  hg19.Condition=c('TAM', 'TAM', 'TAM', 'UNT', 'UNT', 'UNT')
  
  hg19.Rep=c('rep1', 'rep2', 'rep3', 'rep1','rep2', 'rep3')
  
  hg19.SampleInfo=data.frame(condition= hg19.Condition, rep= hg19.Rep)
  
  rownames(hg19.SampleInfo)<-colnames(hg19)
  
  hg19.SampleInfo
  
  ## Obtaining size factors from (prenormalised) hg19 count data in hg19 CGIs for spikein_norm of mm10 counts
  
  hg19.DESEQ2 <-DESeqDataSetFromMatrix(countData=hg19, colData=hg19.SampleInfo, design = ~ rep + condition)
  
  cds_hg19=estimateSizeFactors(hg19.DESEQ2)
  
  size_Factors_hg19=sizeFactors(cds_hg19)
  
  write.table(size_Factors_hg19, paste0('results/SizeFactors/', line_factor, '_', intervals, '_', spikein_norm  ,'_Size_Factors.txt'), quote=FALSE, col.names=FALSE)
  
####################################################################################################################
  
## preparing mm10 counts table and sample info table
  
  print(paste0(line_factor, '_', intervals))
  
  exp=read.table(paste0('data/', CellLine ,'.', chip, '.', intervals, '.Raw.Read.Counts.txt'), sep='\t', header=T)
  
  exp$ID = paste0('Peak_', 1:nrow(exp))
  
  exp=exp[c("ID", paste0(line_factor, "_TAM_rep1"), paste0(line_factor, "_TAM_rep2"), paste0(line_factor, "_TAM_rep3"), paste0(line_factor, "_UNT_rep1"), paste0(line_factor, "_UNT_rep2"), paste0(line_factor, "_UNT_rep3"))]
  
  rownames(exp) = make.names(c(as.character(exp$ID)), unique=TRUE)
  
  exp=exp[,2:ncol(exp)]
  
  exp.Condition=c('TAM', 'TAM', 'TAM', 'UNT', 'UNT', 'UNT')
  
  exp.Rep=c('rep1', 'rep2', 'rep3', 'rep1','rep2', 'rep3')
  
  exp.SampleInfo=data.frame(condition = exp.Condition, rep = exp.Rep)
  
  rownames(exp.SampleInfo) <- colnames(exp)
  
  exp.SampleInfo$condition = factor(exp.SampleInfo$condition, levels = c('UNT', 'TAM'))
  
  exp.SampleInfo
  
  #Performing Differential analysis
  
  exp.DESEQ2 <- DESeqDataSetFromMatrix(countData=exp, colData=exp.SampleInfo, design = ~ rep + condition)
  
  sizeFactors(exp.DESEQ2) <- size_Factors_hg19
  
  exp.DESEQ2.analysis <- DESeq(exp.DESEQ2)
  
  
##################################################################################################
  
## Principal component plot of the samples
  
  se = SummarizedExperiment(log2(counts(exp.DESEQ2, normalized=TRUE) + 1),colData=colData(exp.DESEQ2))
  
  pcaData <- plotPCA(DESeqTransform(se), ntop=1000, intgroup=c("condition", "rep"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  pdf(paste0('plots/PCAplots/', line_factor, '_', intervals, '_', spikein_norm, '_norm.DESeq2_spikein.PCA.pdf'))
  p=ggplot(pcaData, aes(PC1, PC2, color=condition, shape=rep)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  print(p)
  dev.off()
  
###################################################################################
  
  # Extracting comparisons between conditions
  
  
  exp.DESEQ2.results <- results(exp.DESEQ2.analysis, name = 'condition_TAM_vs_UNT',alpha=0.05) #by default MLE
  
  exp.DESEQ2.results.LFC=lfcShrink(exp.DESEQ2.analysis, coef = 'condition_TAM_vs_UNT', res=exp.DESEQ2.results, type = "apeglm")
  
  exp.DESEQ2.results.LFC.old=lfcShrink(exp.DESEQ2.analysis, coef = 'condition_TAM_vs_UNT', res=exp.DESEQ2.results, type = "normal")
  
  exp.DESEQ2.results$LFC=exp.DESEQ2.results.LFC$log2FoldChange
  
  exp.DESEQ2.results$LFC_old=exp.DESEQ2.results.LFC.old$log2FoldChange
  
  
   #Annotate table with normalised counts
  
  normalised.counts = as.data.frame(counts(exp.DESEQ2, normalized = TRUE))
  
  normalised.counts$ID = row.names(normalised.counts)
  
  exp.DESEQ2.results$ID = row.names(exp.DESEQ2.results)
  
  exp.DESEQ2.results.counts=merge(as.data.frame(exp.DESEQ2.results), normalised.counts, by = 'ID')
  
  #Output summary of stats
  
  sink(file=paste0('results/summary/',line_factor, '_', intervals, '_', spikein_norm, ".summary_BE.txt"))
  summary(exp.DESEQ2.results,alpha=0.05)
  sink()
  
  #####Make a dataframe for export
  
  Output.Table <- exp.DESEQ2.results.counts %>% select('ID', everything())
  
  ########Annotate the table with the coordinates of intervals
  
  exp=read.table(paste0('data/', CellLine ,'.', chip, '.', intervals, '.Raw.Read.Counts.txt'), sep='\t', header=T)
  
  exp$ID = paste0('Peak_', 1:nrow(exp))
  
  Output.Table <- merge(exp[, c('chr', 'start', 'end', 'ID')], Output.Table, by='ID' )
  
  ###Sort the table based on p-adj
  
  Output.Table <- Output.Table[order(Output.Table$padj),]
  
  ###Calculate the size of the intervals
  
  Output.Table$Size = Output.Table$end - Output.Table$start + 1
  
  
  ########Export Peak coordinates of significantly differentially enriched peaks
  
  ###Make separate tables for each group of peaks
  
  ###APEGLM LFC
  
  UP <- Output.Table[which(Output.Table$padj<0.05 & Output.Table$LFC > log2(cutoff)),]
  
  DOWN <- Output.Table[which(Output.Table$padj<0.05 & Output.Table$LFC < -log2(cutoff)),]
  
  NonSign <- Output.Table[which(!(Output.Table$ID %in% UP$ID | Output.Table$ID %in% DOWN$ID)),]
  
  ###Save tables
  
  write.table(as.data.frame(UP[, c('chr', 'start', 'end')]),file=paste0('results/DiffPeakSets/', line_factor, '_', intervals, '_', spikein_norm, ".DESeq2.UP.FC", cutoff, ".bed"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  write.table(as.data.frame(DOWN[, c('chr', 'start', 'end')]),file=paste0('results/DiffPeakSets/', line_factor, '_', intervals, '_', spikein_norm, ".DESeq2.DOWN.FC", cutoff, ".bed"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  write.table(as.data.frame(NonSign[, c('chr', 'start', 'end')]),file=paste0('results/DiffPeakSets/', line_factor, '_', intervals, '_', spikein_norm, ".DESeq2.NonSign.FC", cutoff, ".bed"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
###OLD LFC
  
  UP.old <- Output.Table[which(Output.Table$padj<0.05 & Output.Table$LFC_old > log2(cutoff)),]
  
  DOWN.old <- Output.Table[which(Output.Table$padj<0.05 & Output.Table$LFC_old < -log2(cutoff)),]
  
  NonSign.old <- Output.Table[which(!(Output.Table$ID %in% UP.old$ID | Output.Table$ID %in% DOWN.old$ID)),]
  
  ###Save tables
  
  write.table(as.data.frame(UP.old[, c('chr', 'start', 'end')]),file=paste0('results/DiffPeakSets/', line_factor, '_', intervals, '_', spikein_norm, ".DESeq2.UP.old.FC", cutoff, ".bed"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  write.table(as.data.frame(DOWN.old[, c('chr', 'start', 'end')]),file=paste0('results/DiffPeakSets/', line_factor, '_', intervals, '_', spikein_norm, ".DESeq2.DOWN.old.FC", cutoff, ".bed"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  write.table(as.data.frame(NonSign.old[, c('chr', 'start', 'end')]),file=paste0('results/DiffPeakSets/', line_factor, '_', intervals, '_', spikein_norm, ".DESeq2.NonSign.old.FC", cutoff, ".bed"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  
  
  ###Calculate RPKMs
  
  ###UNT
  
  Output.Table[, c(paste0(line_factor, "_UNT_rep1_RPKM"))] <- Output.Table[, c(paste0(line_factor, "_UNT_rep1"))]/((Output.Table[, c("Size")])/1000)
  Output.Table[, c(paste0(line_factor, "_UNT_rep2_RPKM"))] <- Output.Table[, c(paste0(line_factor, "_UNT_rep2"))]/((Output.Table[, c("Size")])/1000)
  Output.Table[, c(paste0(line_factor, "_UNT_rep3_RPKM"))] <- Output.Table[, c(paste0(line_factor, "_UNT_rep3"))]/((Output.Table[, c("Size")])/1000)
  
  Output.Table$UNT_RPKM <- (Output.Table[, c(paste0(line_factor, "_UNT_rep1_RPKM"))]+Output.Table[, c(paste0(line_factor, "_UNT_rep2_RPKM"))]+Output.Table[, c(paste0(line_factor, "_UNT_rep3_RPKM"))])/3
  
  Output.Table$UNT_log2RPKM <- log2(Output.Table$UNT_RPKM)
  
  
  ###TAM
  
  Output.Table[, c(paste0(line_factor, "_TAM_rep1_RPKM"))] <- Output.Table[, c(paste0(line_factor, "_TAM_rep1"))]/((Output.Table[, c("Size")])/1000)
  Output.Table[, c(paste0(line_factor, "_TAM_rep2_RPKM"))] <- Output.Table[, c(paste0(line_factor, "_TAM_rep2"))]/((Output.Table[, c("Size")])/1000)
  Output.Table[, c(paste0(line_factor, "_TAM_rep3_RPKM"))] <- Output.Table[, c(paste0(line_factor, "_TAM_rep3"))]/((Output.Table[, c("Size")])/1000)
  
  Output.Table$TAM_RPKM <- (Output.Table[, c(paste0(line_factor, "_TAM_rep1_RPKM"))]+Output.Table[, c(paste0(line_factor, "_TAM_rep2_RPKM"))]+Output.Table[, c(paste0(line_factor, "_TAM_rep3_RPKM"))])/3
  
  Output.Table$TAM_log2RPKM <- log2(Output.Table$TAM_RPKM)
  
  
  ###Change the column names prior to exporting the final Output tables
  
  Output.Table.renamed = Output.Table
  
  
  colnames(Output.Table.renamed)[which(colnames(Output.Table.renamed) == 'LFC')] = paste0('LFC_', line_factor)
  
  colnames(Output.Table.renamed)[which(colnames(Output.Table.renamed) == 'LFC_old')] = paste0('LFC_old_', line_factor)
  
  colnames(Output.Table.renamed)[which(colnames(Output.Table.renamed) == 'padj')] = paste0('padj_', line_factor)
  
  colnames(Output.Table.renamed)[which(colnames(Output.Table.renamed) == 'pvalue')] = paste0('pvalue_', line_factor)
  
  colnames(Output.Table.renamed)[which(colnames(Output.Table.renamed) == 'UNT_RPKM')] = paste0('UNT_RPKM_', line_factor)
  
  colnames(Output.Table.renamed)[which(colnames(Output.Table.renamed) == 'TAM_RPKM')] = paste0('TAM_RPKM_', line_factor)
  
  colnames(Output.Table.renamed)[which(colnames(Output.Table.renamed) == 'UNT_log2RPKM')] = paste0('UNT_log2RPKM_', line_factor)
  
  colnames(Output.Table.renamed)[which(colnames(Output.Table.renamed) == 'TAM_log2RPKM')] = paste0('TAM_log2RPKM_', line_factor)
  
  
  write.table(as.data.frame(Output.Table.renamed),file=paste0('results/DESeq2OutputTables/', line_factor, '_', intervals, '_', spikein_norm, ".spikenormalised_DESeq2_BE.csv"), quote=FALSE, row.names=FALSE, sep = '\t')
  
  
  ####plot RK plots
  
  ####################plot MA plot with the new apeglm LFC
  
  pdf(paste0('plots/MAplots/',line_factor, '_', intervals, '_', spikein_norm, '_density_scatterplot.apeglm.LFC', cutoff, '.pdf'))
  p=ggplot(data=Output.Table,aes(UNT_log2RPKM,LFC))+theme_bw()+
    theme(aspect.ratio = 1, axis.ticks = element_line(colour = "black", size = 2), axis.text.x = element_text(colour = 'black',size=20), axis.text.y = element_text(colour = 'black', size=20), axis.title.x=element_text(colour = 'black', size=20), axis.title.y=element_text(colour = 'black', size=20), plot.title = element_text(colour = 'black', size=15,hjust = 0.5), legend.title=element_blank(),legend.position=c(0.85,0.2),legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border = element_rect(colour = 'black', size = 2))+
    geom_point(colour="gray57",alpha=0.35,size=1.75, pch=16)+
    geom_point(data=Output.Table[which(Output.Table$padj<0.05 & abs(Output.Table$LFC) > log2(cutoff)),],colour="red",alpha=0.35,size=1.75, pch=16)+
    scale_y_continuous(limits=c(-4,4),breaks = seq(-4,4, by = 1))+
    scale_x_continuous(limits=c(0,12), breaks = seq(0,12, by = 2))+
    xlab("Log2 Mean of UNT RPKM")+ylab("Log2 Fold Change")+ggtitle(paste0(line_factor, ' ChIP-seq', '\n', 'UP=', nrow(UP), ' DOWN=', nrow(DOWN)))+
    geom_hline(yintercept=0,linetype="dashed",lwd=2)+
    theme(legend.position="none")
  print(ggMarginal(p, type="density", margins='y', col='black', size=6, lwd=1.5, fill='gray57'))
  dev.off()
  

  pdf(paste0('plots/MAplots/',line_factor, '_', intervals, '_', spikein_norm, '_density_scatterplot.old.LFC', cutoff, '.pdf'))
  p=ggplot(data=Output.Table,aes(UNT_log2RPKM,LFC_old))+theme_bw()+
    theme(aspect.ratio = 1, axis.ticks = element_line(colour = "black", size = 2), axis.text.x = element_text(colour = 'black',size=20), axis.text.y = element_text(colour = 'black', size=20), axis.title.x=element_text(colour = 'black', size=20), axis.title.y=element_text(colour = 'black', size=20), plot.title = element_text(colour = 'black', size=15,hjust = 0.5), legend.title=element_blank(),legend.position=c(0.85,0.2),legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border = element_rect(colour = 'black', size = 2))+
    geom_point(colour="gray57",alpha=0.35,size=1.75, pch=16)+
    geom_point(data=Output.Table[which(Output.Table$padj<0.05 & abs(Output.Table$LFC_old) > log2(cutoff)),],colour="red",alpha=0.35,size=1.75, pch=16)+
    scale_y_continuous(limits=c(-4,4),breaks = seq(-4,4, by = 1))+
    scale_x_continuous(limits=c(0,12), breaks = seq(0,12, by = 2))+
    xlab("Log2 Mean of UNT RPKM")+ylab("Log2 Fold Change")+ggtitle(paste0(line_factor, ' ChIP-seq', '\n', 'UP=', length(which(Output.Table$padj < 0.05 & Output.Table$LFC_old > log2(cutoff))), ' DOWN=', length(which(Output.Table$padj < 0.05 & Output.Table$LFC_old <  -log2(cutoff)))))+
    geom_hline(yintercept=0,linetype="dashed",lwd=2)+
    theme(legend.position="none")
  print(ggMarginal(p, type="density", margins='y', col='black', size=6, lwd=1.5, fill='gray57'))
  dev.off()
  
}

mapply(FUN=failwith(NULL, deseq), line_factor=line_factor, intervals=intervals, spikein_norm=spikein_norm, cutoff = cutoff, chip = chip, CellLine = CellLine)

