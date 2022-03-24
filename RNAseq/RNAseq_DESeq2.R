##Load required packages

library('DESeq2')
library('ggplot2')
library('ggExtra')
library('plyr')
library('dplyr')
library('reshape2')
library('RColorBrewer')
library('purrr')


##DESeq2 function:

#####Arguments: required to properly read in the files with file names and column names specified in the correct format

print(rnaseq_app) #rnaseq_app: Total or Nuc RNAseq approach

print(cell_line) #cell_line: Cell line in the name of the columns

print(cutoff) #cutoff: Fold change cutoff, for example 1.5

################

deseq_RNAseq=function(rnaseq_app, cell_line, cutoff){
  
  #######Print initialising parametres:
  
  print(paste0('Run DESeq2 on data/dm6.refGene.GeneBody.Uniq.ANNOTATED.', rnaseq_app, 'RNAseq.txt file'))
  
  print(paste0('Run DESeq2 on data/mm10_NonRedundantRefGene_GeneBody.ANNOTATED.', rnaseq_app, 'RNAseq.txt file'))
  
  print(paste0('Run DESeq2 for ', cell_line, ' cell line'))
  
  
  #############################################################
  
  ##Preparing dm6 (spike-in) count table
  
  print(paste0(cell_line, ' ', rnaseq_app, 'RNAseq dm6'))
  
  #Read in spike-in read count table from a set of unique dm6 refGene genes
  
  dm6 <- read.table(file=paste0('data/dm6.refGene.GeneBody.Uniq.ANNOTATED.', rnaseq_app, 'RNAseq.csv'), sep='\t', header=T)
  
  dm6 <- dm6[,c('ID', paste0(cell_line, '_TAM_rep1'), paste0(cell_line, '_TAM_rep2'), paste0(cell_line, '_TAM_rep3'), 
             paste0(cell_line, '_UNT_rep1'), paste0(cell_line, '_UNT_rep2'), paste0(cell_line, '_UNT_rep3'))]
  
  rownames(dm6) <-  make.names(c(as.character(dm6$ID)), unique=TRUE)
  
  dm6 <- dm6[,2:ncol(dm6)]
  
  #Make a DESeq2 sample info table - make sure this matches the actual order of samples
  
  dm6.Condition <- c(rep('TAM',3), rep('UNT',3))
  
  dm6.Rep <- rep(c('rep1', 'rep2', 'rep3'),2)
  
  dm6.SampleInfo <- data.frame(condition = dm6.Condition, rep = dm6.Rep)
  
  rownames(dm6.SampleInfo) <- colnames(dm6)
  
  dm6.SampleInfo$condition  <- factor(dm6.SampleInfo$condition, levels = c('UNT', 'TAM'))
  
  dm6.SampleInfo
  
  ## Obtaining size factors from (prenormalised) dm6 count data for the subsequent normalisation of mm10 counts
  
  dm6.DESEQ2 <- DESeqDataSetFromMatrix(countData = dm6, colData = dm6.SampleInfo, design = ~ rep + condition)
  
  cds_dm6 <- estimateSizeFactors(dm6.DESEQ2)
  
  size_Factors_dm6 <- sizeFactors(cds_dm6)
  
  write.table(size_Factors_dm6, paste0('results/SizeFactors/', cell_line , '.', rnaseq_app, 'RNAseq_Size_Factors_dm6.txt'), quote=FALSE, col.names=FALSE)
  
  ####################################################################################################################
  
  ## Preparing mm10 (experimental) count table
  
  print(paste0(cell_line, ' ', rnaseq_app, 'RNAseq mm10'))
  
  exp <- read.table(paste0('data/mm10_NonRedundantRefGene_GeneBody.ANNOTATED.', rnaseq_app, 'RNAseq.csv'), sep='\t', header=T)
  
  exp <- exp[,c('ID', paste0(cell_line, '_TAM_rep1'), paste0(cell_line, '_TAM_rep2'), paste0(cell_line, '_TAM_rep3'), 
             paste0(cell_line, '_UNT_rep1'), paste0(cell_line, '_UNT_rep2'), paste0(cell_line, '_UNT_rep3'))]
  
  rownames(exp) <- make.names(c(as.character(exp$ID)), unique=TRUE)
  
  ##Preparing mm10 sample info table
  
  exp <- exp[,2:ncol(exp)]
  
  exp.Condition <- c(rep('TAM',3), rep('UNT',3))
  
  exp.Rep <- rep(c('rep1', 'rep2', 'rep3'), times=2)
  
  exp.SampleInfo <- data.frame(condition=exp.Condition, rep = exp.Rep)
  
  rownames(exp.SampleInfo) <- colnames(exp)
  
  exp.SampleInfo$condition  <- factor(exp.SampleInfo$condition, levels = c('UNT', 'TAM'))
  
  exp.SampleInfo
  
  #Performing Differential analysis
  
  exp.DESEQ2 <- DESeqDataSetFromMatrix(countData=exp, colData=exp.SampleInfo, design = ~ rep + condition) #variable of interest (condition) has to be last in the design formula
  
  sizeFactors(exp.DESEQ2) <- size_Factors_dm6
  
  exp.DESEQ2.analysis <- DESeq(exp.DESEQ2)
  
  
  ##################################################################################################
  
  ##Principal component plot of the samples
  
  se <- SummarizedExperiment(log2(counts(exp.DESEQ2, normalized=TRUE) + 1),colData=colData(exp.DESEQ2))
  
  pcaData <- plotPCA(DESeqTransform(se), ntop=1000, intgroup=c('condition', 'rep'), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, 'percentVar'))
  
  pdf(paste0('plots/PCAplots/', cell_line, '.', rnaseq_app, 'RNAseq_dm6norm.PCA.pdf'))
  p=ggplot(pcaData, aes(PC1, PC2, color=condition, shape=rep)) +
    geom_point(size=3) +
    xlab(paste0('PC1: ',percentVar[1],'% variance')) +
    ylab(paste0('PC2: ',percentVar[2],'% variance')) + 
    coord_fixed()
  print(p)
  dev.off()
  
  
  ##############################################################################
  
  ## Extracting comparisons between conditions
  
  #Compare TAM vs UNT
  
  ###Raw LFC
  
  exp.DESEQ2.results <- results(exp.DESEQ2.analysis, 
                                contrast=c('condition','TAM','UNT'),
                                alpha=0.05) 
  
  ###Normal shrinking using contrast LFC -> for differential gene expression analysis of data with global transcriptional changes, apeglm type of shrinking was shown to over-normalize LFC towards LFC = 0
  
  exp.DESEQ2.results.LFC.normal <- lfcShrink(exp.DESEQ2.analysis, 
                                            contrast=c('condition','TAM','UNT'),
                                            type = 'normal')
  
  ###Add shrunk LFC to the results table
  ###log2FoldChange - raw
  ###LFC_normal - normal shrinking
  
  exp.DESEQ2.results$LFC_normal <- exp.DESEQ2.results.LFC.normal$log2FoldChange
  
  #Preparing to save the results table
  
  print('Preparing to save the results table')
  
  #Add ID column
  
  exp.DESEQ2.results$ID <- row.names(exp.DESEQ2.results)
  
  #Annotate results table with counts
  
  normalised.counts <- as.data.frame(counts(exp.DESEQ2, normalized = TRUE))
  
  normalised.counts$ID <- row.names(normalised.counts)
  
  exp.DESEQ2.results.counts <- merge(as.data.frame(exp.DESEQ2.results), 
                                  normalised.counts, 
                                  by = 'ID')
  
  Output.Table <- exp.DESEQ2.results.counts %>% 
    dplyr::select('ID', everything())
  
  #Adding Gene name and coordinates information and ordering the output table by padj
  
  annotate <- read.table(paste0('data/mm10_NonRedundantRefGene_GeneBody.ANNOTATED.', rnaseq_app, 'RNAseq.csv'), sep='\t', header=T)
  
  Output.Table <- merge(annotate[, c('Chr', 'Start', 'Stop', 
                                     'Strand', 'ID', 'Size', 'Gene_name')], 
                        Output.Table, by='ID' )
  
  Output.Table <- Output.Table[order(Output.Table$padj),]
  
  ##Calculate RPKMs, add +1 pseudocount
  
  reps <- colnames(select(Output.Table,contains('Rep')))
  
  Output.Table[,paste0(reps, '_RPKM')] <- lapply(reps, function(x) (Output.Table[,x]+1) / (Output.Table[, 'Size']/1000))
  
  
  ##Calculate mean RPKMs for conditions
  
   Output.Table[, c(paste0(cell_line, '_UNT_RPKM'))] <- rowMeans(Output.Table[, c(paste0(cell_line, '_UNT_rep1_RPKM'), paste0(cell_line, '_UNT_rep2_RPKM'), paste0(cell_line, '_UNT_rep3_RPKM'))])
  
  Output.Table[, c(paste0(cell_line, '_TAM_RPKM'))] <- rowMeans(Output.Table[, c(paste0(cell_line, '_TAM_rep1_RPKM'), paste0(cell_line, '_TAM_rep2_RPKM'), paste0(cell_line, '_TAM_rep3_RPKM'))])
  
  
  ####Save the final output table
  
  write.table(as.data.frame(Output.Table),file=paste0('results/DESeq2OutputTables/', cell_line, '.', rnaseq_app, 'RNAseq_spikenormalised_DESeq2.txt'), quote=FALSE, row.names=FALSE, sep = '\t')
  
  #######################################
  
  print('Saving results')
  
  ###Extract IDs, Gene names and GeneBody coordinates for differentially expressed genes
  
  ###Use normal-shrunk LFC and p-adj for calling significant gene expression changes
  
  UP.normal <- Output.Table[which(Output.Table$padj<0.05 & Output.Table$LFC_normal > log2(cutoff)),]
  
  DOWN.normal <- Output.Table[which(Output.Table$padj<0.05 & Output.Table$LFC_normal < -log2(cutoff)),]
  
  NonSign.normal <- Output.Table[which(!(Output.Table$ID %in% UP.normal$ID | Output.Table$ID %in% DOWN.normal$ID)),]
  
  
  ####Output gene body coordinates for DEGs
  
  write.table(as.data.frame(UP.normal[, c('Chr', 'Start', 'Stop', 'ID', 'ID', 'Strand')]),file=paste0('results/DiffGeneBodyCoordinates/', cell_line, '_', rnaseq_app, 'RNAseq.DESeq2.UP.LFC_normal.', cutoff, '.bed'), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  write.table(as.data.frame(DOWN.normal[, c('Chr', 'Start', 'Stop', 'ID', 'ID', 'Strand')]),file=paste0('results/DiffGeneBodyCoordinates/', cell_line, '_', rnaseq_app, 'RNAseq.DESeq2.DOWN.LFC_normal.', cutoff, '.bed'), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  write.table(as.data.frame(NonSign.normal[, c('Chr', 'Start', 'Stop', 'ID', 'ID', 'Strand')]),file=paste0('results/DiffGeneBodyCoordinates/', cell_line, '_', rnaseq_app, 'RNAseq.DESeq2.NonSign.LFC_normal.', cutoff, '.bed'), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  #####Output Gene IDs for DEGs
  
  write.table(as.data.frame(UP.normal[, c('ID')]),file=paste0('results/DiffRefSeqID/', cell_line, '_', rnaseq_app, 'RNAseq.UP.LFC_normal.', cutoff, '.txt'), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  write.table(as.data.frame(DOWN.normal[, c('ID')]),file=paste0('results/DiffRefSeqID/', cell_line, '_', rnaseq_app, 'RNAseq.DOWN.LFC_normal.', cutoff, '.txt'), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  write.table(as.data.frame(NonSign.normal[, c('ID')]),file=paste0('results/DiffRefSeqID/', cell_line, '_', rnaseq_app, 'RNAseq.NonSign.LFC_normal.', cutoff, '.txt'), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  
  #####Output Gene Named for DEGs
  
  write.table(as.data.frame(UP.normal[, c('Gene_name')]),file=paste0('results/DiffRefSeqID/', cell_line, '_', rnaseq_app, 'RNAseq.Gene_name.UP.LFC_normal.', cutoff, '.txt'), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  write.table(as.data.frame(DOWN.normal[, c('Gene_name')]),file=paste0('results/DiffRefSeqID/', cell_line, '_', rnaseq_app, 'RNAseq.Gene_name.DOWN.LFC_normal.', cutoff, '.txt'), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  write.table(as.data.frame(NonSign.normal[, c('Gene_name')]),file=paste0('results/DiffRefSeqID/', cell_line, '_', rnaseq_app, 'RNAseq.Gene_name.NonSign.LFC_normal.', cutoff, '.txt'), quote=FALSE, row.names=FALSE, col.names=FALSE, sep = '\t')
  
  
  ###Plot MA plots and Vulcano Plots
  
  print('Making plots')
  
  ###MA plots
  
  ####normal LFC
  
  Output.Table[, paste0(cell_line, '_UNT_RPKM_log2')] <- log2(Output.Table[, paste0(cell_line, '_UNT_RPKM')])
  
  pdf(paste0('plots/MAplots/',cell_line, '_', rnaseq_app, 'RNAseq_density_scatterplot.normal.LFC', cutoff, '.pdf'))
  
  p=ggplot(data=Output.Table,aes_string(paste0(cell_line, '_UNT_RPKM_log2'),'LFC_normal'))+
    theme_bw()+
    theme(aspect.ratio = 1, axis.ticks = element_line(colour = 'black', size = 2), 
          axis.text.x = element_text(colour = 'black',size=20), 
          axis.text.y = element_text(colour = 'black', size=20), 
          axis.title.x=element_text(colour = 'black', size=20), 
          axis.title.y=element_text(colour = 'black', size=20), 
          plot.title = element_text(colour = 'black', size=15,hjust = 0.5), 
          legend.title=element_blank(),legend.position=c(0.85,0.2),
          legend.key = element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(), 
          panel.border = element_rect(colour = 'black', size = 2))+
    geom_point(colour='gray57',alpha=0.35,size=1.75, pch=16)+
    geom_point(data=Output.Table[which(Output.Table$padj<0.05 & abs(Output.Table$LFC_normal) > log2(cutoff)),],colour='red',alpha=0.35,size=1.75, pch=16)+
    scale_y_continuous(limits=c(-5,6),breaks = seq(-5,6, by = 1))+
    scale_x_continuous(limits=c(-12,16), breaks = seq(-12,16, by = 2))+
    xlab('Log2 RPKM UNT')+
    ylab('Log2 Fold Change')+
    ggtitle(paste0(cell_line, ' ', rnaseq_app, '\n', 'UP=', nrow(UP.normal), ' DOWN=', nrow(DOWN.normal)))+
    geom_hline(yintercept=0,linetype='dashed',lwd=2)+
    theme(legend.position='none')
  print(ggMarginal(p, type='density', margins='y', col='black', size=6, lwd=1.5, fill='gray57'))
  
  dev.off()
  
  
  ###Volcano plots
  
  Output.Table$log10_padj = -log10(Output.Table$padj)
  
  pdf(paste0('plots/VolcanoPlots/', cell_line, '_', rnaseq_app, 'RNAseq_volcano_scatterplot.normal.LFC', cutoff, '.pdf'))
  
  p=ggplot(data=Output.Table,aes(LFC_normal,log10_padj))+
    theme_bw()+
    theme(aspect.ratio = 1, axis.ticks = element_line(colour = 'black', size = 2), 
          axis.text.x = element_text(colour = 'black',size=20), 
          axis.text.y = element_text(colour = 'black', size=20), 
          axis.title.x=element_text(colour = 'black', size=20), 
          axis.title.y=element_text(colour = 'black', size=20), 
          plot.title = element_text(colour = 'black', size=15,hjust = 0.5), 
          legend.title=element_blank(),legend.position=c(0.85,0.2),
          legend.key = element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(), 
          panel.border = element_rect(colour = 'black', size = 2))+
    geom_point(colour='gray57',alpha=0.4,size=1.5,pch=16)+
    geom_point(data=Output.Table[which(Output.Table$padj<0.05 & abs(Output.Table$LFC_normal)>log2(cutoff)),],colour='red',alpha=0.5,size=1.5,pch=16)+
    scale_x_continuous(limits=c(-5,6), breaks = seq(-5,6, by = 2))+
    xlab('Log2 Fold Change')+
    ylab('-log10 p-adj')+
    ggtitle(paste0(cell_line, ' ', rnaseq_app, '\n', 'UP=', nrow(UP.normal), ' DOWN=', nrow(DOWN.normal)))+
    geom_vline(xintercept=0,linetype='dashed',lwd=2)+
    theme(legend.position='none')
  print(p)
  
  dev.off()
  
}


mapply(FUN = deseq_RNAseq, rnaseq_app = rnaseq_app, cell_line = cell_line, cutoff = cutoff)




