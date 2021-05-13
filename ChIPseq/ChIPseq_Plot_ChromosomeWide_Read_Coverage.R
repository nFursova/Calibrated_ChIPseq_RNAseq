library('plyr')
library('dplyr')
require('ggplot2')
require('reshape2')
library('ggExtra')


####################################################################################################
###############Run ChIPseq_Chromosome_Wide_Read_Coverage.R function on ChIP-seq data ###############
####################################################################################################

setwd('/YOUR_WD/')

####Run ChIPseq_Chromosome_Wide_Read_Coverage.R from the command line in the working directory specified above
####ChIPseq_Chromosome_Wide_Read_Coverage.R script has to be in the PATH
####ChIPseq_Chromosome_Wide_Read_Coverage.R takes two arguments: 1) location of the bigwig file 2) the name of the sample
####The output is saved in the working directory 

#ChIPseq_Chromosome_Wide_Read_Coverage.R 'CellLine_UNT_ChIPseq.bw' 'UNT'

#ChIPseq_Chromosome_Wide_Read_Coverage.R 'CellLine_TAM_ChIPseq.bw' 'TAM'


################################################################################################################################################
#####Plot chromosome density ###################################################################################################################
################################################################################################################################################

####Specify the list of cell lines to plot chromosome density for. before the first '_' in the RData name

cell_lines=c('BAP1ff')

####Specify the list of chromosomes to plot chromosome density for.

chr=c('chr18')

####Specify the list of modifications/factors profiled by ChIP-seq to plot chromosome density for.

factor = c('H2AK119ub1')

####The following function is to plot chromosome density for x cell lines, y chromosomes and z modifications/factors

####Chromosome denisty is plotted in 250 kb bins as was calculated by .

chromosome_density=function(x,y,z){
  
  print(x)
  print(y)
  print(z)
  
  #Read in the data
  
  load(paste0(x, "_UNT_", z, ".RData"))
  UNT=binsTotal
  
  load(paste0(x, "_TAM_", z, ".RData"))
  TAM=binsTotal
  
  plot=data.frame(UNT=UNT[which(seqnames(UNT)==y)]$score, TAM=TAM[which(seqnames(TAM)==y)]$score)
  plot=melt(plot)
  plot$coordinate=rep(seq(from=1, to = length(UNT[which(seqnames(UNT)==y)])*250000, by = 250000),2)
  
  #####Order for plotting: first plot TAM as it has higher values in this case and looks better when it is at the back
  #####This can be changed by changing the order of factor levels
  
  plot$variable=factor(plot$variable, levels=c('TAM', 'UNT'))
  
  #####Calculate the maximum value for the y-scale
  
  ymax = round(max(plot$value), digits = -1)*1.5
  
  #Plot Data
  
  pdf(paste0('ChromosomeDensity_', z, '_',x, '_', y, '_250kb.pdf'), height=6, width=12)
  
  p=ggplot(data = plot, aes(y=value, x=coordinate, fill=variable))+
    theme_bw()+
    theme(axis.line = element_line(colour = 'black', size = 1.5), 
          axis.ticks = element_line(colour = "black", size = 1.5), 
          plot.title = element_text(size=15, hjust = 0.5), 
          axis.text.x = element_text(size=20, colour = "black"), 
          axis.text.y = element_text(size=20, colour="black"), 
          axis.title.x=element_text(size=20, colour="black"), 
          axis.title.y=element_text(size=20, colour="black"),
          legend.title=element_blank(),legend.position=c(0.1, 0.85),
          legend.key.size = unit(1.5, "lines"), legend.text = element_text(size=15), 
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
    scale_x_continuous(expand = c(0,0), name=y, 
                       limits=c(0, length(UNT[which(seqnames(UNT)==y)])*250000+1), 
                       breaks=seq(0, length(UNT[which(seqnames(UNT)==y)])*250000, by = 25000000), 
                       labels=paste0(seq(0, length(UNT[which(seqnames(UNT)==y)])*250000, by = 25000000)/1000000, ' Mb'))+
    scale_y_continuous(expand = c(0,0), name='Read Density', limits=c(0, ymax), breaks = seq(0, ymax, by = 5))+
    scale_fill_manual(values=c('#ed2224', '#29286b'), labels=c('TAM', 'UNT'))+
    scale_alpha_manual(values = c(1, 1), guide='none')+
    ggtitle(paste0(z, ' in ', x))+
    geom_area(position='identity')
  plot(p)
  
  dev.off()
  
}

mapply(FUN=chromosome_density, x=cell_lines, y=chr, z=factor)
