require(extrafont) #optional. Used to import Arial font into Ubuntu RStudio. Found that if using Windows better to remove it together with family="Arial" from ggplot2 parametre specification
require(ggplot2)
require(reshape2)
require(grid)



font_import()  #optional. Used to import Arial font into Ubuntu RStudio. Found that if using Windows better to remove it together with family="Arial" from ggplot2 parametre specification

####################################################################################################
#####Run Chromatin Express function on H2AK119ub1 Native ChIP-seq data #############################
####################################################################################################

setwd('/YOUR_WD/')

####Run ChromatinExpress.R from the command line. 
####The ChromatinExpress.R script has to be in the PATH. 
####Run from the working directory specified above.
####ChromatinExpress.R takes two arguements: 1) the location of the bigwig file with the ChIP-seq data 2) the name of the sample

#ChromatinExpress.R 'UNT.bw' 'PRC1CKO_UNT_H2AK119ub'

#ChromatinExpress.R 'TAM.bw' 'PRC1CKO_TAM_H2AK119ub'


################################################################################################################################################
#####Plot chromosome density at chr18 #############################################################################################
################################################################################################################################################


####Specify the list of cell lines to run the plotting script on. before the first _ in the RData name

cell_lines=c('PRC1CKO')

####Specify the list of chromosomes for which to plot the density plots.

chr=c(rep('chr18', length(cell_lines)))

####The following function is to plot chromosome density for x cell lines and y chromosomes

####Chromosome denisty is plotted in 250 kb bins. Change this by changing the 250000 factor

chromosome_density=function(x,y){
  
  #Read in the data
  
  load(paste0(x, "_UNT_H2AK119ub.RData"))
  UNT=binsTotal
  
  load(paste0(x, "_TAM_H2AK119ub.RData"))
  TAM=binsTotal
  
  plot=data.frame(A=UNT[which(seqnames(UNT)==y)]$score, B=TAM[which(seqnames(TAM)==y)]$score)
  plot=melt(plot)
  plot$coordinate=rep(seq(from=1, to = length(UNT[which(seqnames(UNT)==y)])*250000, by = 250000),2)
  
  ######Filter from top 5 highest values - usually outliers. Substitute with the mean of the two neighbouring 250 kb bins.
  
  filter_UNT=plot$coordinate[order(plot[which(plot$variable=='A'), c('value')], decreasing=TRUE)[1:5]]
  
  filter_coordinates=plot$coordinate[which(plot$coordinate %in% filter_UNT)]
  
  plot$value[which(plot$coordinate %in% filter_UNT)]=(plot$value[which(plot$coordinate %in% (filter_coordinates-250000))]+plot$value[which(plot$coordinate %in% (filter_coordinates+250000))])/2
  
  plot$variable=factor(plot$variable, levels=c('A', 'B'))
  
  #Plot Data
  
  pdf(paste0('ChromosomeDensity_H2AK119ub_',x, '_', y, '_250kb.pdf'), height=6, width=12)
  p=ggplot(data = plot, aes(y=value, x=coordinate, fill=variable, alpha=variable))+
    theme_bw()+
    theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), plot.title=element_text(size=20, family='Arial'), axis.text.x = element_text(size=20, family='Arial', colour = "black"), axis.text.y = element_text(size=20, family='Arial', colour="black"), axis.title.x=element_text(size=20, family='Arial', colour="black"), axis.title.y=element_text(size=20, family='Arial', colour="black"),legend.title=element_blank(),legend.position=c(0.1, 0.85),legend.key.size = unit(1.5, "lines"), legend.text = element_text(size=15, family="Arial"), panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
    scale_x_continuous(expand = c(0,0), name=y, limits=c(0, length(UNT[which(seqnames(UNT)==y)])*250000), breaks=seq(0, length(UNT[which(seqnames(UNT)==y)])*250000, by = 25000000), labels=paste0(seq(0, length(UNT[which(seqnames(UNT)==y)])*250000, by = 25000000)/1000000, ' Mb'))+
    scale_y_continuous(expand = c(0,0), name='Read Density', limits=c(0, 10))+
    scale_fill_manual(values=c( 'midnightblue', 'red'), labels=c('UNT', 'TAM'))+scale_alpha_manual(values = c(1, 1), guide='none')+
    ggtitle(paste0('HAK119ub in ', x))+theme(plot.title = element_text(size=15, family="Arial", hjust = 0.5))+
    geom_area(position='identity')
  
  plot(p)
  
  dev.off()
  
}

mapply(x=cell_lines, y=chr, FUN=chromosome_density)
