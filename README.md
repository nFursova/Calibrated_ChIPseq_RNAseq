# Calibrated-ChIPseq-RNAseq-tools

Scripts developed in Klose Lab to perform analysis of spike-in calibrated ChIP-seq and RNA-seq data. Contributed by Dr Hamish King (drhamishking@gmail.com), Dr Nadezda Fursova (nfursova.msu@gmail.com), and Dr Anne Turberfield. Used for the data analysis in the manuscript "BAP1 constrains pervasive H2AK119ub1 to control the transcriptional potential of the genome".

Please contact nadya.fursova@nih.gov or nfursova.msu@gmail.com with questions or problems.

Requirements for each script are listed individually, but most scripts rely on the use of basic read aligment and processing tools:

- FastQC	https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Samtools	http://samtools.sourceforge.net/ 
- BEDtools	https://bedtools.readthedocs.io/en/latest/
- STAR	https://github.com/alexdobin/STAR
- Sambamba	https://lomereiter.github.io/sambamba/
- WiggleTools	https://github.com/Ensembl/WiggleTools
- UCSC tools http://hgdownload.cse.ucsc.edu/admin/exe/
- DeepTools	https://deeptools.readthedocs.io/en/develop/index.html
- HOMER	http://homer.ucsd.edu/homer/index.html
- R	https://www.r-project.org/
- MultiQC	https://multiqc.info/

Scripts performed on a list of files require the input in tab-limited format as indicated below:

file.bam	SampleName

# Summary of the read alignment and processing pipelines for calibrated ChIP-seq and RNA-seq

RNAseq_Alignment_Pipeline_Summary.txt

ChIPseq_Alignment_Pipeline_Summary.txt

# ChIP-seq

## ChIPseq_PairedEnd_Bowtie2_SpikeIn_Alignment.pl

Used to align paired-end ChIP-seq reads against the concatenated genome (spike-in plus the genome of interest) for one sample. Produces two separate sorted and indexed bam files with uniquely mapped reads for the experimental and spike-in genomes, with PCR duplicates removed.

## ChIPseq_PairedEnd_Bowtie2_SpikeIn_AlignmentFromList.pl

Used to align paired-end ChIP-seq reads against the concatenated genome (spike-in plus the genome of interest) for a list of fastq files. Internally calls  ChIPseq_PairedEnd_Bowtie2_SpikeIn_Alignment.pl. 

## ChIPseq_RNAseq_Downsample_BamFile_FromList.pl

Used to subsample a certain fraction of reads for a list of bam files, with subsampling factors calculated based on the spike-in read normalisation as described in Fursova et al., 2019. 

## ChIPseq_bam2bigwig_MACS2.pl

Generates bigwig tracks from bam files using a MACS2 pileup function.

## ChIPseq_Calculate_Chromosome_Wide_Read_Coverage.R

Used to calculate chromosome-wide read coverage in 250 kb windows from bigwig files.

## ChIPseq_Plot_ChromosomeWide_Read_Coverage.R

Used to plot chromosome-wide read coverage calculated using the ChIPseq_Calculate_Chromosome_Wide_Read_Coverage.R function.

# RNA-seq

## RNAseq_PairedEnd_SpikeIn_Alignment.pl

Used to align paired-end RNA-seq reads against the concatenated genome (spike-in plus the genome of interest) for one sample. Produces two separate sorted and indexed bam files with uniquely mapped reads for the experimental and spike-in genomes, with PCR duplicates removed.

## RNAseq_PairedEnd_SpikeIn_AlignmentFromList.pl

Used to align paired-end RNA-seq reads against the concatenated genome (spike-in plus the genome of interest) for a list of fastq files. Internally calls  RNAseq_PairedEnd_SpikeIn_Alignment.pl. 

## ChIPseq_RNAseq_Downsample_BamFile_FromList.pl

Used to subsample a certain fraction of reads for a list of bam files, with subsampling factors calculated based on the spike-in read normalisation as described in Fursova et al., 2019. 

## RNAseq_PairedEnd_Split_Strand_Bigwig.pl

Used to split bam files into two separate files containing reads mapping either to the forward or reverse strand in order to generate strand-specific bigWig files.

## bed2gff.pl

Annotates BED3 or BED4 file using annotatePeaks.pl from HOMER and converts it to a 9 column .gff file.

## RNAseq_GFF_PairedEnd_RNA_Read_Count.pl

Used to count the number of mapped reads in the bodies of genes specified in a .gff file generated with bed2gff.pl

## GFF2table.pl

Converts a 9 column .gff file into a tab-delimited table for the subsequent analysis in R.

## RNAseq_DESeq2.R

Used to perform differential gene expression analysis, with spike-in based normalisation incorporated into a DESeq2 pipeline.

## Running_DESeq2_analysis.R

Used to call RNAseq_DESeq2.R script to run it for multiple cell lines.

