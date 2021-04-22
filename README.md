# Calibrated-ChIPseq-RNAseq-tools

Scripts developed in Klose Lab to perform analysis of spike-in calibrated ChIP-seq and RNA-seq data. Contributed by Dr Hamish King (h.king@qmul.ac.uk, drhamishking@gmail.com), Dr Nadezda Fursova (nadya.fursova@nih.gov, nfursova.msu@gmail.com), and Dr Anne Turberfield. Used for the data analysis in the manuscript "BAP1 constrains pervasive H2AK119ub1 to control the transcriptional potential of the genome".

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

# Summary of the read alignment and processing pipelines for calibrated RNA-seq and ChIP-seq
RNAseq_Alignment_Pipeline_Summary.txt

ChIPseq_Alignment_Pipeline_Summary.txt
