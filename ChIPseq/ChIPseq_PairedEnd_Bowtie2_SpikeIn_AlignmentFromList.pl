#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;


pod2usage("\nPerforms alignment pipeline on a list of paired-end fastq files. \nUsage: -file <Tab separated file with FASTQ1, FASTQ2 and filenames eg. file_R1.fastq     file_R2.fastq		filename> -genome <mm10,hg19> -spikein <SpikeIn genome>\n") if (($#ARGV<0) && (-t STDIN));

&GetOptions ("file=s"=> \my $file,
             "genome=s"=> \my $genome,
             "spikein=s"=> \my $spikegenome,
             );

open(IN, $file) or die "Could not open $file";

#Specify location of the ChIPseq_PairedEnd_Bowtie2_SpikeIn_Alignment.pl script if not in PATH

while(<IN>){
            chomp;
            my ($fastq1,$fastq2,$name)=split (/\t/,);
            print "$name\n";
            my $command = `ChIPseq_PairedEnd_Bowtie2_SpikeIn_Alignment.pl -fastq1 $fastq1 -fastq2 $fastq2 -genome $genome -spikein $spikegenome -name $name`;

#indexing bam files

`sambamba index -t 56 $name\_UniqMapped_sorted_rmdup.bam`;
`sambamba index -t 56 $name\_$genome.UniqMapped_sorted_rmdup.bam`;
`sambamba index -t 56 $name\_$spikegenome.UniqMapped_sorted_rmdup.bam`;

#counting reads in bam files

my $count = `sambamba view -c -t 56 $name\_UniqMapped_sorted_rmdup.bam`;
print "\n";
print "\tNumber of total uniquely aligning reads is $count";

my $genomecount = `sambamba view -c -t 50 $name\_$genome.UniqMapped_sorted_rmdup.bam`;
print "\tNumber of reads aligning to $genome is $genomecount";

my $spikegenomecount = `sambamba view -c -t 50 $name\_$spikegenome.UniqMapped_sorted_rmdup.bam`;
print "\tNumber of reads aligning to $spikegenome is $spikegenomecount\n";


}

close IN;
