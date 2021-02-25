#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;


pod2usage("\nPerforms alignment pipeline on list of paired-end fastq files. \nUsage: -file <Tab separated file with FASTQ and filenames eg. file.fastq     filename> -genome <mm9,hg19> -spikein <SpikeIn genome>\n") if (($#ARGV<0) && (-t STDIN));

&GetOptions ("file=s"=> \my $file,
             "genome=s"=> \my $genome,
             "spikein=s"=> \my $spikegenome,
             );

my $outfile = "AlignmentStatistics.txt";

open(IN, $file) or die "Could not open $file";
open(OUT, ">$outfile") or die "Could not open $outfile";


while(<IN>){
            chomp;
            my ($fastq1,$fastq2,$name)=split (/\t/,);
            print OUT "$name\n";
            my $command = `RNAseq_PairedEnd.SpikeIn_Alignment_new.pl -fastq1 $fastq1 -fastq2 $fastq2 -genome $genome -spikein $spikegenome -name $name`;

my $count = `sambamba view -t 56 -c $name\_mapped_sorted_rmdup.bam`;
print OUT "\n";
print OUT "\tNumber of total uniquely aligning reads is $count";

my $genomecount = `sambamba view -t 56 -c $name\_$genome\_mapped_sorted_rmdup.bam`;
print OUT "\tNumber of reads aligning to $genome is $genomecount";

my $spikegenomecount = `sambamba view -t 56 -c $name\_$spikegenome\_mapped_sorted_rmdup.bam`;
print OUT "\tNumber of reads aligning to $spikegenome is $spikegenomecount\n\n";

}

close IN;

