#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;


pod2usage("\nPerforms alignment pipeline on list of paired-end fastq files. \nUsage: -file <Tab separated file with FASTQ and filenames eg. file.fastq     filename> -genome <mm9,hg19>\n") if (($#ARGV<0) && (-t STDIN));

&GetOptions ("file=s"=> \my $file,
             "genome=s"=> \my $genome,
             );

open(IN, $file) or die "Could not open $file";

#Specify location of the ATACseq_PairedEnd_Bowtie2_SpikeIn_Alignment.pl script if not in PATH

while(<IN>){
            chomp;
            my ($fastq1,$fastq2,$name)=split (/\t/,);
            print "$name\n";
            my $command = `ATACseq_PairedEnd_Bowtie2_SpikeIn_Alignment.pl -fastq1 $fastq1 -fastq2 $fastq2 -genome $genome -name $name`;

my $count = `samtools view -c $name\_mapped_sorted_rmdup_filtered.bam`;

print "\n";
print "Final read count (of full fragments) is $count\n";

}

close IN;
