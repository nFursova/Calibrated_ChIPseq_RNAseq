#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;

pod2usage("$0: usage is : -fastq1 <FASTQ file> -fastq2 <FASTQ file> -genome  <eg. mm9> -name <filename>") if (($#ARGV<0) && (-t STDIN));

&GetOptions ("fastq1=s"=>\my $fastq1,
             "fastq2=s"=>\my $fastq2,
	     "genome=s"=> \my $genome,
	     "name=s"=> \my $filename,
	     );

#align fastq file with bowtie2
`bowtie2 -p 20 --no-mixed --no-discordant -x /databank/bowtie2/$genome/$genome -1 $fastq1 -2 $fastq2 | grep -v XS: - | samtools view -bhS -F4 - > $filename\_mapped.bam`;

##Sort and remove duplicates 
`samtools sort $filename\_mapped.bam > $filename\_mapped_sorted.bam`;                                     
`samtools rmdup $filename\_mapped_sorted.bam $filename\_mapped_sorted_rmdup.bam`;                   

##filter file with ATAC.Filter.Regions
`intersectBed -v -abam $filename\_mapped_sorted_rmdup.bam -b /data/bioc1336/annotations/mm10.ATAC.Filter.Regions.MERGED.bed > $filename\_mapped_sorted_rmdup_filtered.bam`;
`samtools index $filename\_mapped_sorted_rmdup_filtered.bam $filename\_mapped_sorted_rmdup_filtered.bai`;


`rm $filename\_mapped.bam`;
`rm $filename\_mapped_sorted.bam`;
`rm $filename\_mapped_sorted_rmdup.bam`;


exit;
