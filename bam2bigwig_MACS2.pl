#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;

pod2usage("$0: usage is : -bam <bam file> -readtype <single or paired> -genome <genome> -name <name for file>") if (($#ARGV<0) && (-t STDIN));

&GetOptions ("bam=s"=> \my $bamfile,
             "genome=s"=> \my $genome,
             "readtype=s"=> \my $readtype,
             "name=s"=> \my $filename,
	     );

if ($readtype =~ "Single" || $readtype =~ "single") {
	`macs2 pileup -i $bamfile -f BAM -o $filename.bg`;
	`wigToBigWig -clip $filename.bg /databank/chrom.sizes/$genome.chrom.sizes.txt $filename.bw`;
	`rm $filename.bg`;
}
elsif ($readtype =~ "Paired" || $readtype =~ "paired") {
	`macs2 callpeak -t $bamfile -f BAMPE -g mm --bdg -n $filename --tempdir /data/tmp/`;
	`wigToBigWig -clip $filename\_treat_pileup.bdg /databank/chrom.sizes/$genome.chrom.sizes.txt $filename.bw`;
	`rm $filename\_control_lambda.bdg`;
	`rm $filename\_peaks.narrowPeak`;
	`rm $filename\_peaks.xls`;
	`rm $filename\_summits.bed`;
	`rm $filename\_treat_pileup.bdg`;
}






