#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;

pod2usage("\nAligns paired-end FASTQ data with bowtie2. Uniquely mapped reads are sorted, indexed and PCR duplicates are removed with sambamba. \nUsage: -fastq1 <FASTQ> -fastq2 <FASTQ> -genome  <mm10, hg19> -spikein <genome of spiked-in normalisation control> -name <filename>\n") if (($#ARGV<0) && (-t STDIN));

&GetOptions ("fastq1=s"=>\my $fastq1,
	     "fastq2=s"=>\my $fastq2,
	     "genome=s"=> \my $genome,
	     "spikein=s"=> \my $spikegenome,
	     "name=s"=> \my $filename,
	     );


print "\nAligning to concatenated genome file...\n";

`bowtie2 -p 56 --no-mixed --no-discordant -x /databank/bowtie2/$genome.$spikegenome/$genome.$spikegenome -1 $fastq1 -2 $fastq2 | grep -v XS: - | samtools view -bhS -F4 - > $filename\_UniqMapped.bam`;
`sambamba sort --tmpdir /data/tmp/ -t 56 -m 60G -o $filename\_UniqMapped_sorted.bam $filename\_UniqMapped.bam`;
`sambamba markdup --tmpdir /data/tmp/ -r -t 56 $filename\_UniqMapped_sorted.bam $filename\_UniqMapped_sorted_rmdup.bam`;

print "\nExtracting reads aligning uniquely to $genome.\n";
`samtools view -h $filename\_UniqMapped_sorted_rmdup.bam | grep -v $spikegenome | sed s/$genome\_chr/chr/g | samtools view -bhS - > $filename\_$genome.UniqMapped_sorted_rmdup.bam`;
print "\nExtracting reads aligning uniquely to $spikegenome.\n";
`samtools view -h $filename\_UniqMapped_sorted_rmdup.bam | grep -v $genome | sed s/$spikegenome\_chr/chr/g | samtools view -bhS - > $filename\_$spikegenome.UniqMapped_sorted_rmdup.bam`;


`rm $filename\_UniqMapped.bam`;
`rm $filename\_UniqMapped_sorted.bam`;

