#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;

pod2usage("\nSplits paired-end RNA seq file into strand-specific files and generates bigwigs. \nRequires samtools, sambamba, bedtools, and wigToBigWig from UCSC in PATH. \nUsage: -bam <BAM file> \n") if (($#ARGV<0) && (-t STDIN));

&GetOptions ("bam=s"=> \my $file,
             );


            my $name = $file;
                $name =~ s/.bam//;

##Forward reads
`samtools view -b -f 83 $file > $name\_83.bam`;
`samtools view -b -f 163 $file > $name\_163.bam`;

##Reverse reads
`samtools view -b -f 147 $file > $name\_147.bam`;
`samtools view -b -f 99 $file > $name\_99.bam`;


##Merge reads
`sambamba merge -t 56 $name\_FORWARD.bam $name\_83.bam $name\_163.bam`;
`sambamba merge -t 56 $name\_REVERSE.bam $name\_147.bam $name\_99.bam`;

`sambamba index -t 56 $name\_FORWARD.bam`;
`sambamba index -t 56 $name\_REVERSE.bam`;

#Specify location of YOUR LOCAL chromosome size file that can be downloaded from https://hgdownload-test.gi.ucsc.edu/goldenPath/$YOUR_GENOME_OF_INTEREST/biZips/.

`genomeCoverageBed -bga -split -ibam $name.bam -g chrom.sizes.txt | wigToBigWig stdin chrom.sizes.txt $name.bw`;
`genomeCoverageBed -bga -split -ibam $name\_FORWARD.bam -g chrom.sizes.txt | wigToBigWig stdin chrom.sizes.txt $name\_FORWARD.bw`;
`genomeCoverageBed -bga -split -ibam $name\_REVERSE.bam -g chrom.sizes.txt | wigToBigWig stdin chrom.sizes.txt $name\_REVERSE.bw`;

`rm $name\_83.bam`;
`rm $name\_163.bam`;
`rm $name\_147.bam`;
`rm $name\_99.bam`;


exit;
