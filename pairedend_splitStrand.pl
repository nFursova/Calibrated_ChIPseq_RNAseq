#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;

pod2usage("\nSplits paired-end RNA seq file into strand-specific files.\nUsage: -bam <BAM file> -genome <eg. mm10>\n") if (($#ARGV<0) && (-t STDIN));


&GetOptions ("bam=s"=> \my $file,
	"genome=s"=> \my $genome,
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

`genomeCoverageBed -bga -split -ibam $name.bam -g /databank/chrom.sizes/$genome.chrom.sizes.txt | wigToBigWig stdin /databank/chrom.sizes/$genome.chrom.sizes.txt $name.bw`;
`genomeCoverageBed -bga -split -ibam $name\_FORWARD.bam -g /databank/chrom.sizes/$genome.chrom.sizes.txt | wigToBigWig stdin /databank/chrom.sizes/$genome.chrom.sizes.txt $name\_FORWARD.bw`;
`genomeCoverageBed -bga -split -ibam $name\_REVERSE.bam -g /databank/chrom.sizes/$genome.chrom.sizes.txt | wigToBigWig stdin /databank/chrom.sizes/$genome.chrom.sizes.txt $name\_REVERSE.bw`;

`rm $name\_83.bam`;
`rm $name\_163.bam`;
`rm $name\_147.bam`;
`rm $name\_99.bam`;


exit;
