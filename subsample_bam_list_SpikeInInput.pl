#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;

pod2usage("\nNormalises list of bam files by dm6/mm10 ratio in Input. Requires a list with two tab-separated fields: filename	fraction.  \nUsage: -file <Tab separated file with filenames and downsample fraction>>\n") if (($#ARGV<0) && (-t STDIN));


&GetOptions ("file=s"=> \my $file,
             );

open(IN, $file) or die "Could not open $file";


while(<IN>){
            chomp;
            my ($sample, $fraction)=split (/\t/,);

                my $bam = $sample;

            my $name = $bam;
                $name =~ s/.bam//;

        print "$bam\t$fraction\n";

`sambamba view -h -t 56 -f bam --subsampling-seed=123 -s $fraction $bam -o $name\_normalisedSpikeInInput.bam`;
`sambamba index -t 56 $name\_normalisedSpikeInInput.bam`;

my $count = `sambamba view -c -t 56 $name\_normalisedSpikeInInput.bam`;
print "\n";
print "Read Count after downsampling is $count\n";

}

close IN;

