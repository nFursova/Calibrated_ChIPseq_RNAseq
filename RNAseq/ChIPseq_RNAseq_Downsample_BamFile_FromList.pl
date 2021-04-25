#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;

pod2usage("\nDownsamples a list of bam files. Requires a list with two tab-separated fields: filename	fraction. Filename should be a full file name.\nDesigned to be used with alignment_from_list.pl workflow.") if (($#ARGV<0) && (-t STDIN));


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

`sambamba view -h -t 56 -f bam --subsampling-seed=123 -s $fraction $bam -o $name\_downsampled.bam`;
`sambamba index -t 56 $name\_downsampled.bam`;

my $count = `sambamba view -c -t 56 $name\_downsampled.bam`;
print "\n";
print "Read Count after downsampling is $count\n";

}

close IN;
