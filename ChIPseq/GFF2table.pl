#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;

pod2usage("\nConvert GFF file into tab separated file. \nUsage: -gff <gff file>\n") if (($#ARGV<0) && (-t STDIN));

&GetOptions ("gff=s"=> \my $gff_file,
	     );
             
my $filename = $gff_file;
$filename =~ s/.gff//;

my $outfile = "$filename.txt";

open(IN, $gff_file) or die "Could not open $gff_file";
open(OUT, ">$outfile") or die "Could not open $outfile";


my $firstline = <IN>; 
            chomp $firstline;

            my ($chr, $mig, $peak, $start, $stop, $col6, $strand, $col8, $col9)= split /[\t]+/,$firstline;
            print OUT "Chr\tStart\tStop\tStrand\t";
            
            my @array = split (/[\;]+/,$col9);
            my $length = @array;
            
            foreach(my $i = 0; $i <= $length-1; $i += 1){
                        my $arrayindex = $array[$i];
                        my ($factor,$value) = split (/[\=]+/,$arrayindex);
                        print OUT "$factor\t"
            }
            print OUT "\n";




while(<IN>){
            chomp;

            my ($chr, $mig, $peak, $start, $stop, $col6, $strand, $col8, $col9)= split /[\t]+/,;
            
            print OUT "$chr\t$start\t$stop\t$strand\t";
            
            my @array = split (/[\;]+/,$col9);
            my $length = @array;
            
            foreach(my $i = 0; $i <= $length-1; $i += 1){
                        my $arrayindex = $array[$i];
                        my ($factor,$value) = split (/[\=]+/,$arrayindex);
                        print OUT "$value\t"
            }
            print OUT "\n";
}

close IN;
close OUT;

