#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;

pod2usage("\nAnnotates GFF file with raw read counts or RPKM from a list of bam files. \nUsage: -gff <9column GFF file> -file <Tab delimited file with list of .bam files and labels> -genome <mm10,hg19> -rpkm <True/False - whether to normalise read count to reads per million per kilobase> \n") if (($#ARGV<0) && (-t STDIN));

&GetOptions ("file=s"=> \my $file,
             "genome=s"=> \my $genome,
             "gff=s"=> \my $gff,
			 "rpkm=s"=> \my $rpkmstatement,
	          );

my $gffname = $gff;
   $gffname =~ s/.gff//;           
my $input = "$gffname.ANNOTATED.gff";   
   `cp $gff $input`;
   

open(IN, $file) or die "Could not open $file";

while(<IN>){
chomp;
my ($bamfile, $label)=split (/\t/,);
print "$label\n";

#my $outfile = "$gffname\.tmp.$label\.gff";
my $outfile = "$genome\.$label\.gff";

open(INPUT, $input) or die "Could not open $input";
open(OUT, ">$outfile") or die "Could not open $outfile";


	while (<INPUT>){
                    chomp;
                        my ($col1, $col2, $col3, $col4, $col5, $col6, $col7, $col8, $col9)=split (/\t+/);
                        my $interval = "$col1:$col4-$col5"; 
                            my $sizekb = ($col5 - $col4)/1000;

                if ($col7 =~ /[+]/){
                    my $commandinterval83 = `samtools view -c -f 83 $bamfile $interval` or die "Could not perform interval read count with samtool";
                    my $commandinterval163 = `samtools view -c -f 163 $bamfile $interval` or die "Could not perform interval read count with samtool";
                    my $normalisedintervalreads = (($commandinterval83+$commandinterval163));
                   print OUT "$col1\t$col2\t$col3\t$col4\t$col5\t$col6\t$col7\t$col8\t$col9;$label\_ReadCount=$normalisedintervalreads\n";
                }
                elsif ($col7 =~ /[-]/){
                    my $commandinterval99 = `samtools view -c -f 99 $bamfile $interval` or die "Could not perform interval read count with samtool";
                    my $commandinterval147 = `samtools view -c -f 147 $bamfile $interval` or die "Could not perform interval read count with samtool";
                    my $normalisedintervalreads = (($commandinterval99+$commandinterval147));
                   print OUT "$col1\t$col2\t$col3\t$col4\t$col5\t$col6\t$col7\t$col8\t$col9;$label\_ReadCount=$normalisedintervalreads\n";
                }
	}						
	`mv $outfile $input`;
	close OUT;


close INPUT;

}

close IN;




