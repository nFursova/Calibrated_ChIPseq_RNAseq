#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;

pod2usage("\nAnnotates bed file with HOMER and converts to nine column GFF file. \nUsage: -bed <chr,start,stop> -genome <mm10, hg19>\n") if (($#ARGV<0) && (-t STDIN));

&GetOptions ("bed=s"=> \my $bed_file,
             "genome=s"=> \my $genome,
	     );

my $filename = $bed_file;
   $filename =~ s/.bed//;
my $filepath = `fp $bed_file`;
   $filepath =~ s/.bed//;

my $homerfile = "$filename\_homer.txt";   
my $outfile = "$filename.tmp.gff";

my $homercommand = `annotatePeaks.pl $bed_file $genome > $homerfile`;

open(IN, $homerfile) or die "Could not open $homerfile";
open(OUT, ">$outfile") or die "Could not open $outfile";


while(my $line = <IN>){
            chomp;
            if ($line=~ /PeakID/) { next }                                          #skips header

            else{                      
            my ($peakid, $chr, $start, $stop, $strand, $peakscore, $focusratio, $annotation, $detail_annotation, $distancetoTSS, $nearestpromoterid, $entrezid, $unigeneid, $refseqid, $ensemblid, $gene_name, $gene_alias, $gene_description, $gene_type)=split (/\t+/,$line); #chomps, splits and assigns bed file into three columns
            my $distancetoTSS_value = abs($distancetoTSS);
		my $length = $stop-$start;
		my $start = $start-1;
            print OUT "$chr\tMIG\tPeak\t$start\t$stop\t.\t$strand\t.\tSize=$length;Gene_name=$gene_name;Distance_Nearest_TSS=$distancetoTSS_value;Nearest_Promoter_ID=$nearestpromoterid;Refseq_ID=$refseqid\n";
            }
}

close IN;
close OUT;

`sortBed -i $outfile > $filename.gff`;
`rm $outfile`;
exit

