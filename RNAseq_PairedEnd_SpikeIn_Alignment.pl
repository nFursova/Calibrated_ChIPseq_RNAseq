#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
use strict;

pod2usage("$0: usage is : -fastq1 <FASTQ file> -fastq2 <FASTQ file> -genome  <eg. mm9> -spikein <spikein genome eg. dm6> -name <filename>") if (($#ARGV<0) && (-t STDIN));

&GetOptions ("fastq1=s"=>\my $fastq1,
	     "fastq2=s"=>\my $fastq2,
	     "genome=s"=> \my $genome,
	     "spikein=s"=> \my $spikegenome,
	     "name=s"=> \my $filename,
	     );

my $outfile = "AlignmentStatistics.txt";
open(OUT, ">>$outfile") or die "Could not open $outfile";


##Align reads to rDNA, collect only unmapped reads and convert to fastq.
print "Aligning $filename to ribosomal DNA\n";
`bowtie2 --very-fast -p 56 --no-mixed --no-discordant -x /databank/bowtie2/$genome.$spikegenome.rDNA/$genome.$spikegenome.rDNA -1 $fastq1 -2 $fastq2 | samtools view -bSh -u -f 4 - | bamToFastq -i stdin -fq $filename\_1.fastq -fq2 $filename\_2.fastq`;

print "Finished processing ribosomal DNA\n";


##Align remaining reads with STAR
print "Aligning ribosomal-depleted reads with STAR\n";
	`STAR --genomeDir /databank/STAR/$genome.$spikegenome --readFilesIn $filename\_1.fastq $filename\_2.fastq --runThreadN 40 --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outFileNamePrefix $filename --genomeLoad LoadAndKeep`;
	`samtools view -bhS -q 4 $filename\Aligned.out.sam > $filename\_STARmapped.bam`;
        `sambamba sort -t 56 -m 75000000000 $filename\_STARmapped.bam -o $filename\_STARmapped_sorted.bam`;
print "Finished STAR alignment and filtering\n";
	`head -n 32 $filename\\Log.final.out >> $outfile`;

##Align remaining reads with bowtie2

print "Now aligning with bowtie2\n";
	`bowtie2 --sensitive-local -p 56 --no-mixed --no-discordant --mm -x /databank/bowtie2/$genome.$spikegenome/$genome.$spikegenome -1 $filename\\Unmapped.out.mate1 -2 $filename\\Unmapped.out.mate2 | grep -v XS: - | samtools view -bhS -F4 - > $filename\_Bowtie2mapped.bam`;
         `sambamba sort -t 56 -m 75000000000 $filename\_Bowtie2mapped.bam -o $filename\_Bowtie2mapped_sorted.bam`;
print "Finished bowtie2 alignment\n";



`sambamba merge -t 56  $filename\_mapped.bam $filename\_Bowtie2mapped_sorted.bam $filename\_STARmapped_sorted.bam`;

`sambamba sort -t 56 -m 75000000000 $filename\_mapped.bam -o $filename\_mapped_sorted.bam`;
`sambamba markdup -r -t 56 $filename\_mapped_sorted.bam $filename\_mapped_sorted_rmdup.bam`;

#

	print "\nExtracting reads aligning uniquely to $genome.\n";
`samtools view -@ 56 -h $filename\_mapped_sorted_rmdup.bam | grep -v $spikegenome | sed s/$genome\_chr/chr/g | samtools view -@ 56 -bhS - > $filename\_$genome\_mapped_sorted_rmdup.bam`;
	print "\nExtracting reads aligning uniquely to $spikegenome.\n";
`samtools view -@ 56 -h $filename\_mapped_sorted_rmdup.bam | grep -v $genome | sed s/$spikegenome\_chr/chr/g | samtools view -@ 56 -bhS - > $filename\_$spikegenome\_mapped_sorted_rmdup.bam`;


#indexing bam files

`sambamba index -t 56 $filename\_mapped_sorted_rmdup.bam`;
`sambamba index -t 56 $filename\_$genome\_mapped_sorted_rmdup.bam`;
`sambamba index -t 56 $filename\_$spikegenome\_mapped_sorted_rmdup.bam`;


`rm $filename\_1.fastq`;
`rm $filename\_2.fastq`;
`rm $filename\Aligned.out.sam`;
`rm $filename\\Unmapped.out.mate1`;
`rm $filename\\Unmapped.out.mate2`;
`rm $filename\\Log.progress.out`;
`rm $filename\SJ.out.tab`;


`rm $filename\_Bowtie2mapped.bam`;
`rm $filename\_Bowtie2mapped_sorted.bam`;

`rm $filename\_STARmapped.bam`;
`rm $filename\_STARmapped_sorted.bam`;

`rm $filename\_mapped.bam`;
`rm $filename\_mapped_sorted.bam`;

