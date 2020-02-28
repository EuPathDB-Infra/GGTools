# Written by Elisabetta Manduchi, based on a script by Greg Grant
# University of Pennslyvania, 2010

use strict;

if (@ARGV < 1) {
  die "
Usage: fastq2qualitystats.pl <fastq_file> > <output_file>

Where: <fastq_file> is the fastq sequence file.
This program prints to output the average quality score over all reads for each
read position.";
}

my $string = "\@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

my @a = split(//,$string);

my %qualitychar2score;
for(my $i=0; $i<@a; $i++) {
    $qualitychar2score{$a[$i]} = $i;
}
#foreach my $char (sort {$qualitychar2score{$a} <=> $qualitychar2score{$b}} keys %qualitychar2score) {
#    print "$char: $qualitychar2score{$char}\n";
#}

open(INFILE, $ARGV[0]);
my $num_seqs_total;
my @num_seqs;
my @sum_of_quals;
while(my $line = <INFILE>) {
    $line = <INFILE>;
    $line = <INFILE>;
    $line = <INFILE>;
#    print $line;
    chomp($line);
    my @a = split(//,$line);
    for(my $i=0; $i<@a; $i++) {
      $sum_of_quals[$i] += $qualitychar2score{$a[$i]};
      $num_seqs[$i] += 1;
    }
    $num_seqs_total++;
}

print "POSITION\tFRAC_READS\tAVG_Q\n";
for (my $i=1;$i<=@num_seqs;$i++) {
  printf("%s\t%.3f\t%.1f\n",
    $i,
    $num_seqs[$i-1] / $num_seqs_total,
    $sum_of_quals[$i-1] / $num_seqs[$i-1],
 );
}
