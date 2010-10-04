#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennslyvania, 2010

if (@ARGV < 1) {
  die "
Usage: fastq2qualitystats.pl <fastq_file> > <output_file>

Where: <fastq_file> is the fastq sequence file.
This program prints to output the average quality score over all reads for each
read position.";
}

$string = "\@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh";

@a = split(//,$string);

for($i=0; $i<@a; $i++) {
    $qualitychar2score{$a[$i]} = $i + 64;
}
#foreach $char (sort {$qualitychar2score{$a} <=> $qualitychar2score{$b}} keys %qualitychar2score) {
#    print "$char: $qualitychar2score{$char}\n";
#}

open(INFILE, $ARGV[0]);
$num_seqs = 0;
while($line = <INFILE>) {
    $line = <INFILE>;
    $line = <INFILE>;
    $line = <INFILE>;
    chomp($line);
    @a = split(//,$line);
    $N = @a;
    for($i=0; $i<@a; $i++) {
	$sum_of_quals{$i} = $sum_of_quals{$i} + $qualitychar2score{$a[$i]};
    }
    $num_seqs++;
}
foreach $i (sort {$a <=> $b} keys %sum_of_quals) {
    $ave_of_quals{$i} = $sum_of_quals{$i} / $num_seqs;
}

foreach $i (sort {$a <=> $b} keys %ave_of_quals) {
    print "ave qual position $i: $ave_of_quals{$i}\n";
}
