#!/usr/bin/perl

open(INFILE, $ARGV[0]);
$line = <INFILE>;
chomp($line);
@a = split(/\t/,$line);
$chr_prev = $a[0];
$end_prev = $a[2];
$value_prev = $a[3];
print "$chr_prev\t$a[1]";
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    if(!($a[0] eq $chr_prev && $a[1] == $end_prev && $a[3] == $value_prev)) {
	print "\t$end_prev\t$value_prev\n";
	print "$a[0]\t$a[1]";
    }
    $chr_prev = $a[0];
    $end_prev = $a[2];
    $value_prev = $a[3];
}
print "\t$end_prev\t$value_prev\n";
close(INFILE);
