#!/usr/bin/perl

# Written by Gregory R Grant
# University of Pennsylvania, 2010

if(@ARGV < 2) {
    die "
Usage: get_pbase-tbase.pl <probename2sequence> <seq col>

 This gets the pbase-tbase file needed by create_cdf.pl, this only works for match only arrays.
 <probename2sequence> is the input file mapping probe names to probe sequences
 <seq col> is the column number that has the probe sequence (start counting at zero).
 Note: it is assumed the first column (column zero) has the probe name.

This script outputs to std out.

";
}


$compliment{"A"} = "T";
$compliment{"T"} = "A";
$compliment{"C"} = "G";
$compliment{"G"} = "C";

open(INFILE, $ARGV[0]);
$line = <INFILE>;
$line = <INFILE>;
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    @b = split(//,$a[$ARGV[1]]);
    print "$a[0]";
    print "\t";
    print "$b[12]";
    print "\t";
    print $compliment{$b[12]};
    print "\n";
}
