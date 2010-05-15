#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

$|=1;

if(@ARGV < 2) {
    die "
Usage: get_chr_from_cov.pl <cov file> <chr>

Output is to standard out, redirect to a file using \">\"

";
}

open(INFILE, $ARGV[0]);
# line of input (wig file) looks like this:
# chr9    3342879 3342987 1
$chr = $ARGV[1];
$line = <INFILE>;
print $line;
while($line = <INFILE>) {
    if($line =~ /^$chr\t/) {
	print $line;
    }
}
