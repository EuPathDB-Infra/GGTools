#!/usr/bin/perl

if(@ARGV < 2) {
    die "
Usage fixheader.pl <filename> <header>

Where <filename> is the name of a file and <header> is what you want
to change the first line of <filename> to be.  

This script outputs to standard out.

";
}

print "$ARGV[1]\n";
open(INFILE, $ARGV[0]);
$line = <INFILE>;
while($line = <INFILE>) {
    print $line;
}
close(INFILE);
