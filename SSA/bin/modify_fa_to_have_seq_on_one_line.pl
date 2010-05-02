#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

if(@ARGV < 1) {
    die "
Usage: modify_fa_to_have_seq_on_one_line.pl <fasta_file>

Where <fasta_file> is a fasta file.

This program returns a fasta file that is equivalent to the
input file, but has the sequences for each entry on one line.

This script outputs to standard out.

";
}

open(INFILE, $ARGV[0]);
$flag = 0;
while($line = <INFILE>) {
    if($line =~ />/) {
	if($flag == 0) {
	    print $line;	
	    $flag = 1;
	} else {
	    print "\n$line";
	}
    } else {
	chomp($line);
	print $line;
    }
}
close(INFILE);
