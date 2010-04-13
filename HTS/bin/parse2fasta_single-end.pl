#/usr/bin/perl

# Written by Gregory R Grant
# University of Pennsylvania, 2010

use strict;

if(@ARGV < 1) {
    print "\nUsage: parse2fasta.pl <infile>\n\n";
    print "<infile> can be one file or a list of files separated by spaces.\n\n";
    print "PURPOSE: ";
    print "This program reformats files of reads into the appropriate fasta\nformat needed for the RUM pipeline.  If there are multiple files specified,\nthey are merged into one fasta file with consecutive sequence numbers.\n\n";
    print "NOTE: This script is for single-end reads only, for paired-end reads use\nparse2fasta_paried-end.pl\n\n";
    print "INPUT: ";
    print "Files can be fastq files, or fasta files, or more genreally the input\nfiles should have blocks of N rows for each read, where the read is on the\nsame row of each block.  N can be any positive integer and does not need to\nbe specified.\n\n";
    print "OUTPUT: ";
    print "Output is written to standard out, you should redirect it to a file.\n\n";

    exit();
}

$|=1;
open(INFILE1, $ARGV[0]);
my $cnt = 1;
my $firstNArow = 0;
my $secondNArow = 0;
while(1==1) {  # this loop figures out how many rows per block and which row the sequence is on.
    my $line = <INFILE1>;
    if($line eq '') {
	last;
    }
    chomp($line);
    $line =~ s/\^M$//;
    $line =~ s/[^ACGTN]$//;
    if($line =~ /^(A|C|G|T|N){10}(A|C|G|T|N)+$/) {  # find a line that looks like a line of seq,
                                                    # it's all A's, C's, G's, T's and N's and at
                                                    # least 11 chars long
	$firstNArow = $cnt;
	$line = <INFILE1>;
	chomp($line);
	$line =~ s/\^M$//;       #  This is to get rid of that pesky ^M character if it's there
	$line =~ s/[^ACGTN]$//;  #  In case any other weird charcter at the end, this line gets rid of it
	$cnt++;
	until($line =~ /^(A|C|G|T|N){10}(A|C|G|T|N)+$/) {  # find the next line of seq
	    $line = <INFILE1>;
	    chomp($line);
	    $line =~ s/\^M$//;
	    $line =~ s/[^ACGTN]$//;
	    if($line eq '') {
		last;
	    }
	    $cnt++;
	}
	$secondNArow = $cnt;
	last;
    }
    $cnt++;
}
close(INFILE1);

if($firstNArow == 0) {
    print "\nThis does not appear to be a valid file.\n\n";
    exit();
}

my $block = $secondNArow - $firstNArow;

# The number of rows in each block of rows is $block
# The number of the row in each block that has the sequence is $firstNArow

my $n = @ARGV;
for(my $i=0; $i<$n; $i++) { # loop over all the input files
    open(INFILE, $ARGV[$i]);
    my $cnt = 0;
    my $cnt2 = 1;
    my $line;
    my $linecnt = 0;
    while($line = <INFILE>) {    # this loop writes out the fasta file for file $i
	$linecnt++;
	$cnt++;
	if((($cnt - $firstNArow) % $block) == 0) {
	    print ">seq";
	    print ".$cnt2";
	    print "a\n";
	    chomp($line);
	    my $line_hold = $line;
	    $line =~ s/\^M$//;
	    $line =~ s/[^ACGTN]+$//;
	    if($line =~ /[^ACGTN]/) {
		print STDERR "\nERROR: There's something wrong with line $linecnt in file $ARGV[$i]\nIt should be a line of sequence but it is:\n$line_hold\n\n";
		exit();
	    }
	    print "$line\n";
	    $cnt2++;
	}
    }
    close(INFILE);
}
