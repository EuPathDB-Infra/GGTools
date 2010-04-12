# Written by Gregory R Grant
# University of Pennsylvania, 2010

use strict;
$|=1;

if(@ARGV != 2) {
    print STDERR "\nUsage: parse2fasta_paired-end.pl <foward reads file> <reverse reads file>\n\n";
    print "PURPOSE: ";
    print "This program reformats files of reads into the appropriate fasta format\nneeded for the RUM pipeline.  One file is output with the forward reads indicated\nby an 'a' and the reverse by a 'b'.\n\n";
    print "NOTE: This script is for paired-end reads only, for single-end reads use\nparse2fasta_single-end.pl\n\n";
    print "INPUT: ";
    print "Files can be fastq files, or fasta files, or more genreally the input\nfiles should have blocks of N rows for each read, where the read is on the\nsame row of each block.  N can be any positive integer and does not need to\nbe specified.\n\n";
    print "OUTPUT: ";
    print "Output is written to standard out, you should redirect it to a file.\n\n";
    exit(0);
}

open(INFILE1, $ARGV[0]);
my $cnt = 1;
my $firstNArow = 0;
my $secondNArow = 0;
while(1==1) {  # this loop figures out how many rows per block and which row the sequence is on.
    my $line = <INFILE1>;
    chomp($line);
    $line =~ s/\^M$//;       #  This is to get rid of that pesky ^M character if it's there
    $line =~ s/[^ACGTN]$//;  #  In case any other weird charcter at the end, this line gets rid of it
    if($line =~ /^(A|C|G|T|N){10}(A|C|G|T|N)+$/) {  # find a line that looks like a line of seq,
                                                    # it's all A's, C's, G's, T's and N's and at
                                                    # least 11 chars long
	$firstNArow = $cnt;
	$line = <INFILE1>;
	chomp($line);
	$line =~ s/\^M$//;
	$line =~ s/[^ACGTN]$//;
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

if($firstNArow == 0) {
    print "\nThis does not appear to be a valid file.\n\n";
    exit();
}

my $block = $secondNArow - $firstNArow;

# The number of rows in each block of rows is $block
# The number of the row in each block that has the sequence is $firstNArow

close(INFILE1);
open(INFILE1, $ARGV[0]);
open(INFILE2, $ARGV[1]);
$cnt = 0;
my $cnt2 = 1;
while(my $line1 = <INFILE1>) {    # this loop writes out the fasta file for file $i
    my $line2 = <INFILE2>;
    $cnt++;
    if((($cnt - $firstNArow) % $block) == 0) {
	print ">seq";
	print ".$cnt2";
	print "a\n";
	chomp($line1);
	$line1 =~ s/\^M$//;
	$line1 =~ s/[^ACGTN]+$//;
	print "$line1\n";
	print ">seq";
	print ".$cnt2";
	print "b\n";
	chomp($line2);
	$line2 =~ s/\^M$//;
	$line2 =~ s/[^ACGTN]+$//;
	print "$line2\n";
	$cnt2++;
    }
}
close(INFILE);
