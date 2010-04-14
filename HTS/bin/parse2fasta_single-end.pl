# Written by Gregory R Grant
# University of Pennsylvania, 2010

use strict;

if(@ARGV < 1) {
    print "\nUsage: parse2fasta.pl <infile> [options]\n\n";
    print "<infile> can be one file or a list of files separated by spaces.\n\n";
    print "Options:\n";
    print "     -firstrow n  : n is the first row that has sequence on it.\n";
    print "     -secondrow n : n is the second row that has sequence on it.\n";
    print "PURPOSE: ";
    print "This program reformats files of reads into the appropriate fasta\nformat needed for the RUM pipeline.  If there are multiple files specified,\nthey are merged into one fasta file with consecutive sequence numbers.\n\n";
    print "NOTE: This script is for single-end reads only, for paired-end reads use\nparse2fasta_paried-end.pl\n\n";
    print "INPUT: ";
    print "Files can be fastq files, or fasta files, or more genreally the input\nfiles should have blocks of N rows for each read, where the read is on the\nsame row of each block.  N can be any positive integer and does not need to\nbe specified.\n\n";
    print "OUTPUT: ";
    print "Output is written to standard out, you should redirect it to a file.\n\n";

    exit();
}

my $firstNArow = 0;
my $secondNArow = 0;
my $userparamsgiven = 0;
for(my $i=1; $i<@ARGV; $i++) {
    my $optionrecognized = 1;
    if($ARGV[$i] eq "-firstrow") {
	$firstNArow = $ARGV[$i+1] - 1;
	$optionrecognized = 0;
	$i++;
	$userparamsgiven = 1;
    }
    if($ARGV[$i] eq "-secondrow") {
	$secondNArow = $ARGV[$i+1] - 1;
	$optionrecognized = 0;
	$i++;
    }
    if($optionrecognized == 1) {
	die "\nERROR: option '$ARGV[$i-1]' not recognized.  Must be 'single' or 'paired'.\n";
    }
}
if(($firstNArow =~ /\S/ && !($secondNArow =~ /\S/)) && ($secondNArow =~ /\S/ && !($firstNArow =~ /\S/))) {
    die "\nERROR: you must set *both* -firstrow and -secondrow, or neither\n";
}

$|=1;
open(INFILE1, $ARGV[0]);
my $cnt = 0;
my @linearray;
my $line;
if($userparamsgiven == 0) {
    while(1==1) {  # this loop figures out how many rows per block and which row the sequence is on.
	$line = <INFILE1>;
	if($line eq '') {
	    last;
	}
	chomp($line);
	$line =~ s/\^M$//;
	$line =~ s/[^ACGTN]$//;
	if($line =~ /^(A|C|G|T|N){10}(A|C|G|T|N)+$/) {
	    $linearray[$cnt] = 1;
	} else {
	    $linearray[$cnt] = 0;
	}
	$cnt++;
	if($cnt > 20000) {
	    last;
	}
    }
    close(INFILE1);
    
    my $k;
    my $i;
    my $j;
    my $flag;
    for($k=0; $k<10; $k++) {
	for($i=1; $i<10; $i++) {
	    $flag = 0;
	    for($j=0; $j<$cnt/20; $j++) {
		my $x = $k+$i*$j;
#	    print "k=$k, i=$i, j=$j\n";
#	    print "linearray[$x] = $linearray[$x]\n";
		if($linearray[$k+$i*$j] == 0) {
		    $flag = 1;
		}
	    }
	    if($flag == 0) {
		$firstNArow = $k;
		$secondNArow = $k+$i;
		$k=10;
		$i=10;
	    }
	}
	if($k==9 && $flag == 0) {
	    die "\nError: canont determine which lines have the sequence.\n\n";
	}
    }
}
#print "firstNArow = $firstNArow\n";
#print "secondNArow = $secondNArow\n";

if($firstNArow == 0) {
    print "\nThis does not appear to be a valid file.\n\n";
    exit();
}

my $block = $secondNArow - $firstNArow;

# The number of rows in each block of rows is $block
# The number of the row in each block that has the sequence is $firstNArow

my $n = @ARGV;
my $linecnt = 0;
for(my $i=0; $i<$n; $i++) { # loop over all the input files
    open(INFILE, $ARGV[$i]);
    $cnt = 0;
    my $cnt2 = 1;
    $line;
    $linecnt = 0;
    while($line = <INFILE>) {    # this loop writes out the fasta file for file $i
	$linecnt++;
	if((($cnt - $firstNArow) % $block) == 0) {
	    print ">seq";
	    print ".$cnt2";
	    print "a\n";
	    chomp($line);
	    my $line_hold = $line;
	    $line =~ s/\^M$//;
	    $line =~ s/[^ACGTN]+$//;
	    if($line =~ /[^ACGTN]/ || !($line =~ /\S/)) {
		print STDERR "\nERROR: There's something wrong with line $linecnt in file $ARGV[$i]\nIt should be a line of sequence but it is:\n$line_hold\n\n";
		exit();
	    }
	    print "$line\n";
	    $cnt2++;
	}
	$cnt++;
    }
    close(INFILE);
    if($linecnt % $block != 0) {
	print STDERR "\nWarning: the last block of lines in file $ARGV[$i] is not the right size.\n\n";
    }
}
