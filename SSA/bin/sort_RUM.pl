#!/usr/bin/perl

# Written by Gregory R. Grant
# Universiity of Pennsylvania, 2010

if(@ARGV < 2) {
    die "
Usage: sort_RUM.pl <RUM file> <sorted outfile>

This script sorts a RUM output file by sequence number.  It keeps
consistent pairs together.
";
}

$|=1;
open(INFILE, $ARGV[0]);
$seqnum_prev = 0;
open(OUTFILE1, ">RUM_sorted_temp1");
open(OUTFILE2, ">RUM_unsorted_temp1");
while($line = <INFILE>) {
    $line =~ /^seq.(\d+)/;
    $seqnum = $1;
    if($seqnum >= $seqnum_prev) {
	print OUTFILE1 $line;
	$seqnum_prev = $seqnum;
    } else {
	print OUTFILE2 $line;
	$still_unsorted_flag = 1;
    }
}
close(OUTFILE1);
close(OUTFILE2);
close(INFILE);

$still_unsorted_flag = 1;
while($still_unsorted_flag == 1) {
    $still_unsorted_flag = 0;
    open(INFILE, "RUM_unsorted_temp1");
    $seqnum_prev = 0;
    open(OUTFILE1, ">RUM_sorted_temp2");
    open(OUTFILE2, ">RUM_unsorted_temp2");
    while($line = <INFILE>) {
	$line =~ /^seq.(\d+)/;
	$seqnum = $1;
	if($seqnum >= $seqnum_prev) {
	    print OUTFILE1 $line;
	    $seqnum_prev = $seqnum;
	} else {
	    print OUTFILE2 $line;
	    $still_unsorted_flag = 1;
	}
    }
    close(OUTFILE1);
    close(OUTFILE2);
    close(INFILE);
    `mv RUM_unsorted_temp2 RUM_unsorted_temp1`;
    merge();
    print "done merging...\n";
}
$sortedfile = $ARGV[1];
`mv RUM_sorted_temp1 $sortedfile`;
unlink("RUM_unsorted_temp1");

sub merge () {
    print STDERR "merging...\n";
    open(INFILE1, "RUM_sorted_temp1");
    open(INFILE2, "RUM_sorted_temp2");
    open(OUTFILE, ">RUM_sorted_temp3");
    $flag = 0;
    $line1 = <INFILE1>;
    chomp($line1);
    $line1 =~ /^seq.(\d+)/;
    $seqnum1 = $1;
    $line2 = <INFILE2>;
    chomp($line2);
    $line2 =~ /^seq.(\d+)/;
    $seqnum2 = $1;
    while($flag == 0) {
	while($seqnum1 <= $seqnum2 && $line1 ne '') {
	    print OUTFILE "$line1\n";
	    $line1 = <INFILE1>;
	    chomp($line1);
	    $line1 =~ /^seq.(\d+)/;
	    $seqnum1 = $1;
	    if($line1 eq '') {
		if($line2 =~ /\S/) {
		    chomp($line2);
		    print OUTFILE "$line2\n";
		}
		while($line2 = <INFILE2>) {
		    print OUTFILE $line2;		    
		}
	    }
	}
	while($seqnum2 <= $seqnum1 && $line2 ne '') {
	    print OUTFILE "$line2\n";
	    $line2 = <INFILE2>;
	    chomp($line2);
	    $line2 =~ /^seq.(\d+)/;
	    $seqnum2 = $1;
	    if($line2 eq '') {
		if($line1 =~ /\S/) {
		    chomp($line1);
		    print OUTFILE "$line1\n";
		}
		while($line1 = <INFILE1>) {
		    print OUTFILE $line1;
		}
	    }
	}
	if($line1 eq '' && $line2 eq '') {
	    $flag = 1;
	}
    }
    `mv RUM_sorted_temp3 RUM_sorted_temp1`;
    unlink("RUM_sorted_temp2");
}
