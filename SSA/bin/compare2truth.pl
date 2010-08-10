#!/bin/perl

if(@ARGV < 2) { 
    die "
Usage: perl compare2truth.pl <truth file> <sam file>

The <truth file> is output from reads_simulator.pl

The <sam file> is the alignments of the same reads output from 
               an alignment algoirthm.

This script computes the per-base false alignment rates.

";
}

open(INFILE1, $ARGV[0]);
open(INFILE2, $ARGV[1]);

$head = `head -1 $ARGV[0]`;
chomp($head);
$head =~ /seq.(\d+)/;
$first_seq_num = $1;

$head =~ /\t([^\t]+)$/;
$readlength = length($1);
print "readlength = $readlength\n";

$tail = `tail -1 $ARGV[0]`;
$tail =~ /seq.(\d+)/;
$last_seq_num = $1;

$head = `head -1 $ARGV[1]`;
$head =~ /seq.(\d+)/;
if($first_seq_num ne $1) {
    die "
Error: both files must start and end with the same sequence number
and must have an entry for every sequence number in between.

";
}

$tail = `tail -1 $ARGV[0]`;
$tail =~ /seq.(\d+)/;
if($last_seq_num ne $1) {
    die "
Error: both files must start and end with the same sequence number
and must have an entry for every sequence number in between.

";
}

$linenum = 0;
$total_number_of_bases_of_reads = 0;
$total_number_of_bases_aligned_correctly = 0;
$total_number_of_bases_aligned_incorrectly = 0;
$total_number_of_bases_aligned_ambiguously = 0;
$total_number_of_bases_unaligned = 0;

$flag = 0;
for($seqnum=$first_seq_num; $seqnum<=$last_seq_num; $seqnum++) {
    $truth = <INFILE1>;
    chomp($truth);
    $total_number_of_bases_of_reads = $total_number_of_bases_of_reads + $readlength;
    $linenum++;
    $sam = <INFILE2>;
    chomp($sam);
    if($sam =~ /IH:i:(\d+)/) {
	$num_alignments = $1;
	if($num_alignments > 1) {
	    $total_number_of_bases_aligned_ambiguously = $total_number_of_bases_aligned_ambiguously + $readlength * 2;
	    for($i=0; $i<$num_alignments*2-1; $i++) {
		$sam = <INFILE2>;	    
	    }
	    $truth = <INFILE1>;	    
	    $total_number_of_bases_of_reads = $total_number_of_bases_of_reads + $readlength;
	    $linenum++;
	    next;
	}
    } else { # this means the entire read did not align
	$total_number_of_bases_unaligned = $total_number_of_bases_unaligned + $readlength;
	if($flag == 0) {
	    $seqnum--;
	    $flag = 1;
	} else {
	    $flag = 0;
	}
	next;
    }
    $truth =~ /seq.(\d+.)/;
    $sn1 = $1;
    $sam =~ /seq.(\d+.)/;
    $sn2 = $1;
    $sn = $sn1;
    $sn =~ s/.$//;
    if($sn1 ne $sn2 || $sn != $seqnum) {
	die "\tError: Check line $linenum.  It doesn't refer to the same sequence in both files\n\n";
    }
    @a = split(/\t/, $truth);
    $truth_chr = $a[1];
    $truth_start = $a[2];
    $truth_cigar = $a[3];
    @a = split(/\t/, $sam);
    $sam_chr = $a[2];
    $sam_start = $a[3];
    $sam_cigar = $a[5];

    if($sam_cigar eq "*") {
	# this means the entire read did not align
	$total_number_of_bases_unaligned = $total_number_of_bases_unaligned + $readlength;
	if($flag == 0) {
	    $seqnum--;
	    $flag = 1;
	} else {
	    $flag = 0;
	}
	next;
    }

    $cigar_string_temp = $sam_cigar;
    $number_bases_of_this_read_aligned = 0;
    while($cigar_string_temp =~ /^(\d+)([^\d])/) {
	$num = $1;
	$type = $2;
	if($type eq 'S') {
	    $total_number_of_bases_unaligned = $total_number_of_bases_unaligned + $num;
	}
	if($type eq 'M') {	
	    $number_bases_of_this_read_aligned = $number_bases_of_this_read_aligned + $num;
	}
	$cigar_string_temp =~ s/^\d+[^\d]//;
    }

    if($truth_chr ne $sam_chr) {
	$total_number_of_bases_aligned_incorrectly = $total_number_of_bases_aligned_incorrectly + $number_bases_of_this_read_aligned;
    } else {
	$pos_on_genome = $sam_start;
	$pos_on_read = 0;
	$cigar_string_temp = $sam_cigar;
	undef @location;
	undef @truelocation;
	while($cigar_string_temp =~ /^(\d+)([^\d])/) {
	    $num = $1;
	    $type = $2;
	    if($type eq 'S') {
		for($i=0; $i<$num; $i++) {
		    $location[$pos_on_read] = "x";
		    $pos_on_read++;
		}
	    }
	    if($type eq 'M') {
		for($i=0; $i<$num; $i++) {
		    $location[$pos_on_read] = $position_on_genome;
		    $pos_on_read++;
		    $pos_on_genome++;
		}
	    }
	    if($type eq 'I') {
		for($i=0; $i<$num; $i++) {
		    $location[$pos_on_read] = "i";
		    $pos_on_read++;
		}
	    }
	    $cigar_string_temp =~ s/^\d+[^\d]//;
	}
	$pos_on_genome = $truth_start;
	$pos_on_read = 0;
	$cigar_string_temp = $truth_cigar;
	while($cigar_string_temp =~ /^(\d+)([^\d])/) {
	    $num = $1;
	    $type = $2;
	    if($type eq 'M') {
		for($i=0; $i<$num; $i++) {
		    $truelocation[$pos_on_read] = $position_on_genome;
		    $pos_on_read++;
		    $pos_on_genome++;
		}
	    }
	    if($type eq 'I') {
		for($i=0; $i<$num; $i++) {
		    $truelocation[$pos_on_read] = "i";
		    $pos_on_read++;
		}
	    }
	    $cigar_string_temp =~ s/^\d+[^\d]//;
	}
	for($pos_on_read=0; $pos_on_read<@truelocation; $pos_on_read++) {
	    if($truelocation[$pos_on_read] eq $location[$pos_on_read]) {
		$total_number_of_bases_aligned_correctly++;
	    } else {
		if($location[$pos_on_read] ne "x") {
		    $total_number_of_bases_aligned_incorrectly++;
		}
	    }
	}
    }
    if($flag == 0) {
	$seqnum--;
	$flag = 1;
    } else {
	$flag = 0;
    }
}

print "total_number_of_bases_of_reads = $total_number_of_bases_of_reads\n";

$percent_bases_aligned_correctly = int($total_number_of_bases_aligned_correctly / $total_number_of_bases_of_reads * 10000) / 100;;
#print "total_number_of_bases_aligned_correctly = $total_number_of_bases_aligned_correctly\n";
print "% bases aligned correctly: $percent_bases_aligned_correctly%\n";

$percent_bases_aligned_incorrectly = int($total_number_of_bases_aligned_incorrectly / $total_number_of_bases_of_reads * 10000) / 100;;
#print "total_number_of_bases_aligned_incorrectly = $total_number_of_bases_aligned_incorrectly\n";
print "% bases aligned incorrectly: $percent_bases_aligned_incorrectly%\n";

$percent_bases_aligned_ambiguously = int($total_number_of_bases_aligned_ambiguously / $total_number_of_bases_of_reads * 10000) / 100;
#print "total_number_of_bases_aligned_ambiguously = $total_number_of_bases_aligned_ambiguously\n";
print "% bases aligned ambiguously: $percent_bases_aligned_ambiguously%\n";

$percent_bases_unaligned = int($total_number_of_bases_unaligned / $total_number_of_bases_of_reads * 10000) / 100;
#print "total_number_of_bases_unaligned = $total_number_of_bases_unaligned\n";
print "% bases unaligned: $percent_bases_unaligned%\n";

$total_num_unique_aligners = $total_number_of_bases_aligned_correctly + $total_number_of_bases_aligned_incorrectly;
$accuracy_on_unique_aligners = int($total_number_of_bases_aligned_correctly / $total_num_unique_aligners * 10000) / 100;
print "% unique aligners correct: $accuracy_on_unique_aligners\n";
