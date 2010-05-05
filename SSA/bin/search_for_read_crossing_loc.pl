#!/usr/bin/perl

if(@ARGV < 4 && $ARGV != 2) {
    die "
Usage: search_bed_for_read_crossing_loc.pl <bed file> <<spansfile>|<chr> <loc1> <loc2>>

Either give one span on the command line as <chr> <loc1> <loc2>, or give
the name of a file of spans.

Returns all lines of the bed file with span that overlaps the span
given by the command line parameters '<chr>:<loc1>-<loc2>'

File must have chr, start and stop as three consecutive fields.
They do not have to be the first three fields as long as they are 
the first that follow the pattern: '\\t([^\\t]+)\\t(\\d+)\\t(\\d+)\\t'

";
}

if(!(-e $ARGV[1])) {
    $chr = $ARGV[1];
    $loc1 = $ARGV[2];
    $loc2 = $ARGV[3];
    open(INFILE, $ARGV[0]);
    while($line = <INFILE>) {
	chomp($line);
	@a = split(/\t/,$line);
	$line =~ /\t([^\t]+)\t(\d+)\t(\d+)\t/;
	$CHR = $1;
	$START = $2;
	$END = $3;
	if($chr eq $CHR && $START <= $loc2 && $loc1 <= $END) {
	    print "$line\n";
	}
    }
    close(INFILE);
} else {
    open(SPANSFILE, $ARGV[1]);
    $cnt=0;
    while($line = <SPANSFILE>) {
	chomp($line);
	@a = split(/\t/,$line);
	$chr[$cnt] = $a[0];
	$loc1[$cnt] = $a[1];
	$loc2[$cnt] = $a[2];
	$cnt++;
    }
    close(SPANSFILE);
    open(INFILE, $ARGV[0]);
    while($line = <INFILE>) {
	chomp($line);
	@a = split(/\t/,$line);
	$line =~ /\t([^\t]+)\t(\d+)\t(\d+)\t/;
	$CHR = $1;
	$START = $2;
	$END = $3;
	for($i=0; $i<@chr; $i++) {
	    if($chr[$i] eq $CHR && $START <= $loc2[$i] && $loc1[$i] <= $END) {
		print "$line\n";
	    }
	}
    }
    close(INFILE);
}
