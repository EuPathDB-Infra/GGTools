#!/usr/bin/perl

if(@ARGV < 4 && $ARGV != 2) {
    die "
Usage: search_for_read_crossing_loc.pl <data file> <<spansfile>|<chr> <start> <end>>

Either give one span on the command line as <chr> <start> <end>, or give
the name of a file of spans.

Returns all lines of the data file whose span overlaps the span
given by the command line parameters '<chr>:<start>-<end>'.
Endpoints are considered included in all the spans.

Lines in the data file must have chr, start and end as three consecutive fields
separted by tabs. They do not have to be the first three fields as long as they
are the first that follow the pattern: '([^\\t]+)\\t(\\d+)\\t(\\d+)\\t'

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
	$line =~ /([^\t]+)\t(\d+)\t(\d+)\t/;
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
	$line =~ /([^\t]+)\t(\d+)\t(\d+)\t/;
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
