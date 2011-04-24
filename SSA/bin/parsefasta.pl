#!/usr/bin/perl
use strict;

if(@ARGV<3) {
    die "
Usage: parsefasta.pl <infile> <num chunks> <reads out>

";
}

my $infile = $ARGV[0];
my $numchunks = $ARGV[1];
my $reads_out = $ARGV[2];

my $paired = "true";
my $infile1;
my $infile2;

if($infile =~ /,,,/) {
    $infile =~ /^(.*),,,(.*)$/;
    $infile1 = $1;
    $infile2 = $2
} else {
    $infile1 = $infile;
    $paired = "false";
}
open(INFILE1, $infile1) or die "\nERROR: in script parsefasta.pl: cannot open '$infile1' for reading\n";
if($paired eq "true") {
    open(INFILE2, $infile2) or die "\nERROR: in script parsefasta.pl: cannot open '$infile2' for reading\n";
}

my $filesize = -s $infile1;

# put something here for the case the file is less than 10,000 lines (or 2,500 entries)

my $FL = `head -10000 $infile1 | wc -l`;
chomp($FL);
$FL =~ s/[^\d]//gs;

my $s1 = `head -$FL $infile1`;
my $s2 = `tail -$FL $infile1`;
my $totalsize = length($s1) + length($s2);
my $recordsize = $totalsize / $FL;
my $numrecords = int($filesize / $recordsize);
my $numrecords_per_chunk = int($numrecords / $numchunks);

my $seq_counter = 0;
my $endflag = 0;
open(ROUTALL, ">$reads_out");
if($paired eq "false") {
    for(my $chunk=1; $chunk<=$numchunks; $chunk++) {
	my $reads_file = $reads_out . ".$chunk";
	if($endflag == 1) {
	    $chunk = $numchunks;
	    next;
	}
	open(ROUT, ">$reads_file");
	if($chunk == $numchunks) {
	    # just to make sure we get everything in the last chunk
	    $numrecords_per_chunk = $numrecords_per_chunk * 100; 
	}
	for(my $i=0; $i<$numrecords_per_chunk; $i++) {
	    $seq_counter++;
	    my $line = <INFILE1>;
	    chomp($line);
	    my $line_hold = $line;
	    $line = <INFILE1>;
	    chomp($line);
	    if($line eq '' && $line_hold ne '') {
		print STDERR "ERROR: in script parsefasta.pl: something is wrong, the file seems to end with an incomplete record...\n";
		exit(0);
	    }
	    if($line eq '') {
		$i = $numrecords_per_chunk;
		$endflag = 1;
		next;
	    }
	    print ROUT ">seq.$seq_counter";
	    print ROUTALL ">seq.$seq_counter";
	    print ROUT "a\n";
	    print ROUTALL "a\n";
	    $line =~ s/\./N/g;
	    $line = uc $line;
	    print ROUT "$line\n";
	    print ROUTALL "$line\n";
	}
	close(ROUT);
    }
}

if($paired eq "true") {
    for(my $chunk=1; $chunk<=$numchunks; $chunk++) {
	my $reads_file = $reads_out . ".$chunk";
	if($endflag == 1) {
	    $chunk = $numchunks;
	    next;
	}
	open(ROUT, ">$reads_file");
	if($chunk == $numchunks) {
	    # just to make sure we get everything in the last chunk
	    $numrecords_per_chunk = $numrecords_per_chunk * 100; 
	}
	for(my $i=0; $i<$numrecords_per_chunk; $i++) {
	    $seq_counter++;
	    my $line = <INFILE1>;
	    chomp($line);
	    my $line_hold = $line;
	    $line = <INFILE1>;
	    chomp($line);
	    if($line eq '' && $line_hold ne '') {
		print STDERR "ERROR: in script parsefasta.pl: something is wrong, the forward file seems to end with an incomplete record...\n";
		exit(0);
	    }
	    if($line eq '') {
		$i = $numrecords_per_chunk;
		$endflag = 1;
		next;
	    }
	    print ROUT ">seq.$seq_counter";
	    print ROUTALL ">seq.$seq_counter";
	    print ROUT "a\n";
	    print ROUTALL "a\n";
	    $line =~ s/\./N/g;
	    $line = uc $line;
	    print ROUT "$line\n";
	    print ROUTALL "$line\n";

	    $line = <INFILE2>;
	    chomp($line);
	    $line_hold = $line;
	    $line = <INFILE2>;
	    chomp($line);
	    if($line eq '' && $line_hold ne '') {
		$i = $numrecords_per_chunk;
		print STDERR "ERROR: in script parsefasta.pl: something is wrong, the reverse file seems to end with an incomplete record...\n";
		exit(0);
	    }
	    print ROUT ">seq.$seq_counter";
	    print ROUTALL ">seq.$seq_counter";
	    print ROUT "b\n";
	    print ROUTALL "b\n";
	    $line =~ s/\./N/g;
	    print ROUT "$line\n";
	    print ROUTALL "$line\n";
	}
	close(ROUT);
    }
}

close(INFILE1);
if($paired eq "true") {
    close(INFILE2);
}
close(ROUTALL);
