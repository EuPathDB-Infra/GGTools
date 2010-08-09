#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

$|=1;

if(@ARGV < 6) {
    die "
---------------------------------------------------------------------
*** NOTE: This script is under construction and is not working yet...
---------------------------------------------------------------------

Usage: make_RUM_junctions_file.pl <rum_unique> <rum_nu> <bowtie_unique> <bowtie_nu> <junctions outfile> <genome seq> [options]

Options:
   -faok  : the fasta file already has sequence all on one line
   -minintron n : the size of the smallest intron allowed 0<n (default = 15 bp)

This script finds the junctions in the RUM_Unique and RUM_NU files
and reports them to a junctions file that can be uploaded to the UCSC
browser.

";
}

$rumU = $ARGV[0];
$rumNU = $ARGV[1];
$bowtieU = $ARGV[2];
$bowtiemNU = $ARGV[3];
$outfile = $ARGV[4];

open(INFILE, $bowtieU);
while($line = <INFILE>) {
    $line =~ /seq.(\d+.)/;
    $bowtie{$1} = 0;
}
close(INFILE);
open(INFILE, $bowtieNU);
while($line = <INFILE>) {
    $line =~ /seq.(\d+.)/;
    $bowtie{$1} = 0;
}
close(INFILE);

$faok = "false";
$minintron = 15;

for($i=6; $i<@ARGV; $i++) {
    $optionrecognized = 0;
    if($ARGV[$i] eq "-faok") {
	$faok = "true";
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-minintron") {
	$minintron = $ARGV[$i+1];
	if(!($minintron =~ /^\d+$/)) {
	    die "\nError: -minintron must be an integer greater than zero, you gave '$minintron'.\n\n";
	} else {
	    die "\nError: -minintron must be an integer greater than zero, you gave '$minintron'.\n\n";
	}
	$i++;
	$optionrecognized = 1;
    }
    if($optionrecognized == 0) {
	die "\nERROR: option '$ARGV[$i]' not recognized\n";
    }
}

if($faok eq "false") {
    print STDERR "Modifying genome fa file\n";
    $r = int(rand(1000));
    $f = "temp_" . $r . ".fa";
    `perl modify_fa_to_have_seq_on_one_line.pl $ARGV[5] > $f`;
    open(GENOMESEQ, $f);
} else {
    open(GENOMESEQ, $ARGV[5]);
}

$FLAG = 0;
while($FLAG == 0) {
    undef %CHR2SEQ;
    $sizeflag = 0;
    $totalsize = 0;
    while($sizeflag == 0) {
	$line = <GENOMESEQ>;
	if($line eq '') {
	    $FLAG = 1;
	    $sizeflag = 1;
	} else {
	    chomp($line);
	    $line =~ />(.*):1-(\d+)_strand=./;
	    $chr = $1;
	    $ref_seq = <GENOMESEQ>;
	    chomp($ref_seq);
	    $CHR2SEQ{$chr} = $ref_seq;
	    $totalsize = $totalsize + length($ref_seq);
	    if($totalsize > 1000000000) {  # don't store more than 1 gb of sequence in memory at once...
		$sizeflag = 1;
	    }
	}
    }
    &getjunctions();

}
close(GENOMESEQ);

sub getjunctions () {
    open(INFILE, $rumU);
    while($line = <INFILE>) {
	$flag = 0;
	if(!($line =~ /, /)) {
	    next;
	}
	chomp($line);
	@a = split(/\t/,$line);
	if(!(defined $CHR2SEQ{$a[1]})) {
	    next;
	}
	while($seq =~ /\+/) {
	    $seq =~ s/\+[^+]\+//;
	}
	$strand = $a[3];
	$chr = $a[1];
	@SPANS = split(/, /,$a[2]);
	@SEQ = split(/:/, $seq);
	for($i=0; $i<@SPANS-1; $i++) {
	    @c1 = split(/-/,$SPANS[$i]);
	    @c2 = split(/-/,$SPANS[$i+1]);
	    $elen1 = $c1[1] - $c1[0] + 1;
	    $elen2 = $c2[1] - $c2[0] + 1;
	    $ilen = $c2[0] - $c1[1] + 1;
	    $istart = $c1[1]+1;
	    $iend = $c2[0]-1;
	    $intron = $chr . ":" . $istart . "-" . $iend;
	    if(defined $bowtie{$a[0]}) {
		$hqualU{$intron}++;
	    } elsif($ilen >= $minintron) {
		$SEQ[$i] =~ /.$/;
		$leftexon_lastbase = $1;
		$SEQ[$i+1] =~ /^./;
		$rightexon_firstbase = $1;
		$intron_firstbase = substr($CHR2SEQ{$chr}, $istart-1, 1);
		$intron_lastbase = substr($CHR2SEQ{$chr}, $iend-1, 1);
		if($leftexon_lastbase eq $intron_lastbase) {
		    $istart_alt = $istart-1;
		    $iend_alt = $iend-1;
		    $altintron = $chr . ":" . $istart_alt . "-" . $iend_alt;
		    $amb{$intron}=1;  # amb for ambiguous
		    $amb{$altintron}=1;
		}
		if($rightexon_firstbase eq $intron_firtbase) {
		    $istart_alt = $istart+1;
		    $iend_alt = $iend+1;
		    $altintron = $chr . ":" . $istart_alt . "-" . $iend_alt;
		    $amb{$intron}=1;  # amb for ambiguous
		    $amb{$altintron}=1;
		}
		if($elen1 <= 6 || $elen2 <=6) {
		    $weakevidence{$junction}++;
		} else {
		    $strongevidence{$junction}++;
		}
	    }
	}
    }
    close(INFILE);

    open(INFILE, $rumNU);
    while($line = <INFILE>) {
	$flag = 0;
	if(!($line =~ /, /)) {
	    next;
	}
	chomp($line);
	@a = split(/\t/,$line);
	if(!(defined $CHR2SEQ{$a[1]})) {
	    next;
	}
	while($seq =~ /\+/) {
	    $seq =~ s/\+[^+]\+//;
	}
	$strand = $a[3];
	$chr = $a[1];
	@SPANS = split(/, /,$a[2]);
	@SEQ = split(/:/, $seq);
	for($i=0; $i<@SPANS-1; $i++) {
	    @c1 = split(/-/,$SPANS[$i]);
	    @c2 = split(/-/,$SPANS[$i+1]);
	    $elen1 = $c1[1] - $c1[0] + 1;
	    $elen2 = $c2[1] - $c2[0] + 1;
	    $ilen = $c2[0] - $c1[1] + 1;
	    $istart = $c1[1]+1;
	    $iend = $c2[0]-1;
	    $junctions = $chr . ":" . $istart . "-" . $iend;
	    if(defined $bowtie{$a[0]}) {
		$hqualNU{$junction}++;
	    } elsif($ilen >= $minintron) {
		if($elen1 <= 6 || $elen2 <=6) {
		    $lqualNU{$junction}++;
		} else {
		    $hqualNU{$junction}++;
		}
	    }
	}
    }
    close(INFILE);

    open(OUTFILE, ">>$outfile");
    close(OUTFILE);
}

#	    substr($CHR2SEQ{$chr}, XXX, 1);
