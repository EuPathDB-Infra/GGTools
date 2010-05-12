#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

$|=1;

if(@ARGV < 3) {
    die "
Usage: finalcleanup.pl <rum_unique> <rum_nu> <genome seq> [options]

Options:
   -faok  : the fasta file already has sequence all on one line

This is the last script run on the RUM_Unique and RUM_NU files to clean
up things like mismatches at the ends of alignments.  This is necessary
because Bowtie extends alignments to terminal mismatches and does not
report the location of the mismatches in the output.

";

}
$faok = "false";
for($i=3; $i<@ARGV; $i++) {
    $optionrecognized = 0;
    if($ARGV[$i] eq "-faok") {
	$faok = "true";
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
    `perl modify_fa_to_have_seq_on_one_line.pl $ARGV[1] > $f`;
    open(GENOMESEQ, $f);
} else {
    open(GENOMESEQ, $ARGV[2]);
}

while($line = <GENOMESEQ>) {
    chomp($line);
    $line =~ />(.*):1-(\d+)_strand=./;
    $chr = $1;
    print STDERR "chr=$chr\n";
    $ref_seq = <GENOMESEQ>;
    chomp($ref_seq);
    open(INFILE, $ARGV[0]);
    while($line = <INFILE>) {
	chomp($line);
	@a = split(/\t/,$line);
	if($a[1] eq $chr) {
	    @b = split(/, /, $a[2]);
	    $SEQ = "";
	    for($i=0; $i<@b; $i++) {
		@c = split(/-/,$b[$i]);
		$len = $c[1] - $c[0] + 1;
		$start = $c[0] - 1;
		$SEQ = $SEQ . substr($ref_seq, $start, $len);
	    }
	    $return = &trimleft($SEQ, $a[3], $a[2]);
	    print "$a[0]\t$chr\t$return\n";
	}
    }
}
close(GENOMESEQ);

#print "-----------\nSEQ= $SEQ\n";
#print "a[2]=$a[2]\n";
#print "a[3]=$a[3]\n";
#print "x=$x\n";

sub removefirst () {
    ($n_1, $spans_1, $seq_1) = @_;
    $seq_1 =~ s/://g;
    @a_1 = split(/, /, $spans_1);
    $length_1 = 0;
    @b_1 = split(/-/,$a_1[0]);
    $length_1 = $b_1[1] - $b_1[0] + 1;
    if($length_1 <= $n_1) {
	$m_1 = $n_1 - $length_1;
	$spans2_1 = $spans_1;
	$spans2_1 =~ s/^\d+-\d+, //;
	for($j_1=0; $j_1<$length_1; $j_1++) {
	    $seq_1 =~ s/^.//;
	}
	$return = removefirst($m_1, $spans2_1, $seq_1);
	return $return;
    } else {
	for($j_1=0; $j_1<$n_1; $j_1++) {
	    $seq_1 =~ s/^.//;
	}
	$spans_1 =~ /^(\d+)-/;
	$start_1 = $1 + $n_1;
	$spans_1 =~ s/^(\d+)-/$start_1-/;
	return $spans_1 . "\t" . $seq_1;
    }
}

sub removelast () {
    ($n_1, $spans_1, $seq_1) = @_;
    $seq_1 =~ s/://g;
    @a_1 = split(/, /, $spans_1);
    @b_1 = split(/-/,$a_1[@a_1-1]);
    $length_1 = $b_1[1] - $b_1[0] + 1;
    if($length_1 <= $n_1) {
	$m_1 = $n_1 - $length_1;
	$spans2_1 = $spans_1;
	$spans2_1 =~ s/, \d+-\d+$//;
	for($j_1=0; $j_1<$length_1; $j_1++) {
	    $seq_1 =~ s/.$//;
	}
	$return = removelast($m_1, $spans2_1, $seq_1);
	return $return;
    } else {
	for($j_1=0; $j_1<$n_1; $j_1++) {
	    $seq_1 =~ s/.$//;
	}
	$spans_1 =~ /-(\d+)$/;
	$end_1 = $1 - $n_1;
	$spans_1 =~ s/-(\d+)$/-$end_1/;
	return $spans_1 . "\t" . $seq_1;
    }
}

sub trimleft () {
    ($seq1_2, $seq2_2, $spans_2) = @_;
    # seq2_2 is the one that gets modified and returned

    $seq1_2 =~ s/://g;
    $seq1_2 =~ /^(.)(.)(.)(.)/;
    $genomebase_2[0] = $1;
    $genomebase_2[1] = $2;
    $genomebase_2[2] = $3;
    $genomebase_2[3] = $4;
    $seq2_2 =~ s/://g;
    $seq2_2 =~ /^(.)(.)(.)(.)/;
    $readbase_2[0] = $1;
    $readbase_2[1] = $2;
    $readbase_2[2] = $3;
    $readbase_2[3] = $4;
    $mismatch_count_2 = 0;
    for($j_2=0; $j_2<2; $j_2++) {
	if($genomebase_2[$j_2] eq $readbase_2[$j_2]) {
	    $equal_2[$j_2] = 1;
	} else {
	    $equal_2[$j_2] = 0;
	    $mismatch_count_2++;
	}
    }
    if($mismatch_count_2 == 0) {
	return $spans_2 . "\t" . $seq2_2;
    }
    if($mismatch_count_2 == 1 && $equal_2[0] == 0) {
	&removefirst(1, $spans_2, $seq2_2) =~ /^(.*)\t(.*)/;
	$spans_new_2 = $1;
	$seq2_new_2 = $2;
	$seq1_2 =~ s/^.//;
	$return = &trimleft($seq1_2, $seq2_new_2, $spans_new_2);
	return $return;
    }
    if($equal[1] == 0 || $mismatch_count_2 == 2) {
	&removefirst(2, $spans_2, $seq2_2) =~ /^(.*)\t(.*)/;
	$spans_new_2 = $1;
	$seq2_new_2 = $2;
	$seq1_2 =~ s/^..//;
	$return = &trimleft($seq1_2, $seq2_new_2, $spans_new_2);
	return $return;
    }
}
