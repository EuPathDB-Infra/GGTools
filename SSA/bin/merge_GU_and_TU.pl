#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

$|=1;

if(@ARGV < 7) {
    die "
Usage: merge_GU_and_TU.pl <GU infile> <TU infile> <GNU infile> <TNU infile> <BowtieUnique outfile> <CNU outfile> <type>

Where:   <GU infile> is the file of unique mappers that is output from the
                     script make_GU_and_GNU.pl

         <TU infile> is the file of unique mappers that is output from the
                     script make_TU_and_TNU.pl

         <GNU infile> is the file of non-unique mappers that is output from
                      the script make_GU_and_GNU.pl

         <TNU infile> is the file of non-unique mappers that is output from
                      the script make_TU_and_TNU.pl

         <BowtieUnique outfile> is the name of the file of unique mappers to be output

         <CNU outfile> is the name of the file of non-unique mappers to be output

         <type> is 'single' for single-end reads, or 'paired' for paired-end reads

  Options:
         -maxpairdist N : N is an integer greater than zero representing
                          the furthest apart the forward and reverse reads
                          can be.  They could be separated by an exon/exon
                          junction so this number can be as large as the largest
                          intron.  Default value = 500,000

";
}

$infile1 = $ARGV[2];
$infile2 = $ARGV[3];
$infile3 = $ARGV[0];
$infile4 = $ARGV[1];
$outfile1 = $ARGV[4];
$outfile2 = $ARGV[5];
$type = $ARGV[6];
$typerecognized = 1;
if($type eq "single") {
    $paired_end = "false";
    $typerecognized = 0;
}
if($type eq "paired") {
    $paired_end = "true";
    $typerecognized = 0;
}
if($typerecognized == 1) {
    die "\nERROR: type '$type' not recognized.  Must be 'single' or 'paired'.\n";
}

$max_distance_between_paired_reads = 500000;
for($i=7; $i<@ARGV; $i++) {
    $optionrecognized = 0;
    if($ARGV[$i] eq "-maxpairdist") {
	$i++;
	$max_distance_between_paired_reads = $ARGV[$i];
	$optionrecognized = 1;
    }

    if($optionrecognized == 0) {
	die "\nERROR: option '$ARGV[$i-1] $ARGV[$i]' not recognized\n";
    }
}

open(INFILE, $infile4) or die "\nERROR: Cannot open file '$infile4' for reading\n";
$cnt = 0;
while($line = <INFILE>) {
    if($line =~ /seq.\d+a/ || $line =~ /seq.\d+b/) {
	chomp($line);
	@a = split(/\t/,$line);
	$span = $a[2];
	if(!($span =~ /,/)) {
	    $cnt++;
	    @b = split(/-/,$span);
	    $length = $b[1] - $b[0] + 1;
	    if($length > $readlength) {
		$readlength = $length;
		$cnt = 0;
	    }
	    if($cnt > 50000) { # it checked 50,000 lines without finding anything larger than the last time
		# readlength was changed, so it's most certainly found the max.
		# Went through this to avoid the user having to input the readlength.
		last;
	    }
	}
    }
}
$cnt = 0;
open(INFILE, $infile3) or die "\nERROR: Cannot open file '$infile3' for reading\n";
while($line = <INFILE>) {
    if($line =~ /seq.\d+a/ || $line =~ /seq.\d+b/) {
	chomp($line);
	@a = split(/\t/,$line);
	$span = $a[2];
	if(!($span =~ /,/)) {
	    $cnt++;
	    @b = split(/-/,$span);
	    $length = $b[1] - $b[0] + 1;
	    if($length > $readlength) {
		$readlength = $length;
		$cnt = 0;
	    }
	    if($cnt > 50000) { # it checked 50,000 lines without finding anything larger than the last time
		# readlength was changed, so it's most certainly found the max.
		# Went through this to avoid the user having to input the readlength.
		last;
	    }
	}
    }
}
close(INFILE);
$cnt = 0;
open(INFILE, $infile1) or die "\nERROR: Cannot open file '$infile1' for reading\n";
while($line = <INFILE>) {
    if($line =~ /seq.\d+a/ || $line =~ /seq.\d+b/) {
	chomp($line);
	@a = split(/\t/,$line);
	$span = $a[2];
	if(!($span =~ /,/)) {
	    $cnt++;
	    @b = split(/-/,$span);
	    $length = $b[1] - $b[0] + 1;
	    if($length > $readlength) {
		$readlength = $length;
		$cnt = 0;
	    }
	    if($cnt > 50000) { # it checked 50,000 lines without finding anything larger than the last time
		# readlength was changed, so it's most certainly found the max.
		# Went through this to avoid the user having to input the readlength.
		last;
	    }
	}
    }
}
close(INFILE);
$cnt = 0;
open(INFILE, $infile2) or die "\nERROR: Cannot open file '$infile3' for reading\n";
while($line = <INFILE>) {
    if($line =~ /seq.\d+a/ || $line =~ /seq.\d+b/) {
	chomp($line);
	@a = split(/\t/,$line);
	$span = $a[2];
	if(!($span =~ /,/)) {
	    $cnt++;
	    @b = split(/-/,$span);
	    $length = $b[1] - $b[0] + 1;
	    if($length > $readlength) {
		$readlength = $length;
		$cnt = 0;
	    }
	    if($cnt > 50000) { # it checked 50,000 lines without finding anything larger than the last time
		# readlength was changed, so it's most certainly found the max.
		# Went through this to avoid the user having to input the readlength.
		last;
	    }
	}
    }
}
close(INFILE);

print "readlength = $readlength\n";

if($readlength < 80) {
    $min_overlap = 35;
} else {
    $min_overlap = 45;
}
if($min_overlap >= .8 * $readlength) {
    $min_overlap = int(.6 * $readlength);
}

open(INFILE, $infile1) or die "\nERROR: Cannot open file '$infile1' for reading\n";
while($line = <INFILE>) {
    $line =~ /^seq.(\d+)/;
    $ambiguous_mappers{$1}++;
}
close(INFILE);
open(INFILE, $infile2) or die "\nERROR: Cannot open file '$infile2' for reading\n";
while($line = <INFILE>) {
    $line =~ /^seq.(\d+)/;
    $ambiguous_mappers{$1}++;
}
close(INFILE);
open(INFILE1, $infile3) or die "\nERROR: Cannot open file '$infile3' for reading\n";
open(INFILE2, $infile4) or die "\nERROR: Cannot open file '$infile4' for reading\n";
open(OUTFILE1, ">$outfile1") or die "\nERROR: Cannot open file '$outfile1' for writing\n";
open(OUTFILE2, ">$outfile2") or die "\nERROR: Cannot open file '$outfile2' for writing\n";

$num_lines_at_once = 10000;
$linecount = 0;
$FLAG = 1;
$line_prev = <INFILE2>;
chomp($line_prev);
while($FLAG == 1) {
    undef %hash1;
    undef %hash2;
    undef %allids;
    $linecount = 0;
    until($linecount == $num_lines_at_once) {
	$line=<INFILE1>;
	if(!($line =~ /\S/)) {
	    $FLAG = 0;
	    $linecount = $num_lines_at_once;
	}
	else {
	    chomp($line);
	    @a = split(/\t/,$line);
	    $a[0] =~ /seq.(\d+)/;
	    $id = $1;
	    $last_id = $id;
	    $allids{$id}++;
	    if($a[0] =~ /a$/ || $a[0] =~ /b$/) {
		$hash1{$id}[0]++;
		$hash1{$id}[$hash1{$id}[0]]=$line;
	    }
	    else {
		$hash1{$id}[0]=-1;
		$hash1{$id}[1]=$line;
	    }
	    if($paired_end eq "true") {
		# this makes sure we have read in both a and b reads, this approach might cause a problem
		# if no, or very few, b reads mapped at all.
		if( (($linecount == ($num_lines_at_once - 1)) && !($a[0] =~ /a$/)) || ($linecount < ($num_lines_at_once - 1)) ) {
		    $linecount++;
		}
	    }
	    else {
		if( ($linecount == ($num_lines_at_once - 1)) || ($linecount < ($num_lines_at_once - 1)) ) {
		    $linecount++;
		}
	    }
	}
    }
    $line = $line_prev;
    @a = split(/\t/,$line);
    $a[0] =~ /seq.(\d+)/;
    $prev_id = $id;
    $id = $1;
    if($prev_id eq $id) {
	$FLAG2 = 0;
    }
    $FLAG2 = 1;
    until($id > $last_id || $FLAG2 == 0) {
	$allids{$id}++;
	if($a[0] =~ /a$/ || $a[0] =~ /b$/) {
	    $hash2{$id}[0]++;
	    $hash2{$id}[$hash2{$id}[0]]=$line;
	}
	else {
	    $hash2{$id}[0]=-1;
	    $hash2{$id}[1]=$line;
	}
	$line=<INFILE2>;
	chomp($line);
	if(!($line =~ /\S/)) {
	    $FLAG2 = 0;
	}
	else {
	    @a = split(/\t/,$line);
	    $a[0] =~ /seq.(\d+)/;
	    $id = $1;
	}
    }
    if($FLAG2 == 1) {
	$line_prev = $line;
    }
    foreach $id (sort {$a <=> $b} keys %allids) {
	if($ambiguous_mappers{$id}+0 > 0) {
	    next;
	}
	$hash1{$id}[0] = $hash1{$id}[0] + 0;
	$hash2{$id}[0] = $hash2{$id}[0] + 0;
	# MUST DO 15 CASES IN TOTAL:
	# THREE CASES:
	if($hash1{$id}[0] == 0) {
	    if($hash2{$id}[0] == -1) {
		print OUTFILE1 "$hash2{$id}[1]\n";
	    }
	    else {
		for($i=0; $i<$hash2{$id}[0]; $i++) {
		    print OUTFILE1 "$hash2{$id}[$i+1]\n";
		}
	    }
	}
	# THREE CASES
	if($hash2{$id}[0] == 0) {
	    if($hash1{$id}[0] == -1) {
		print OUTFILE1 "$hash1{$id}[1]\n";
	    }
	    else {
		for($i=0; $i<$hash1{$id}[0]; $i++) {
		    print OUTFILE1 "$hash1{$id}[$i+1]\n";
		}
	    }
	}
	# ONE CASE
	if($hash1{$id}[0] == -1 && $hash2{$id}[0] == -1) {
	    undef @spans;
	    @a1 = split(/\t/,$hash1{$id}[1]);
	    @a2 = split(/\t/,$hash2{$id}[1]);
	    $spans[0] = $a1[2];
	    $spans[1] = $a2[2];
	    $str = intersect(\@spans, $a1[3]);
	    $str =~ /^(\d+)/;
	    $length_overlap = $1;
	    if(($length_overlap > $min_overlap) && ($a1[1] eq $a2[1])) {
		print OUTFILE1 "$hash2{$id}[1]\n";
	    }
	    else {
		print OUTFILE2 "$hash1{$id}[1]\n";
		print OUTFILE2 "$hash2{$id}[1]\n";
	    }
	}
	# ONE CASE
	if($hash1{$id}[0] == 1 && $hash2{$id}[0] == 1) {
	    # If single-end then this is the only case where $hash1{$id}[0] > 0 and $hash2{$id}[0] > 0
	    if((($hash1{$id}[1] =~ /seq.\d+a/) && ($hash2{$id}[1] =~ /seq.\d+a/)) || (($hash1{$id}[1] =~ /seq.\d+b/) && ($hash2{$id}[1] =~ /seq.\d+b/))) {
		undef @spans;
		@a1 = split(/\t/,$hash1{$id}[1]);
		@a2 = split(/\t/,$hash2{$id}[1]);
		$spans[0] = $a1[2];
		$spans[1] = $a2[2];
		$str = intersect(\@spans, $a1[3]);
		$str =~ /^(\d+)/;
		$length_overlap = $1;
		if(($length_overlap > $min_overlap) && ($a1[1] eq $a2[1])) {
                    # preference TU
		    print OUTFILE1 "$hash2{$id}[1]\n";
		}
		else {
		    if($paired_end eq "false") {
			print OUTFILE2 "$hash1{$id}[1]\n";			
			print OUTFILE2 "$hash2{$id}[1]\n";			
		    }
		}
	    }
	    if((($hash1{$id}[1] =~ /seq.\d+a/) && ($hash2{$id}[1] =~ /seq.\d+b/)) || (($hash1{$id}[1] =~ /seq.\d+b/) && ($hash2{$id}[1] =~ /seq.\d+a/))) {
		@a = split(/\t/,$hash1{$id}[1]);
		$aspans = $a[2];
		$a[2] =~ /^(\d+)[^\d]/;
		$astart = $1;
		$a[2] =~ /[^\d](\d+)$/;
		$aend = $1;
		$chra = $a[1];
		$aseq = $a[3];
		$seqnum = $a[0];
		$seqnum =~ s/a$//;
		$seqnum =~ s/b$//;
		@a = split(/\t/,$hash2{$id}[1]);
		$bspans = $a[2];
		$a[2] =~ /^(\d+)[^\d]/;
		$bstart = $1;
		$a[2] =~ /[^\d](\d+)$/;
		$bend = $1;
		$chrb = $a[1];
		$bseq = $a[3];
		if(($chra eq $chrb) && ($aend < $bstart-1) && ($bstart - $aend < $max_distance_between_paired_reads)) {
		    if($hash1{$id}[1] =~ /a\t/) {
			print OUTFILE1 "$hash1{$id}[1]\n$hash2{$id}[1]\n";
		    }
		    else {
			print OUTFILE1 "$hash2{$id}[1]\n$hash1{$id}[1]\n";
		    }
		}
		if(($chra eq $chrb) && ($bend < $astart-1) && ($astart - $bend < $max_distance_between_paired_reads)) {
		    if($hash1{$id}[1] =~ /a\t/) {
			print OUTFILE1 "$hash1{$id}[1]\n$hash2{$id}[1]\n";
		    }
		    else {
			print OUTFILE1 "$hash2{$id}[1]\n$hash1{$id}[1]\n";
		    }
		}
		$Eflag =0;
		if(($chra eq $chrb) && ($aend >= $bstart-1) && ($astart <= $bstart)) {
		    $overlap = $aend - $bstart + 1;
		    $seq_merged = $aseq;
		    $bseq =~ s/://g;
		    @s = split(//,$bseq);
		    for($i=$overlap; $i<@s; $i++) {
			$seq_merged = $seq_merged . $s[$i];
		    }
		    $spans_merged = $bspans;
		    $spans_merged =~ s/^$bstart/$astart/;
		    $seq_j = addJunctionsToSeq($seq_merged, $spans_merged);
		    print OUTFILE1 "$seqnum\t$chra\t$spans_merged\t$seq_j\n";
		    $Eflag =1;
		}
		if(($chra eq $chrb) && ($bend >= $astart-1) && ($bstart <= $astart) && ($Eflag == 0)) {
		    $overlap = $bend - $astart + 1;
		    $bseq =~ s/://g;
		    $seq_merged = $bseq;
		    @s = split(//,$aseq);
		    for($i=$overlap; $i<@s; $i++) {
			$seq_merged = $seq_merged . $s[$i];
		    }
		    $spans_merged = $bspans;
		    $spans_merged =~ s/$bend$/$aend/;
		    $seq_j = addJunctionsToSeq($seq_merged, $spans_merged);
		    print OUTFILE1 "$seqnum\t$chra\t$spans_merged\t$seq_j\n";
		}
	    }
	}
	# ONE CASE
	if($hash1{$id}[0] == 2 && $hash2{$id}[0] == 2) {
	    undef @spansa;
	    undef @spansb;
	    @a = split(/\t/,$hash1{$id}[1]);
	    $chr1 = $a[1];
	    $spansa[0] = $a[2];
	    $seqa = $a[3];
	    @a = split(/\t/,$hash1{$id}[2]);
	    $spansb[0] = $a[2];
	    $seqb = $a[3];
	    @a = split(/\t/,$hash2{$id}[1]);
	    $chr2 = $a[1];
	    $spansa[1] = $a[2];
	    @a = split(/\t/,$hash2{$id}[2]);
	    $spansb[1] = $a[2];
	    $str = intersect(\@spansa, $seqa);
	    $str =~ /^(\d+)/;
	    $length_overlap1 = $1;
	    $str = intersect(\@spansb, $seqb);
	    $str =~ /^(\d+)/;
	    $length_overlap2 = $1;
	    if(($length_overlap1 > $min_overlap) && ($length_overlap2 > $min_overlap) && ($chr1 eq $chr2)) {
		print OUTFILE1 "$hash2{$id}[1]\n";
		print OUTFILE1 "$hash2{$id}[2]\n";
	    }
	    else {
		print OUTFILE2 "$hash1{$id}[1]\n";
		print OUTFILE2 "$hash1{$id}[2]\n";
		print OUTFILE2 "$hash2{$id}[1]\n";
		print OUTFILE2 "$hash2{$id}[2]\n";
	    }
	}	
	# NINE CASES DONE
	# ONE CASE
	if($hash1{$id}[0] == -1 && $hash2{$id}[0] == 2) {
	    print OUTFILE2 "$hash1{$id}[1]\n";
	    print OUTFILE2 "$hash2{$id}[1]\n";
	    print OUTFILE2 "$hash2{$id}[2]\n";
	}
	# ONE CASE
	if($hash1{$id}[0] == 2 && $hash2{$id}[0] == -1) {
	    undef @spans;
	    @a = split(/\t/,$hash1{$id}[1]);
	    $chr1 = $a[1];
	    $spans[0] = $a[2];
	    $seq = $a[3];
	    @a = split(/\t/,$hash2{$id}[1]);
	    $chr2 = $a[1];
	    $spans[1] = $a[2];
	    if($chr1 eq $chr2) {
		$str = intersect(\@spans, $seq);
		$str =~ /^(\d+)/;
		$overlap1 = $1;
		@a = split(/\t/,$hash1{$id}[2]);
		$spans[0] = $a[2];
		$str = intersect(\@spans, $seq);
		$str =~ /^(\d+)/;
		$overlap2 = $1;
	    }
	    if($overlap1 >= $min_overlap && $overlap2 >= $min_overlap) {
		print OUTFILE1 "$hash2{$id}[1]\n";
	    }
	    else {
		print OUTFILE2 "$hash1{$id}[1]\n";
		print OUTFILE2 "$hash1{$id}[2]\n";
		print OUTFILE2 "$hash2{$id}[1]\n";
	    }
	}
	# ELEVEN CASES DONE
	if($hash1{$id}[0] == -1 && $hash2{$id}[0] == 1) {
	    print OUTFILE1 "$hash1{$id}[1]\n";
	}
	if($hash1{$id}[0] == 1 && $hash2{$id}[0] == -1) {
	    print OUTFILE1 "$hash2{$id}[1]\n";
	}
	if($hash1{$id}[0] == 1 && $hash2{$id}[0] == 2) {
	    print OUTFILE1 "$hash2{$id}[1]\n";
	    print OUTFILE1 "$hash2{$id}[2]\n";
	}	
	if($hash1{$id}[0] == 2 && $hash2{$id}[0] == 1) {
	    print OUTFILE1 "$hash1{$id}[1]\n";
	    print OUTFILE1 "$hash1{$id}[2]\n";
	}	
	# ALL FIFTEEN CASES DONE
    }
}


sub intersect () {
    ($spans_ref, $seq) = @_;
    @spans = @{$spans_ref};
    $num = @spans;
    undef %chash;
    for($s=0; $s<$num; $s++) {
	@a = split(/, /,$spans[$s]);
	for($i=0;$i<@a;$i++) {
	    @b = split(/-/,$a[$i]);
	    for($j=$b[0];$j<=$b[1];$j++) {
		$chash{$j}++;
	    }
	}
    }
    $spanlength = 0;
    $flag = 0;
    $maxspanlength = 0;
    $maxspan_start = 0;
    $maxspan_end = 0;
    $prevkey = 0;
    for $key (sort {$a <=> $b} keys %chash) {
	if($chash{$key} == $num) {
	    if($flag == 0) {
		$flag = 1;
		$span_start = $key;
	    }
	    $spanlength++;
	}
	else {
	    if($flag == 1) {
		$flag = 0;
		if($spanlength > $maxspanlength) {
		    $maxspanlength = $spanlength;
		    $maxspan_start = $span_start;
		    $maxspan_end = $prevkey;
		}
		$spanlength = 0;
	    }
	}
	$prevkey = $key;
    }
    if($flag == 1) {
	if($spanlength > $maxspanlength) {
	    $maxspanlength = $spanlength;
	    $maxspan_start = $span_start;
	    $maxspan_end = $prevkey;
	}
    }
    if($maxspanlength > 0) {
	@a = split(/, /,$spans[0]);
	@b = split(/-/,$a[0]);
	$i=0;
	until($b[1] >= $maxspan_start) {
	    $i++;
	    @b = split(/-/,$a[$i]);
	}
	$prefix_size = $maxspan_start - $b[0];  # the size of the part removed from spans[0]
	for($j=0; $j<$i; $j++) {
	    @b = split(/-/,$a[$j]);
	    $prefix_size = $prefix_size + $b[1] - $b[0] + 1;
	}
	@s = split(//,$seq);
	$newseq = "";
	for($i=$prefix_size; $i<$prefix_size + $maxspanlength; $i++) {
	    $newseq = $newseq . $s[$i];
	}
	$flag = 0;
	$i=0;
	@b = split(/-/,$a[0]);
	until($b[1] >= $maxspan_start) {
	    $i++;
	    @b = split(/-/,$a[$i]);
	}
	$newspans = $maxspan_start;
	until($b[1] >= $maxspan_end) {
	    $newspans = $newspans . "-$b[1]";
	    $i++;
	    @b = split(/-/,$a[$i]);
	    $newspans = $newspans . ", $b[0]";
	}
	$newspans = $newspans . "-$maxspan_end";
	$off = "";
	for($i=0; $i<$prefix_size; $i++) {
	    $off = $off . " ";
	}
	return "$maxspanlength\t$newspans\t$newseq";
    }
    else {
	return "0";
    }
}

sub addJunctionsToSeq () {
    ($seq, $spans) = @_;
    @s = split(//,$seq);
    @b = split(/, /,$spans);
    $seq_out = "";
    $place = 0;
    for($j=0; $j<@b; $j++) {
	@c = split(/-/,$b[$j]);
	$len = $c[1] - $c[0] + 1;
	if($seq_out =~ /\S/) {
	    $seq_out = $seq_out . ":";
	}
	for($k=0; $k<$len; $k++) {
	    $seq_out = $seq_out . $s[$place];
	    $place++;
	}
    }
    return $seq_out;
}
