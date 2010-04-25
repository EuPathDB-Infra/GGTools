#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

$|=1;

if(@ARGV < 7) {
    die "
Usage: merge_BowtieUnique_and_BlatUnique.pl <Bowtie Unique infile> <BLAT Unique infile> <Bowtie NU infile> <BLAT NU infile> <RUM Unique outfile> <RUM NU outfile> <type>

Where:   

         <Bowtie Unique infile> is the file of unique mappers that is output from
                                the script merge_GU_and_TU.pl

         <BLAT Unique infile> is the file of unique mappers that is output from the
                              script parse_blat_out.pl

         <Bowtie NU infile> is the file of non-unique mappers that is output from
                            the script merge_GNU_and_TNU_and_CNU.pl.

         <Blat NU infile> is the file of non-unique mappers that is output from
                          the script parse_blat_out.pl

         <RUM_Unique outfile> is the name of the file of unique mappers to be output

         <RUM_NU outfile> is the name of the file of non-unique mappers to be output

         <type> is 'single' for single-end reads, or 'paired' for paired-end reads

  Options:
         -maxpairdist N : N is an integer greater than zero representing
                          the furthest apart the forward and reverse reads
                          can be.  They could be separated by an exon/exon
                          junction so this number can be as large as the largest
                          intron.  Default value = 500,000

";
}

# get readlength from bowtie unique/nu, if both empty then get max in blat unique/nu

open(INFILE, $ARGV[0]) or die "\nError: unable to open file '$ARGV[0]' for reading\n\n";
$readlength = 0;
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
close(INFILE);
open(INFILE, $ARGV[1]) or die "\nError: unable to open file '$ARGV[1]' for reading\n\n";
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
	    if($cnt > 50000) {
		last;
	    }
	}
    }
}
close(INFILE);
open(INFILE, $ARGV[2]) or die "\nError: unable to open file '$ARGV[2]' for reading\n\n";
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
	    if($cnt > 50000) {
		last;
	    }
	}
    }
}
close(INFILE);
open(INFILE, $ARGV[3]) or die "\nError: unable to open file '$ARGV[3]' for reading\n\n";
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
	    if($cnt > 50000) {
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

$f0 = $ARGV[3];
open(INFILE, $f0) or die "\nError: unable to open file '$f0' for reading\n\n";

while($line = <INFILE>) {
    $line =~ /^seq.(\d+)/;
    $id = $1;
    if($line =~ /seq.\d+a/) {
	$blat_ambiguous_mappers_a{$id}++;
    }
    if($line =~ /seq.\d+b/) {
	$blat_ambiguous_mappers_b{$id}++;
    }
    if(!($line =~ /(seq.\d+a|seq.\d+b)/)) {
	$blat_ambiguous_mappers_a{$id}++;
	$blat_ambiguous_mappers_b{$id}++;
    }
}
close(INFILE);

open(OUTFILE2, ">>$f0") or die "\nError: unable to open file '$f0' for writing\n\n";
                          # The only things we're going to add to BlatNU.chunk are the reads that
                          # are single direction only mappers in BowtieUnique that are also single
                          # direction only mappers in BlatNU, but the two mappings disagree.
                          # Also, do not write these to RUM_Unique.
$f1 = $ARGV[2];
open(INFILE, $f1) or die "\nError: unable to open file '$f1' for reading\n\n";
while($line = <INFILE>) {
    $line =~ /^seq.(\d+)/;
    $bowtie_ambiguous_mappers{$1}++;
}
close(INFILE);
$f2 = $ARGV[0];
open(INFILE1, $f2) or die "\nError: unable to open file '$f2' for reading\n\n";
$f3 = $ARGV[1];
open(INFILE2, $f3) or die "\nError: unable to open file '$f3' for reading\n\n";
$f4 = $ARGV[4];
open(OUTFILE1, ">$f4") or die "\nError: unable to open file '$f4' for writing\n\n";

$max_distance_between_paired_reads = 500000;
$num_lines_at_once = 10000;
$linecount = 0;
$FLAG = 1;
$FLAG2 = 1;
$line_prev = <INFILE2>;
chomp($line_prev);
$last_id = 10**14;
while($FLAG == 1 || $FLAG2 == 1) {
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
	$FLAG3++;
	if($FLAG3 > 1) {
	    $FLAG2=0;
	}
    }
    until($id > $last_id || $FLAG3 > 1 || $id eq '') {
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
	    $FLAG3 = 2;
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
    if($FLAG2 == 0) {
	$line_prev = "";
    }
    foreach $id (sort {$a <=> $b} keys %allids) {
	if($bowtie_ambiguous_mappers{$id}+0 > 0) {
	    next;
	}
	if($blat_ambiguous_mappers_a{$id}+0 > 0 && $blat_ambiguous_mappers_b{$id}+0 > 0) {
	    next;
	}
	$hash1{$id}[0] = $hash1{$id}[0] + 0;
	$hash2{$id}[0] = $hash2{$id}[0] + 0;
	if(($blat_ambiguous_mappers_a{$id}+0 > 0) && ($hash1{$id}[0]+0 == 1) && ($hash1{$id}[1] =~ /seq.\d+b/)) {
	    # ambiguous foward in in BlatNU, single reverse in BowtieUnique.  See if there is
            # a consistent pairing so we can keep the pair, otherwise this read is considered unmappable
            # (not to be confused with ambiguous)
	    $line1 = $hash1{$id}[1];
	    $str = $id . "a";
	    $x = `grep $str $f0`;
	    chomp($x);
	    @a3 = split(/\n/,$x);
	    $numjoined=0;
	    for($ii=0; $ii<@a3; $ii++) {
		$line2 = $a3[$ii];
		$joined = joinifpossible($line1, $line2);
		if($joined =~ /\S/) {
		    $numjoined++;
		    $joinedsave = $joined;
		}
	    }
	    if($numjoined == 1) {  # if numjoined > 1 should probably intersect them to see if there's a 
                                   # salvagable intersection
		print OUTFILE1 "$joinedsave";
	    }
	    $remove_from_BlatNU{$id}++;
	    next;
	}
	if(($blat_ambiguous_mappers_b{$id}+0 > 0) && ($hash1{$id}[0]+0 == 1) && ($hash1{$id}[1] =~ /seq.\d+a/)) {
	    # ambiguous reverse in in BlatNU, single forward in BowtieUnique.  See if there is
            # a consistent pairing so we can keep the pair, otherwise this read is considered unmappable
            # (not to be confused with ambiguous)
	    $line1 = $hash1{$id}[1];
	    $str = $id . "b";
	    $x = `grep $str $f0`;
	    chomp($x);
	    @a3 = split(/\n/,$x);
	    $numjoined=0;
	    for($ii=0; $ii<@a3; $ii++) {
		$line2 = $a3[$ii];
		$joined = joinifpossible($line1, $line2);
		if($joined =~ /\S/) {
		    $numjoined++;
		    $joinedsave = $joined;
		}
	    }
	    if($numjoined == 1) {  # if numjoined > 1 should probably intersect them to see if there's a 
                                   # salvagable intersection
		print OUTFILE1 "$joinedsave";
	    }
	    $remove_from_BlatNU{$id}++;
	    next;
	}

# Kept for debugging
#	print "hash1{$id}[0]=$hash1{$id}[0]\n";
#	print "hash2{$id}[0]=$hash2{$id}[0]\n";

# %hash1 is bowtie
# %hash2 is blat
# These can have values -1, 0, 1, 2
# All combinations possible except (0,0), so 15 total:
# case -1: both forward and reverse reads mapped, consistently, and overlapped so were joined
# case  0: neither read mapped
# case  1: only one of the forward or reverse mapped
# case  2: both forward and reverse reads mapped, consistently, but did not overlap so were not joined

	# THREE CASES:
	if($hash1{$id}[0] == 0) {
	    if($hash2{$id}[0] == -1) {
		print OUTFILE1 "$hash2{$id}[1]\n";
	    }
	    else {
		for($i=0; $i<$hash2{$id}[0]; $i++) {    # this is in BlatUnique and not in BowtieUnique
                                                        # don't need to check if this is in BlatNU since
		                                        # it can't be in both BlatNU and BlatUnique
		    print OUTFILE1 "$hash2{$id}[$i+1]\n";
		}
	    }
	}
	# THREE CASES:
	if($hash2{$id}[0] == 0) {
	    if($hash1{$id}[0] == -1) {
		print OUTFILE1 "$hash1{$id}[1]\n";
	    }
	    if($hash1{$id}[0] == 2) {
		for($i=0; $i<$hash1{$id}[0]; $i++) {
		    print OUTFILE1 "$hash1{$id}[$i+1]\n";
		}
	    }
	    if($hash1{$id}[0] == 1) {  # this is a one-direction only mapper in BowtieUnique and
                                       # nothing in BlatUnique, so much check it's not in BlatNU
		if($blat_ambiguous_mappers_a{$id}+0 == 0 && $hash1{$id}[1] =~ /seq.\d+a/) {
		    print OUTFILE1 "$hash1{$id}[1]\n";
		}
		if($blat_ambiguous_mappers_b{$id}+0 == 0 && $hash1{$id}[1] =~ /seq.\d+b/) {
		    print OUTFILE1 "$hash1{$id}[1]\n";
		}
	    }
	}
	# ONE CASE:
	if($hash1{$id}[0] == -1 && $hash2{$id}[0] == -1) { # Preference the bowtie mapping.
	    print OUTFILE1 "$hash1{$id}[1]\n";             # This case should actually not happen because we
                                                           # should only send to blat those things which didn't
                                                           # have consistent bowtie maps.
	}
	# ONE CASE:
	if($hash1{$id}[0] == 1 && $hash2{$id}[0] == 1) {
	    if((($hash1{$id}[1] =~ /seq.\d+a/) && ($hash2{$id}[1] =~ /seq.\d+a/)) || (($hash1{$id}[1] =~ /seq.\d+b/) && ($hash2{$id}[1] =~ /seq.\d+b/))) {
	    # If single-end then this is the only case where $hash1{$id}[0] != 0 and $hash2{$id}[0] != 0
		undef @spans;
		@a1 = split(/\t/,$hash1{$id}[1]);
		@a2 = split(/\t/,$hash2{$id}[1]);
		$spans[0] = $a1[2];
		$spans[1] = $a2[2];
		$l1 = spansTotalLength($spans[0]);
		$l2 = spansTotalLength($spans[1]);
		$F=0;
		if($l1 > $l2+3) {
		    print OUTFILE1 "$hash1{$id}[1]\n";  # preference bowtie
		    $F=1;
		}
		if($l2 > $l1+3) {
		    print OUTFILE1 "$hash2{$id}[1]\n";  # preference blat
		    $F=1;
		}
		$str = intersect(\@spans, $a1[3]);
		$str =~ /^(\d+)/;
		$length_overlap = $1;
		if($F == 0) {
		    if(($length_overlap > $min_overlap) && ($a1[1] eq $a2[1])) {
			print OUTFILE1 "$hash1{$id}[1]\n";  # preference bowtie (so no worries about insertions)
		    }
		    else {
# AMBIGUOUS, OUTPUT TO NU FILE
			if(($hash1{$id}[1] =~ /\S/) && ($hash2{$id}[1] =~ /\S/)) {
			    print OUTFILE2 "$hash1{$id}[1]\n";
			    print OUTFILE2 "$hash2{$id}[1]\n";
			}
		    }
		}
	    }
	    if((($hash1{$id}[1] =~ /seq.\d+a/) && ($hash2{$id}[1] =~ /seq.\d+b/)) || (($hash1{$id}[1] =~ /seq.\d+b/) && ($hash2{$id}[1] =~ /seq.\d+a/))) {
		# This is the tricky case where there's a unique forward bowtie mapper and a unique reverse
                # blat mapper, or convsersely.  Must check for consistency.  This cannot be in BlatNU so don't
                # have to worry about that here.
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
# if they overlap, can't merge properly if there's an insertion, so chop it out,
# save it and put it back in before printing the next two if's do the chopping...
		if($aseq =~ /\+/) {
		    $aseq =~ /(.*)(\+.*\+)(.*)/;
		    $astem = $1;
		    $a_insertion = $2;
		    $apost = $3;
		    $aseq =~ s/\+.*\+//;
		    if(!($a_insertion =~ /\S/)) {
			print STDERR "Something is wrong, here 1.07\n";
		    }
		}
		if($bseq =~ /\+/) {
		    $bseq =~ /(.*)(\+.*\+)(.*)/;
		    $bstem = $1;
		    $b_insertion = $2;
		    $bpost = $3;
		    $bseq =~ s/\+.*\+//;
		    if(!($b_insertion =~ /\S/)) {
			print STDERR "Something is wrong, here 1.21\n";
		    }
		}
		$dflag = 0;
		if(($chra eq $chrb) && ($aend >= $bstart-1) && ($astart <= $bstart) && ($aend <= $bend)) {
		    # they overlap and forward is to left of reverse
		    $spans_merged = merge($aspans,$bspans);
		    $merged_length = spansTotalLength($spans_merged);
		    $aseq =~ s/://g;
		    $seq_merged = $aseq;
		    @s = split(//,$aseq);
		    $bsize = $merged_length - @s;
		    $bseq =~ s/://g;
		    @s = split(//,$bseq);
		    $add = "";
		    for($i=@s-1; $i>=@s-$bsize; $i--) {
			$add = $s[$i] . $add;
		    }
		    $seq_merged = $seq_merged . $add;
		    if($a_insertion =~ /\S/) { # put back the insertions, if any...
			$seq_merged =~ s/$astem/$astem$a_insertion/;
		    }
		    if($b_insertion =~ /\S/) {
			$str_temp = $b_insertion;
			$str_temp =~ s/\+/\\+/g;
			if(!($seq_merged =~ /$str_temp$bpost/)) {
			    $seq_merged =~ s/$bpost/$b_insertion$bpost/;
			}
		    }
		    $seq_j = addJunctionsToSeq($seq_merged, $spans_merged);
		    print OUTFILE1 "$seqnum\t$chra\t$spans_merged\t$seq_j\n";
		    $dflag = 1;
		}
		if(($chra eq $chrb) && ($bend >= $astart-1) && ($bstart <= $astart) && ($bend <= $aend) && ($dflag == 0)) {
		    # they overlap and reverse is to left of forward
		    $spans_merged = merge($bspans,$aspans);
		    $merged_length = spansTotalLength($spans_merged);
		    $bseq =~ s/://g;
		    $seq_merged = $bseq;
		    @s = split(//,$bseq);
		    $asize = $merged_length - @s;
		    $aseq =~ s/://g;
		    @s = split(//,$aseq);
		    $add = "";
		    for($i=@s-1; $i>=@s-$asize; $i--) {
			$add = $s[$i] . $add;
		    }
		    $seq_merged = $seq_merged . $add;
		    if($a_insertion =~ /\S/) { # put back the insertions, if any...
			$seq_merged =~ s/$apost/$a_insertion$apost/;
		    }
		    if($b_insertion =~ /\S/) {
			$str_temp = $b_insertion;
			$str_temp =~ s/\+/\\+/g;
			if(!($seq_merged =~ /$bstem$str_temp/)) {
			    $seq_merged =~ s/$bstem/$bstem$b_insertion/;
			}
		    }
		    $seq_j = addJunctionsToSeq($seq_merged, $spans_merged);
		    print OUTFILE1 "$seqnum\t$chra\t$spans_merged\t$seq_j\n";
		}
	    }
	}
	# ONE CASE
	if($hash1{$id}[0] == 2 && $hash2{$id}[0] == 2) { # preference bowtie
	    print OUTFILE1 "$hash1{$id}[1]\n";
	    print OUTFILE1 "$hash1{$id}[2]\n";
	}	
	# NINE CASES DONE
	# ONE CASE
	if($hash1{$id}[0] == -1 && $hash2{$id}[0] == 2) { # preference bowtie
	    print OUTFILE1 "$hash1{$id}[1]\n";
	}
	# ONE CASE
	if($hash1{$id}[0] == 2 && $hash2{$id}[0] == -1) { # preference bowtie
	    print OUTFILE1 "$hash1{$id}[1]\n";
	    print OUTFILE1 "$hash1{$id}[2]\n";
	}
	# ELEVEN CASES DONE
	if($hash1{$id}[0] == 1 && $hash2{$id}[0] == 2) {
	    print OUTFILE1 "$hash2{$id}[1]\n";
	    print OUTFILE1 "$hash2{$id}[2]\n";
	}	
	if($hash1{$id}[0] == 2 && $hash2{$id}[0] == 1) {
	    print OUTFILE1 "$hash1{$id}[1]\n";
	    print OUTFILE1 "$hash1{$id}[2]\n";
	}	
	if($hash1{$id}[0] == -1 && $hash2{$id}[0] == 1) {
	    print OUTFILE1 "$hash1{$id}[1]\n";
	}
	if($hash1{$id}[0] == 1 && $hash2{$id}[0] == -1) {
	    print OUTFILE1 "$hash2{$id}[1]\n";
	}
	# ALL FIFTEEN CASES DONE
    }
}
close(INFILE1);
close(INFILE2);
close(OUTFILE1);
close(OUTFILE2);

# now need to remove the stuff in %remove_from_BlatNU from BlatNU
$filename = $ARGV[5];
open(INFILE, $f0) or die "\nError: unable to open file '$f0' for reading\n\n";
open(OUTFILE, ">$filename") or die "\nError: unable to open file '$filename' for writing\n\n";
while($line = <INFILE>) {
    $line =~ /seq.(\d+)/;
    if($remove_from_BlatNU{$1}+0==0) {
	print OUTFILE $line;
    }
}
close(INFILE);
open(INFILE, $f1) or die "\nError: unable to open file '$f1' for reading\n\n";
                  # now append BowtieNU to get the full NU file
while($line = <INFILE>) {
    print OUTFILE $line;
}
close(INFILE);
close(OUTFILE);

sub joinifpossible () {
    ($line1, $line2) = @_;
    @a = split(/\t/,$line1);
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
    @a = split(/\t/,$line2);
    $bspans = $a[2];
    $a[2] =~ /^(\d+)[^\d]/;
    $bstart = $1;
    $a[2] =~ /[^\d](\d+)$/;
    $bend = $1;
    $chrb = $a[1];
    $bseq = $a[3];
    $returnstring = "";
    if(($chra eq $chrb) && ($aend < $bstart-1) && ($bstart - $aend < $max_distance_between_paired_reads)) {
	if($line1 =~ /a\t/) {
	    $returnstring = $returnstring . "$line1\n$line2\n";
	}
	else {
	    $returnstring = $returnstring . "$line2\n$line1\n";
	}
    }
    if(($chra eq $chrb) && ($bend < $astart-1) && ($astart - $bend < $max_distance_between_paired_reads)) {
	if($line1 =~ /a\t/) {
	    $returnstring = $returnstring . "$line1\n$line2\n";
	}
	else {
	    $returnstring = $returnstring . "$line2\n$line1\n";
	}
    }
# if they overlap, can't merge properly if there's an insertion, so chop it out,
# save it and put it back in before printing the next two if's do the chopping...
    if($aseq =~ /\+/) {
	$aseq =~ /(.*)(\+.*\+)(.*)/;
	$astem = $1;
	$a_insertion = $2;
	$apost = $3;
	$aseq =~ s/\+.*\+//;
	if(!($a_insertion =~ /\S/)) {
	    $returnstring = $returnstring . "Something is wrong, here 1.07\n";
	}
    }
    if($bseq =~ /\+/) {
	$bseq =~ /(.*)(\+.*\+)(.*)/;
	$bstem = $1;
	$b_insertion = $2;
	$bpost = $3;
	$bseq =~ s/\+.*\+//;
	if(!($b_insertion =~ /\S/)) {
	    $returnstring = $returnstring . "Something is wrong, here 1.21\n";
	}
    }
    $dflag = 0;
    if(($chra eq $chrb) && ($aend >= $bstart-1) && ($astart <= $bstart) && ($aend <= $bend)) {
	# they overlap and forward is to left of reverse
	$spans_merged = merge($aspans,$bspans);
	$merged_length = spansTotalLength($spans_merged);
	$aseq =~ s/://g;
	$seq_merged = $aseq;
	@s = split(//,$aseq);
	$bsize = $merged_length - @s;
	$bseq =~ s/://g;
	@s = split(//,$bseq);
	$add = "";
	for($i=@s-1; $i>=@s-$bsize; $i--) {
	    $add = $s[$i] . $add;
	}
	$seq_merged = $seq_merged . $add;
	if($a_insertion =~ /\S/) { # put back the insertions, if any...
	    $seq_merged =~ s/$astem/$astem$a_insertion/;
	}
	if($b_insertion =~ /\S/) {
	    $str_temp = $b_insertion;
	    $str_temp =~ s/\+/\\+/g;
	    if(!($seq_merged =~ /$str_temp$bpost/)) {
		$seq_merged =~ s/$bpost/$b_insertion$bpost/;
	    }
	}
	$seq_j = addJunctionsToSeq($seq_merged, $spans_merged);
	$returnstring = $returnstring . "$seqnum\t$chra\t$spans_merged\t$seq_j\n";
	$dflag = 1;
    }
    if(($chra eq $chrb) && ($bend >= $astart-1) && ($bstart <= $astart) && ($bend <= $aend) && ($dflag == 0)) {
	# they overlap and reverse is to left of forward
	$spans_merged = merge($bspans,$aspans);
	$merged_length = spansTotalLength($spans_merged);
	$bseq =~ s/://g;
	$seq_merged = $bseq;
	@s = split(//,$bseq);
	$asize = $merged_length - @s;
	$aseq =~ s/://g;
	@s = split(//,$aseq);
	$add = "";
	for($i=@s-1; $i>=@s-$asize; $i--) {
	    $add = $s[$i] . $add;
	}
	$seq_merged = $seq_merged . $add;
	if($a_insertion =~ /\S/) { # put back the insertions, if any...
	    $seq_merged =~ s/$apost/$a_insertion$apost/;
	}
	if($b_insertion =~ /\S/) {
	    $str_temp = $b_insertion;
	    $str_temp =~ s/\+/\\+/g;
	    if(!($seq_merged =~ /$bstem$str_temp/)) {
		$seq_merged =~ s/$bstem/$bstem$b_insertion/;
	    }
	}
	$seq_j = addJunctionsToSeq($seq_merged, $spans_merged);
	$returnstring = $returnstring . "$seqnum\t$chra\t$spans_merged\t$seq_j\n";
    }
    return $returnstring;
}

sub merge () {
    ($aspans2, $bspans2) = @_;
    undef @astarts2;
    undef @aends2;
    undef @bstarts2;
    undef @bends2;
    @a = split(/, /, $aspans2);
    for($i=0; $i<@a; $i++) {
	@b = split(/-/,$a[$i]);
	$astarts2[$i] = $b[0];
	$aends2[$i] = $b[1];
    }
    @a = split(/, /, $bspans2);
    for($i=0; $i<@a; $i++) {
	@b = split(/-/,$a[$i]);
	$bstarts2[$i] = $b[0];
	$bends2[$i] = $b[1];
    }
    if($aends2[@aends2-1] + 1 < $bstarts2[0]) {
	$merged_spans = $aspans2 . ", " . $bspans2;
    }
    if($aends2[@aends2-1] + 1 == $bstarts2[0]) {
	$aspans2 =~ s/-\d+$//;
	$bspans2 =~ s/^\d+-//;
	$merged_spans = $aspans2 . "-" . $bspans2;
    }
    if($aends2[@aends2-1] + 1 > $bstarts2[0]) {
	$merged_spans = $aspans2;
	for($i=0; $i<@bstarts2; $i++) {
	    if($aends2[@aends2-1] >= $bstarts2[$i] && ($aends2[@aends2-1] <= $bstarts2[$i+1] || $i == @bstarts2-1)) {
		$merged_spans =~ s/-\d+$//;
		$merged_spans = $merged_spans . "-" . $bends2[$i];
		for($j=$i+1; $j<@bstarts2; $j++) {
		    $merged_spans = $merged_spans . ", $bstarts2[$j]-$bends2[$j]";
		}
	    }
	}
    }
    return $merged_spans;
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
    $seq =~ s/://g;
    @s = split(//,$seq);
    @b = split(/, /,$spans);
    $seq_out = "";
    $place = 0;
    for($j=0; $j<@b; $j++) {
	@c = split(/-/,$b[$j]);
	$len = $c[1] - $c[0] + 1;
	if($seq_out =~ /\S/) { # to avoid putting a colon at the beginning
	    $seq_out = $seq_out . ":";
	}
	for($k=0; $k<$len; $k++) {
	    if($s[$place] eq "+") {
		$seq_out = $seq_out . $s[$place];
		$place++;
		until($s[$place] eq "+") {
		    $seq_out = $seq_out . $s[$place];
		    $place++;
		    if($place > @s-1) {
			last;
		    }
		}
		$k--;
	    }
	    $seq_out = $seq_out . $s[$place];
	    $place++;
	}
    }
    return $seq_out;
}

sub spansTotalLength () {
    ($spans) = @_;
    @a = split(/, /,$spans);
    $length = 0;
    for($i=0; $i<@a; $i++) {
	@b = split(/-/,$a[$i]);
	$length = $length + $b[1] - $b[0] + 1;
    }
    return $length;
}