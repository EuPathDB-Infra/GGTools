#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

$|=1;

if(@ARGV < 5) {
    die "
Usage: rum2sam.pl <rum unique file> <rum nu file> <reads file> <quals file> <sam outfile> [options]

<reads file> and <qual file> are the files coming from parse2fasta.pl and fastq2qualities.pl
if you don't have the qualities just put 'none' for the <quals file> argument.

";
}

$rum_unique_file = $ARGV[0];
$rum_nu_file = $ARGV[1];
$reads_file = $ARGV[2];
$qual_file = $ARGV[3];
$quals = "true";
if($ARGV[3] =~ /none/ || $ARGV[3] =~ /.none./) {
    $quals = "false";
}
$sam_outfile = $ARGV[4];

open(INFILE, $reads_file);
$line = <INFILE>;
chomp($line);
$line =~ /seq.(\d+)/;
$firstseqnum = $1;
$line = <INFILE>;
chomp($line);
$readlength = length($line);
if($quals eq "false") {
    $QUAL = "";
    for($i=0; $i<$readlength; $i++) {
	$QUAL = $QUAL . "I";
    }
}
$line = <INFILE>;
chomp($line);
$line =~ /seq.\d+(.)/;
$type = $1;
if($type eq 'b') {
    $paired = "true";
} else {
    $paired = "false";
}
close(INFILE);
$x = `tail -2 $reads_file | head -1`;
$x =~ /seq.(\d+)/;
$lastseqnum = $1;

$bitflag[0] = "the read is paired in sequencing";
$bitflag[1] = "the read is mapped in a proper pair";
$bitflag[2] = "the query sequence itself is unmapped";
$bitflag[3] = "the mate is unmapped";
$bitflag[4] = "strand of the query";
$bitflag[5] = "strand of the mate";
$bitflag[6] = "the read is the first read in a pair";
$bitflag[7] = "the read is the second read in a pair";
$bitflag[8] = "the alignment is not primary";
$bitflag[9] = "the read fails platform/vendor quality checks";
$bitflag[10] = "the read is either a PCR duplicate or an optical duplicate";

open(RUMU, $rum_unique_file) or die "\nError: cannot open the file '$rum_unique_file' for reading\n\n";
open(RUMNU, $rum_nu_file) or die "\nError: cannot open the file '$rum_nu_file' for reading\n\n";
open(READS, $reads_file) or die "\nError: cannot open the file '$reads_file' for reading\n\n";

# checkin that the first line in RUMU really looks like it should:
$line = <RUMU>;
close(RUMU);
@a = split(/\t/,$line);
$flag = 0;
if(!($a[0] =~ /^seq.\d+[ab]?/)) {
    $flag = 1;
}
if($a[2] =~ /[^\d-, ]/) {
    $flag = 1;
}
if(!($a[3] eq "+" || $a[3] eq "-")) {
    $flag = 1;
}
if(!($a[4] =~ /^[ACGTN:+]+$/)) {
    $flag = 1;
}
if($flag == 1) {
    die "\nError: the first line of the file '$rum_unique_file' is misformatted,\nit does not look like a RUM output file.\n";
}
close(RUMU);
open(RUMU, $rum_unique_file) or die "\nError: cannot open the file '$rum_unique_file' for reading\n\n";

if($quals eq "true") {
    open(QUALS, $qual_file);
}
open(SAM, ">$sam_outfile");

for($seqnum = $firstseqnum; $seqnum <= $lastseqnum; $seqnum++) {
    $forward_read = <READS>;
    $forward_read = <READS>;
    chomp($forward_read);
    if($paired eq "true") {
	$reverse_read = <READS>;
	$reverse_read = <READS>;
	chomp($reverse_read);
    }
    if($quals eq "true") {
	$forward_qual = <QUALS>;
	$forward_qual = <QUALS>;
	chomp($forward_qual);
	if($paired eq "true") {
	    $reverse_qual = <QUALS>;
	    $reverse_qual = <QUALS>;
	    chomp($reverse_qual);
	}
    } else {
	$forward_qual = $QUAL;
	$reverse_qual = $QUAL;
    }

    $flag = 0;
    $rum_u_forward = "";
    $rum_u_reverse = "";
    $rum_u_joined = "";
    $unique_mapper_found = "false";
    $non_unique_mappers_found = "false";
    while($flag == 0) {
	$line = <RUMU>;
	chomp($line);
	$type = "";
	if($line =~ /seq.(\d+)(.)/) {
	    $sn = $1;
	    $type = $2;
	}
	if($sn == $seqnum && $type eq "a") {
	    $rum_u_forward = $line;
	    $unique_mapper_found = "true";
	}
	if($sn == $seqnum && $type eq "b") {
	    $rum_u_reverse = $line;
	    $unique_mapper_found = "true";
	}
	if($sn == $seqnum && $type eq "\t") {
	    $rum_u_joined = $line;
	    $unique_mapper_found = "true";
	}
	if($sn > $seqnum) {
	    $len = -1 * (1 + length($line));
	    seek(RUMU, $len, 1);
	    $flag = 1;
	}
	if($line eq '') {
	    $flag = 1;
	}
    }

    if($unique_mapper_found eq "true") {

	# SET THE BITSCORE

	$bitscore_f = 0;
	$bitscore_r = 0;
	if($paired eq "true") {
	    $bitscore_f = 65;
	    $bitscore_r = 129;
	    if(!($rum_u_joined =~ /\S/)) {
		if(!($rum_u_forward =~ /\S/) && !($rum_u_reverse =~ /\S/)) {
		    $bitscore_r = $bitscore_r + 12;
		    $bitscore_f = $bitscore_f + 12;
		}
		if($rum_u_forward =~ /\S/ && !($rum_u_reverse =~ /\S/)) {
		    $bitscore_r = $bitscore_r + 4;
		    $bitscore_f = $bitscore_f + 8;
		}
		if($rum_u_reverse =~ /\S/ && !($rum_u_forward =~ /\S/)) {
		    $bitscore_f = $bitscore_f + 4;
		    $bitscore_r = $bitscore_r + 8;
		}
	    }
	} else {
	    $bitscore = 0;
	}
	if(($rum_u_forward =~ /\S/ && $rum_u_reverse =~ /\S/) || $rum_u_joined =~ /\S/) {
	    $bitscore_f = $bitscore_f + 2;
	    $bitscore_r = $bitscore_r + 2;
	}

	$joined = "false";
	if($rum_u_joined =~ /\S/) {
	    # FORWARD AND REVERSE MAPPED, AND THEY ARE JOINED, GATHER INFO
	    $joined = "true";
#	    print "---------------\nrum_u_joined = $rum_u_joined\n";
	    undef @piecelength;
	    @ruj = split(/\t/,$rum_u_joined);
	    $ruj[4] =~ s/://g;
	    @PL = split(/\+/,$ruj[4]);
	    $piecelength[0] = length($PL[0]);
	    for($pl=1; $pl<@PL; $pl++) {
		$piecelength[$pl] = length($PL[$pl]) + $piecelength[$pl-1];
	    }
	    @ruj = split(/\t/,$rum_u_joined);
	    if($ruj[3] eq "-") {
		$upstream_read = $reverse_read;
	    } else {
		$upstream_read = $forward_read;
	    }
	    $x = $upstream_read;
	    $ruj[4] =~ s/://g;
	    $ruj[4] =~ s/\+//g;
	    $y = $ruj[4];
	    $prefix_offset_upstream = 0;
	    $L = length($x);
	    $count=0;
	    $suffix_offset_upstream = 0;
	    $LEN = 0;
	    $LENflag = 0;
	    $LEN_current_best=0;
	    while($LENflag == 0) {
		$LENflag = 1;
		until($y =~ /^$x/) {
		    $x =~ s/^.//;
		    $count++;
		    $prefix_offset_upstream++;
		}
		$LEN = $L - $count;
		if($LEN >= $LEN_current_best) {
		    $suffix_offset_upstream_current_best = $suffix_offset_upstream;
		    $prefix_offset_upstream_current_best = $prefix_offset_upstream;
		    $LEN_current_best = $LEN;
		}
		if($LEN < $readlength / 2) {
		    $LENflag = 0;
		    $x = $upstream_read;
		    $suffix_offset_upstream++;
		    for($j=0; $j<$suffix_offset_upstream; $j++) {
			$x =~ s/.$//;
		    }
		    $prefix_offset_upstream = 0;
		    $count = 0;
		    $L = length($x);
		    if($L < 1) {
			last;
		    }
		}
	    }

	    $prefix_offset_upstream = $prefix_offset_upstream_current_best;
	    $suffix_offset_upstream = $suffix_offset_upstream_current_best;
    
	    $UR = $upstream_read;
	    $replace = "";
	    for($i=0; $i<$suffix_offset_upstream; $i++) {
		$UR =~ s/.$//;
		$replace = $replace . "X";
	    }
	    $UR2 = $UR;
	    $UR = $UR . $replace;

	    $plen = $readlength - $prefix_offset_upstream - $suffix_offset_upstream;
	    $pl=0;
	    $RC = 0;
	    $matchlength = $piecelength[0];
#	    print "suffix_offset_upstream = $suffix_offset_upstream\n";
#	    print "UR2 = $UR2\n";
	    while($piecelength[$pl] + $prefix_offset_upstream < $readlength - $suffix_offset_upstream) {
		$plen = $plen - ($piecelength[$pl+1] - $piecelength[$pl]);
		if($piecelength[$pl+1] > $readlength) { # insertion went past the end of the read,
		                                        # so overcorrected, this fixes that
		    $plen = $plen + ($piecelength[$pl+1] - $readlength);
		}
		substr($UR2, $piecelength[$pl]+$RC+$prefix_offset_upstream, 0, "+");
		$RC++;
		$XX = $piecelength[$pl+1]+$RC+$prefix_offset_upstream;
		$YY = length($UR2);
		if($XX <= $YY) {
		    substr($UR2, $piecelength[$pl+1]+$RC+$prefix_offset_upstream, 0, "+");
		} else { # individual alignments don't have insertions at the ends,
                         # so removing it because it'll mess things up downstream
		    if($UR2 =~ s/\+([^\+]+)$//) {
			$suffix_offset_upstream = $suffix_offset_upstream + length($1);
		    }
		}
		if($UR2 =~ s/\+([^\+]+)\+$//) { # just in case there's still an insertion at the end...
		    $suffix_offset_upstream = $suffix_offset_upstream + length($1);
		}
		$RC++;
		$pl=$pl+2;
		$matchlength = $matchlength + $piecelength[$pl] - $piecelength[$pl-1];
	    }
#	    print "UR2 = $UR2\n";
#	    print "suffix_offset_upstream = $suffix_offset_upstream\n";
	    for($i=0; $i<$prefix_offset_upstream; $i++) {
		$UR2 =~ s/^.//;
	    }
	    $upstream_spans = &getprefix($ruj[2], $plen);

	    if($ruj[3] eq "-") {
		$downstream_read = reversecomplement($forward_read);
		$bitscore_f = $bitscore_f + 16;
		$bitscore_r = $bitscore_r + 32;
	    } else {
		$downstream_read = reversecomplement($reverse_read);
		$bitscore_r = $bitscore_r + 16;
		$bitscore_f = $bitscore_f + 32;
	    }
	    $x = $downstream_read;
	    $y = $ruj[4];
	    $suffix_offset_downstream = 0;
	    $L = length($x);
	    $count=0;
	    $prefix_offset_downstream = 0;
	    $LEN = 0;
	    $LENflag = 0;
	    $LEN_current_best=0;
	    while($LENflag == 0) {
		$LENflag = 1;
		until($y =~ /$x$/ || length($x)==0) {
		    $x =~ s/.$//;
		    $count++;
		    $suffix_offset_downstream++;
		}
		$LEN = $L - $count;
		if($LEN >= $LEN_current_best) {
		    $suffix_offset_downstream_current_best = $suffix_offset_downstream;
		    $prefix_offset_downstream_current_best = $prefix_offset_downstream;
		    $LEN_current_best = $LEN;
		}
		if($LEN < $readlength / 2) {
		    $LENflag = 0;
		    $x = $downstream_read;
		    $prefix_offset_downstream++;
		    for($j=0; $j<$prefix_offset_downstream; $j++) {
			$x =~ s/^.//;
		    }
		    $suffix_offset_downstream = 0;
		    $count = 0;
		    $L = length($x);
		    if($L < 1) {
			last;
		    }
		}
	    }

	    $prefix_offset_downstream = $prefix_offset_downstream_current_best;
	    $suffix_offset_downstream = $suffix_offset_downstream_current_best;

	    $DR = $downstream_read;
	    $replace = "";
	    for($i=0; $i<$prefix_offset_downstream; $i++) {
		$DR =~ s/^.//;
		$replace = $replace . "X";
	    }
	    $DR2 = $DR;
	    $DR = $replace . $DR;
#	    print "DR2=$DR2\n";

	    $offset = length($ruj[4]) + $prefix_offset_upstream + $suffix_offset_downstream - length($DR);

	    $OFFSET = $readlength - length($ruj[4]) - $suffix_offset_downstream;
	    $P = "";
	    if($OFFSET < 0) {
		$OFFSET = 0;
	    }
	    for($i=0; $i<$OFFSET; $i++) {
		$P = $P . " ";
	    }
	    $plen = $readlength - $prefix_offset_downstream - $suffix_offset_downstream;
#	    print "plen=$plen\n";

# 	    print "\n$upstream_read\n\n";
#	    print $P;
#	    print "$UR\n";
#	    print $P;
#	    for($i=0; $i<$prefix_offset_upstream; $i++) {
#		print " ";
#	     }
#	    print "$ruj[4]\n";
#	    for($i=0; $i<$offset; $i++) {
#		print " ";
#	    }
#	    print "$DR\n";
#	    print "\n$downstream_read\n\n";



#	    print "piecelength[0] = $piecelength[0]\n";
#	    print "piecelength[1] = $piecelength[1]\n";
#	    print "piecelength[2] = $piecelength[2]\n";
#	    print "piecelength[3] = $piecelength[3]\n";
#	    print "piecelength[4] = $piecelength[4]\n";

	    $RC = 0;
	    $pl=0;
#	    print "offset=$offset\n";
#	    print "prefix_offset_downstream=$prefix_offset_downstream\n";
#	    print "prefix_offset_upstream=$prefix_offset_upstream\n";
	    until($piecelength[$pl] > $offset + $prefix_offset_downstream || $pl >= @piecelength) {
		$pl++;
	    }
#	    print "pl=$pl\n";
	    # the first three if's here deal with the case that there's an insertion right at 
	    # the begginging of the downstream read, either starting at the start of the read,
	    # or ending just before it, or overlapping the end.
	    if($pl == 0 && $piecelength[0] == $offset) {
		substr($DR2, $piecelength[$pl+1]-$piecelength[$pl], 0, "+");
		$DR2 = "+" . $DR2;
	    } elsif(($pl == 1 && $piecelength[1] == $offset) || ($pl >= @piecelength)) {
		# do nothing
	    }
	    elsif($pl % 2 == 1) {
		substr($DR2, $piecelength[$pl]-$offset+$RC-$prefix_offset_downstream+$prefix_offset_upstream, 0, "+");
		$pl--;
		$DR2 = "+" . $DR2;
	    } else {
#		print "pl=$pl\n";
#		print "DR2=$DR2\n";
		while($piecelength[$pl] >= $offset + $prefix_offset_downstream -$prefix_offset_upstream && $pl < @piecelength-1) {
		    $plen = $plen - ($piecelength[$pl+1] - $piecelength[$pl]);
#		    print "x:plen=$plen\n";
		    $ABC = $piecelength[$pl]-$offset+$RC-$prefix_offset_downstream+$prefix_offset_upstream;
		    $XYZ = $piecelength[$pl+1]-$offset+$RC-$prefix_offset_downstream+$prefix_offset_upstream;
		    substr($DR2, $piecelength[$pl]-$offset+$RC-$prefix_offset_downstream+$prefix_offset_upstream, 0, "+");
		    $RC++;
		    substr($DR2, $piecelength[$pl+1]-$offset+$RC-$prefix_offset_downstream+$prefix_offset_upstream, 0, "+");
		    $RC++;
		    $pl=$pl+2;
		}
	    }
	    $DR2 =~ s/^\+([^\+]+)\+//; # individual alignments don't have insertions at the ends,
	                               # so removing it because it'll mess things up downstream
	    $prefix_offset_downstream = $prefix_offset_downstream + length($1);
#	    print "plen=$plen\n";
	    $plen = $plen - length($1);
#	    print "plen=$plen\n";

	    for($i=0; $i<$suffix_offset_downstream; $i++) {
		$DR2 =~ s/.$//;
	    }

	    $downstream_spans = &getsuffix($ruj[2], $plen);

	    if($ruj[3] eq "+") {

		$rum_u_forward = $seqnum . "a\t$ruj[1]\t" . $upstream_spans . "\t+\t" . $UR2;
		$rum_u_reverse = $seqnum . "b\t$ruj[1]\t" . $downstream_spans . "\t+\t" . $DR2;
	    }
	    if($ruj[3] eq "-") {

		$rum_u_forward = $seqnum . "a\t$ruj[1]\t" . $downstream_spans . "\t-\t" . $DR2;
		$rum_u_reverse = $seqnum . "b\t$ruj[1]\t" . $upstream_spans . "\t-\t" . $UR2;
	    }

	}

	if($rum_u_forward =~ /\S/) {

	    # COLLECT INFO FROM FORWARD RUM RECORD
	    # note: this might be a joined read for which the surrogate forward was created above

	    undef @piecelength;
#	    if($rum_u_joined =~ /\S/) {
# 		print "rum_u_forward = $rum_u_forward\n\n";
#	    }
	    @ruf = split(/\t/,$rum_u_forward);
	    $ruf[4] =~ s/://g;
	    @PL = split(/\+/,$ruf[4]);
	    $piecelength[0] = length($PL[0]);
	    $insertions_total_length = 0;
	    for($pl=1; $pl<@PL; $pl=$pl+2) {
		$insertions_total_length = $insertions_total_length + length($PL[$pl]);
	    }
	    for($pl=1; $pl<@PL; $pl++) {
		$piecelength[$pl] = length($PL[$pl]) + $piecelength[$pl-1];
	    }
	    $ruf[4] =~ s/\+//g;
	    $rum_u_forward_length = length($ruf[4]);
	    if($ruf[3] eq "-") {
		$forward_read = reversecomplement($forward_read);
		if(!($rum_u_joined =~ /\S/)) {
		    $bitscore_f = $bitscore_f + 16;
		    $bitscore_r = $bitscore_r + 32;
		    if(!($rum_u_reverse =~ /\S/)) {
			$bitscore_r = $bitscore_r + 16;
		    }
		}
	    }
	    if($rum_u_joined =~ /\S/) {
		if($ruf[3] eq "+") {
		    $prefix_offset_forward = $prefix_offset_upstream;
		    $suffix_offset_forward = $suffix_offset_upstream;
		} else {
		    $prefix_offset_forward = $prefix_offset_downstream;
		    $suffix_offset_forward = $suffix_offset_downstream;
		}
	    } else {
		$prefix_offset_forward = 0;
		if($rum_u_forward_length < $readlength) {
		    $x = $forward_read;
		    $y = $ruf[4];
		    until($x =~ /^$y/) {
			$x =~ s/^.//;
			$prefix_offset_forward++;
		    }
		}
	    }

	    $CIGAR_f = "";
	    $insertions_finished = 0;
	    if($prefix_offset_forward > 0) {
		$CIGAR_f = $prefix_offset_forward . "S";
	    }
	    @aspans = split(/, /,$ruf[2]);
	    @C1 = split(/-/,$aspans[0]);
	    $L = $C1[1] - $C1[0] + 1;
	    $running_length = 0;
	    # code for insertions follows
	    if($L > $piecelength[0]) {
		$pref_length = $piecelength[0] - $running_length;
		$insertion_length = $piecelength[1] - $piecelength[0];
		$suff_length = $L - $piecelength[0];
		$CIGAR_f = $CIGAR_f . $pref_length . "M" . $insertion_length . "I" . $suff_length . "M";
		$running_length = $running_length + $insertion_length;
		$insertions_finished++;
	    } else {
		$CIGAR_f = $CIGAR_f . $L . "M";
	    }
	    $running_length = $running_length + $L;
	    for($i=1; $i<@aspans; $i++) {
		@C2 = split(/-/,$aspans[$i]);
		$skipped = $C2[0] - $C1[1] - 1;
		if($skipped >= 15) {
		    $CIGAR_f = $CIGAR_f . $skipped . "N";
		} else {
		    $CIGAR_f = $CIGAR_f . $skipped . "D";
		}
		$L = $C2[1] - $C2[0] + 1;
		# code for insertions follows
		if($running_length+$L > $piecelength[$insertions_finished*2]) {
		    $pref_length = $piecelength[$insertions_finished*2] - $running_length;
		    $insertion_length = $piecelength[$insertions_finished*2+1] - $piecelength[$insertions_finished*2];
		    $running_length = $running_length + $insertion_length;
		    $suff_length = $running_length + $L - $piecelength[$insertions_finished*2+1];
		    $CIGAR_f = $CIGAR_f . $pref_length . "M" . $insertion_length . "I" . $suff_length . "M";
		    $insertions_finished++;
		} else {
		    $CIGAR_f = $CIGAR_f . $L . "M";
		}
		$running_length = $running_length + $L;
		$C1[0] = $C2[0];
		$C1[1] = $C2[1];
	    }

	    $right_clip_size_f = $readlength - $running_length - $prefix_offset_forward;
	    if($right_clip_size_f > 0) {
		if($rum_u_forward =~ /\+$/) {
		    $CIGAR_f = $CIGAR_f . $right_clip_size_f . "I";
		} else {
		    $CIGAR_f = $CIGAR_f . $right_clip_size_f . "S";
		}
	    }
	}


	if($rum_u_reverse =~ /\S/) {

	    # COLLECT INFO FROM REVERSE RUM RECORD
	    # note: this might be a joined read for which the surrogate forward was created above

	    undef @piecelength;
#	    if($rum_u_joined =~ /\S/) {
#		print "rum_u_reverse = $rum_u_reverse\n\n";
#	    }
	    @rur = split(/\t/,$rum_u_reverse);
	    $rur[4] =~ s/://g;
	    @PL = split(/\+/,$rur[4]);
	    $piecelength[0] = length($PL[0]);
	    $insertions_total_length = 0;
	    for($pl=1; $pl<@PL; $pl=$pl+2) {
		$insertions_total_length = $insertions_total_length + length($PL[$pl]);
	    }
	    for($pl=1; $pl<@PL; $pl++) {
		$piecelength[$pl] = length($PL[$pl]) + $piecelength[$pl-1];
	    }
	    $rur[4] =~ s/\+//g;
	    $rum_u_reverse_length = length($rur[4]);
	    if($rur[3] eq "+") {
		$reverse_read = reversecomplement($reverse_read);
		if(!($rum_u_joined =~ /\S/)) {
		    $bitscore_r = $bitscore_r + 16;
		    $bitscore_f = $bitscore_f + 32;
		    if(!($rum_u_forward =~ /\S/)) {
			$bitscore_f = $bitscore_f + 16;
		    }
		}
	    }

#	    print "$reverse_read\n";
	    if($rum_u_joined =~ /\S/) {
		if($ruf[3] eq "+") {
		    $prefix_offset_reverse = $prefix_offset_downstream;
		    $suffix_offset_reverse = $suffix_offset_downstream;
		} else {
		    $prefix_offset_reverse = $prefix_offset_upstream;
		    $suffix_offset_reverse = $suffix_offset_upstream;
		}
	    } else {
		$prefix_offset_reverse = 0;
		if($rum_u_reverse_length < $readlength) {
		    $x = $reverse_read;
		    $y = $rur[4];
		    until($x =~ /^$y/) {
			$x =~ s/^.//;
			$prefix_offset_reverse++;
#		    print " ";
		    }
		}
	    }
#	    print "$rur[4]\n\n";

	    $CIGAR_r = "";
	    $insertions_finished = 0;
	    if($prefix_offset_reverse > 0) {
		$CIGAR_r = $prefix_offset_reverse . "S";
	    }
	    @bspans = split(/, /,$rur[2]);
	    @C1 = split(/-/,$bspans[0]);
	    $L = $C1[1] - $C1[0] + 1;
	    $running_length = 0;
	    # code for insertions follows
#	    print "piecelength[$insertions_finished*2] = $piecelength[$insertions_finished*2]\n";
	    if($L > $piecelength[0]) {
		$pref_length = $piecelength[0] - $running_length;
		$insertion_length = $piecelength[1] - $piecelength[0];
		$suff_length = $L - $piecelength[0];
		$CIGAR_r = $CIGAR_r . $pref_length . "M" . $insertion_length . "I" . $suff_length . "M";
		$running_length = $running_length + $insertion_length;
		$insertions_finished++;
	    } else {
		$CIGAR_r = $CIGAR_r . $L . "M";
	    }
#	    print "running_length = $running_length\n";
	    $running_length = $running_length + $L;
#	    print "running_length = $running_length\n";
	    for($i=1; $i<@bspans; $i++) {
		@C2 = split(/-/,$bspans[$i]);
		$skipped = $C2[0] - $C1[1] - 1;
		if($skipped >= 15) {
		    $CIGAR_r = $CIGAR_r . $skipped . "N";
		} else {
		    $CIGAR_r = $CIGAR_r . $skipped . "D";
		}
		$L = $C2[1] - $C2[0] + 1;
		# code for insertions follows
#		print "piecelength[$insertions_finished*2] = $piecelength[$insertions_finished*2]\n";
		if($running_length+$L > $piecelength[$insertions_finished*2]) {
		    $pref_length = $piecelength[$insertions_finished*2] - $running_length;
		    $insertion_length = $piecelength[$insertions_finished*2+1] - $piecelength[$insertions_finished*2];
		    $running_length = $running_length + $insertion_length;
		    $suff_length = $running_length + $L - $piecelength[$insertions_finished*2+1];
		    $CIGAR_r = $CIGAR_r . $pref_length . "M" . $insertion_length . "I" . $suff_length . "M";
		    $insertions_finished++;
		} else {
		    $CIGAR_r = $CIGAR_r . $L . "M";
		}
		$running_length = $running_length + $L;
		$C1[0] = $C2[0];
		$C1[1] = $C2[1];
	    }
#	    print "readlength = $readlength\n";
#	    print "running_length = $running_length\n";
#	    print "prefix_offset_reverse = $prefix_offset_reverse\n";
	    $right_clip_size_r = $readlength - $running_length - $prefix_offset_reverse;
	    if($right_clip_size_r > 0) {
		if($rum_u_reverse =~ /\+$/) {
		    $CIGAR_r = $CIGAR_r . $right_clip_size_r . "I";
		} else {
		    $CIGAR_r = $CIGAR_r . $right_clip_size_r . "S";
		}
	    }
	}

#	print "bitscore_f = $bitscore_f\n";
#	print "bitscore_r = $bitscore_r\n";

# COMPUTE IDIST

	$idist_f = 0;
	$idist_r = 0;

	if($ruf[2] =~ /^(\d+)-/) {
	    $start_forward = $1;
	} else {
	    $start_forward = "*";
	}
	if($ruf[2] =~ /-(\d+)$/) {
	    $end_forward = $1;
	} else {
	    $end_forward = "*";
	}
	if($rur[2] =~ /^(\d+)-/) {
	    $start_reverse = $1;
	} else {
	    $start_reverse = "*";
	}
	if($rur[2] =~ /-(\d+)$/) {
	    $end_reverse = $1;
	} else {
	    $end_reverse = "*";
	}
	if($rum_u_forward =~ /\S/ && !($rum_u_reverse =~ /\S/)) {
	    $start_reverse = $start_forward;
	    $end_reverse = $start_forward;
	}
	if($rum_u_reverse =~ /\S/ && !($rum_u_forward =~ /\S/)) {
	    $start_forward = $start_reverse;
	    $end_forward = $start_reverse;
	}
	if($rum_u_forward =~ /\S/ && $rum_u_reverse =~ /\S/) {
	    if($ruf[3] eq "+") {
		$idist_f = $end_reverse - $start_forward;
	    } else {
		$idist_f = $end_forward - $start_reverse;
	    }
	    $idist_r = -1 * $idist_f;
	}


# PRINTING OUT SAM RECORD STARTS HERE

# FORWARD:

	$forward_record = "";
	$forward_record = $forward_record . "seq.$seqnum";
	$forward_record = $forward_record . "a\t$bitscore_f";

	if(!($rum_u_forward =~ /\S/) && $rum_u_reverse =~ /\S/) { # forward unmapped, reverse mapped
	    $forward_record = $forward_record . "\t$rur[1]\t$start_reverse\t255\t*\t=\t$start_reverse\t0\t$forward_read\t$forward_qual";
	}
	if($rum_u_forward =~ /\S/ || $rum_u_joined =~ /\S/) { # forward mapped
	    $forward_record = $forward_record . "\t$ruf[1]\t$start_forward\t255\t$CIGAR_f\t";
	    if($paired eq "true") {
		if($rum_u_reverse =~ /\S/) { # paired and reverse mapped
		    $forward_record = $forward_record . "=\t$start_reverse\t$idist_f\t$forward_read\t$forward_qual";
		} else { # reverse didn't map
		    $forward_record = $forward_record . "=\t$start_forward\t0\t$forward_read\t$forward_qual";
		}
	    } else { # not paired end
		$forward_record = $forward_record . "*\t0\t0\t$forward_read\t$forward_qual";
	    }
	}
	if($joined eq "true") {
	    $forward_record = $forward_record . "\tOL:A:T";
	} else {
	    $forward_record = $forward_record . "\tOL:A:F";
	}
	$forward_record = $forward_record . "\n";
#	print SAM $forward_record;
	print $forward_record;

# REVERSE

	if($paired eq "true") {
	    $reverse_record = "";
	    $reverse_record = $reverse_record . "seq.$seqnum";
	    $reverse_record = $reverse_record . "b\t$bitscore_r";
	    if(!($rum_u_reverse =~ /\S/) && $rum_u_forward =~ /\S/) {  # reverse unmapped, forward mapped
		$reverse_record = $reverse_record . "\t$ruf[1]\t$start_forward\t255\t*\t=\t$start_forward\t0\t$reverse_read\t$reverse_qual";
	    }
	    if($rum_u_reverse =~ /\S/ || $rum_u_joined =~ /\S/) { # reverse mapped
		$reverse_record = $reverse_record . "\t$rur[1]\t$start_reverse\t255\t$CIGAR_r\t=";
		if($rum_u_forward =~ /\S/) { # forward mapped
		    $reverse_record = $reverse_record . "\t$start_forward\t$idist_r\t$reverse_read\t$reverse_qual";
		} else { # forward didn't map
		    $reverse_record = $reverse_record . "\t$start_reverse\t0\t$reverse_read\t$reverse_qual";		
		}
	    }
	    if($joined eq "true") {
		$reverse_record = $reverse_record . "\tOL:A:T";
	    } else {
		$reverse_record = $reverse_record . "\tOL:A:F";
	    }
	    $reverse_record = $reverse_record . "\n";
#	    print SAM $reverse_record;
	    print $reverse_record;
	}
    }

    if($unique_mapper_found eq "false" && $non_unique_mappers_found eq "false") {
	# neither forward nor reverse map
	if($paired eq "false") {
	    $record = "seq.$seqnum";
	    $record = $record . "a\t4\t*\t0\t255\t*\t*\t0\t0\t$forward_read\t$forward_qual\n";
#	    print SAM $record;
	    print $record;
	} else {
	    $record = "seq.$seqnum";
	    $record = $record . "a\t77\t*\t0\t255\t*\t*\t0\t0\t$forward_read\t$forward_qual\n";
	    $record = $record . "seq.$seqnum";
	    $record = $record . "b\t141\t*\t0\t255\t*\t*\t0\t0\t$reverse_read\t$reverse_qual\n";
#	    print SAM $record;
	    print $record;
	}
    }
}

sub getsuffix () {
    ($spans, $suffixlength) = @_;

    $prefixlength = &spansTotalLength($spans) - $suffixlength;
    $newspans = "";
    @OS = split(/, /, $spans);
    $running_length=0;
    for($os=0; $os<@OS; $os++) {
	@B = split(/-/, $OS[$os]);
	$running_length = $running_length + $B[1] - $B[0] + 1;
	if($running_length > $prefixlength) {
	    $STRT = $B[1] - ($running_length - $prefixlength) + 1;
	    $newspans = $STRT . "-" . $B[1];
	    $BB = $B[1];
	    $spans = $spans . ", ";
	    $spans =~ s/^.*-$BB, //;
	    if($spans =~ /\S/) {
		$newspans = $newspans . ", " . $spans;
	    }
	    $newspans =~ s/^\s*,\s*//;
	    $newspans =~ s/\s*,\s*$//;
	    return $newspans;
	}
    }
}

sub getprefix () {
    ($spans, $prefixlength) = @_;

    $newspans = "";
    @OS = split(/, /, $spans);
    $running_length=0;
    for($os=0; $os<@OS; $os++) {
	@B = split(/-/, $OS[$os]);
	$running_length = $running_length + $B[1] - $B[0] + 1;
	if($running_length >= $prefixlength) {
	    $END = $B[1] - ($running_length - $prefixlength);
	    if($newspans =~ /\S/) {
		$newspans =  $newspans . ", " . $B[0] . "-" . $END;
	    } else {
		$newspans = $B[0] . "-" . $END;
	    }
	    $newspans =~ s/^\s*,\s*//;
	    $newspans =~ s/\s*,\s*$//;
	    return $newspans;
	} else {
	    if($newspans =~ /\S/) {
		$newspans = $newspans . ", " . $B[0] . "-" . $B[1];
	    } else {
		$newspans = $B[0] . "-" . $B[1];
	    }
	}
    }
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

sub reversecomplement () {
    ($sq) = @_;
    @A = split(//,$sq);
    $rev = "";
    for($i=@A-1; $i>=0; $i--) {
	$flag = 0;
	if($A[$i] eq 'A') {
	    $rev = $rev . "T";
	    $flag = 1;
	}
	if($A[$i] eq 'T') {
	    $rev = $rev . "A";
	    $flag = 1;
	}
	if($A[$i] eq 'C') {
	    $rev = $rev . "G";
	    $flag = 1;
	}
	if($A[$i] eq 'G') {
	    $rev = $rev . "C";
	    $flag = 1;
	}
	if($flag == 0) {
	    $rev = $rev . $A[$i];
	}
    }
    return $rev;
}