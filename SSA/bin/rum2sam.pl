#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

$|=1;

if(@ARGV < 5) {
    die "
Usage: rum2sam.pl <rum unique file> <rum nu file> <reads file> <quals file> <sam outfile> [options]

<reads file> and <qual file> are the files coming from parse2fasta.pl and fastq2qualities.pl

";
}

$rum_unique_file = $ARGV[0];
$rum_nu_file = $ARGV[1];
$reads_file = $ARGV[2];
$qual_file = $ARGV[3];
$sam_outfile = $ARGV[4];

open(INFILE, $reads_file);
$line = <INFILE>;
chomp($line);
$line =~ /seq.(\d+)/;
$firstseqnum = $1;
$line = <INFILE>;
chomp($line);
$readlength = length($line);
print "readlength = $readlength\n";
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

open(RUMU, $rum_unique_file);
open(RUMNU, $rum_nu_file);
open(READS, $reads_file);
open(QUALS, $qual_file);
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
    $forward_qual = <QUALS>;
    $forward_qual = <QUALS>;
    chomp($forward_qual);
    if($paired eq "true") {
	$reverse_qual = <QUALS>;
	$reverse_qual = <QUALS>;
	chomp($reverse_qual);
    }
#    print "------\n$seqnum\n";
#    print "forward_read = $forward_read\n";
#    print "forward_qual = $forward_qual\n";
#    print "reverse_read = $reverse_read\n";
#    print "reverse_qual = $reverse_qual\n";

    $flag = 0;
    $rum_u_forward = "";
    $rum_u_reverse = "";
    $rum_u_joined = "";
    $unique_mapper_found = "false";
    while($flag == 0) {
	$line = <RUMU>;
	chomp($line);
	$line =~ /seq.(\d+)(.)/;
	$sn = $1;
	$type = $2;
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
    }
    if($unique_mapper_found eq "true") {
	if($rum_u_forward =~ /\S/) {
#	    print "rum_u_forward = $rum_u_forward\n";
	    @ruf = split(/\t/,$rum_u_forward);
	    $ruf[4] =~ s/://g;
	    $ruf[4] =~ s/\+//g;
	    $rum_u_forward_length = length($ruf[4]);
	    if($ruf[3] eq "-") {
		$forward_read = reversecomplement($forward_read);
	    }
#	    print "$forward_read\n";
	    if($rum_u_forward_length < $readlength) {
		$prefix_offset_forward = 0;
		$x = $forward_read;
		$y = $ruf[4];
		until($x =~ /^$y/) {
		    $x =~ s/^.//;
		    $prefix_offset_forward++;
#		    print " ";
		}
	    }
#	    print "$ruf[4]\n";
#	    print "rum_u_forward_length = $rum_u_forward_length\n";
	}
	if($rum_u_reverse =~ /\S/) {
#	    print "rum_u_reverse = $rum_u_reverse\n";
	    @rur = split(/\t/,$rum_u_reverse);
	    $rur[4] =~ s/://g;
	    $rur[4] =~ s/\+//g;
	    $rum_u_reverse_length = length($rur[4]);
	    if($rur[3] eq "+") {
		$reverse_read = reversecomplement($reverse_read);
	    }
#	    print "$reverse_read\n";
	    if($rum_u_reverse_length < $readlength) {
		$prefix_offset_reverse = 0;
		$x = $reverse_read;
		$y = $rur[4];
		until($x =~ /^$y/) {
		    $x =~ s/^.//;
		    $prefix_offset_reverse++;
#		    print " ";
		}
	    }
#	    print "$rur[4]\n";
#	    print "rum_u_reverse_length = $rum_u_reverse_length\n";
	}
	if($rum_u_joined =~ /\S/) {
	    print "--------\n$seqnum\n";
	    print "rum_u_joined = $rum_u_joined\n";
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
	    until($y =~ /^$x/) {
		$x =~ s/^.//;
		$count++;
		$prefix_offset_upstream++;
		if($L - $count < 15) {
		    $x = $upstream_read;
		    $suffix_offset_upstream++;
		    for($j=0; $j<$suffix_offset_upstream; $j++) {
			$x =~ s/.$//;
		    }
		    $prefix_offset_upstream = 0;
		    $count = 0;
		    $L = length($x);
		}
	    }	    
	    $UR = $upstream_read;
	    $replace = "";
	    for($i=0; $i<$suffix_offset_upstream; $i++) {
		$UR =~ s/.$//;
		$replace = $replace . "X";
	    }
	    $UR = $UR . $replace;


	    if($ruj[3] eq "-") {
		$downstream_read = reversecomplement($forward_read);
	    } else {
		$downstream_read = reversecomplement($reverse_read);
	    }
	    $x = $downstream_read;
	    $y = $ruj[4];
	    $suffix_offset_downstream = 0;
	    $L = length($x);
	    $count=0;
	    $prefix_offset_downstream = 0;
	    until($y =~ /$x$/) {
#		print "x=$x\n";
#		print "y=$y\n";
		$x =~ s/.$//;
		$count++;
		$suffix_offset_downstream++;
		if($L - $count < 15) {
		    $x = $downstream_read;
		    $prefix_offset_downstream++;
		    for($j=0; $j<$prefix_offset_downstream; $j++) {
			$x =~ s/^.//;
		    }
		    $suffix_offset_downstream = 0;
		    $count = 0;
		    $L = length($x);
		}
	    }	    
	    $DR = $downstream_read;
	    $replace = "";
	    for($i=0; $i<$prefix_offset_downstream; $i++) {
		$DR =~ s/^.//;
		$replace = $replace . "X";
	    }
	    $DR = $replace . $DR;
	    $offset = length($ruj[4]) + $prefix_offset_upstream + $suffix_offset_downstream - length($DR);
	    print "\n$upstream_read\n\n";
	    $OFFSET = $prefix_offset_upstream + $readlength - length($ruj[4]) - $suffix_offset_downstream;
	    $P = "";
	    for($i=0; $i<$OFFSET; $i++) {
		$P = $P . " ";
	    }
	    print $P;
	    print "$UR\n";
	    print $P;
	    for($i=0; $i<$prefix_offset_upstream; $i++) {
		print " ";
	    }
	    print "$ruj[4]\n";
	    for($i=0; $i<$offset; $i++) {
		print " ";
	    }
	    print "$DR\n";
	    print "\n$downstream_read\n";
	}
    }
}

while($flag == 1) {
    @a = split(/\t/,$line1);
    if($line1 =~ /seq.(\d+)a/) {
	$seqnum_a = $1;
	$line2 = <INFILE>;
	chomp($line2);
	@b = split(/\t/,$line2);	
	if($line2 =~ /seq.(\d+)b/) {
	    $seqnum_q = $1;
	    if($seqnum_a == $seqnum_q) { # consistent pair found
		$seqnum_b = $seqnum_q;
		$seq_a = $a[4];
		$seq_a =~ s/://g;
		$seq_a =~ s/\+//g;
		$seq_b = $b[4];
		$seq_b =~ s/\+//g;
		$a_insertion = "false";
		if($a[4] =~ /\+/) {
		    $a_insertion = "true";
		    if($a[4] =~ /^([^+]*)\+([^+]*)\+/) {
			$a_prefix = $1;
			$a_insert = $2;
			$a_prefix_length = length($a_prefix);
			$a_insert_length = length($a_insert);
		    }			
		}
		$b_insertion = "false";
		if($b[4] =~ /\+/) {
		    $b_insertion = "true";
		    if($b[4] =~ /^([^+]*)\+([^+]*)\+/) {
			$b_prefix = $1;
			$b_insert = $2;
			$b_prefix_length = length($b_prefix);
			$b_insert_length = length($b_insert);
		    }			
		}
		$alignment_length_a = length($seq_a);
		$alignment_length_b = length($seq_b);
		if($alignment_length_a == $readlength && $alignment_length_b == $readlength) {
		    $strand = $a[3];
		    $BITFLAGa = 67;
		    $BITFLAGb =131;
		    if($strand eq '-') {
			$BITFLAGa = $BITFLAGa + 16;
			$BITFLAGb = $BITFLAGb + 32;
		    } else {
			$BITFLAGa = $BITFLAGa + 32;
			$BITFLAGb = $BITFLAGb + 16;
		    }
		    $a[2] =~ /^(\d+)/;
		    $forward_start = $1;
		    $b[2] =~ /^(\d+)/;
		    $reverse_start = $1;
		    print "$a[0]\t$BITFLAGa\t$a[1]\t$forward_start\t255\t";
		    @aspans = split(/, /,$a[2]);
		    @C1 = split(/-/,$aspans[0]);
		    $L = $C1[1] - $C1[0] + 1;
		    $CIGAR = $L . "M";
		    $running_length = $L;
		    for($i=1; $i<@aspans; $i++) {
			@C2 = split(/-/,$aspans[$i]);
			$skipped = $C2[0] - $C1[1] - 1;
			if($skipped >= 15) {
			    $CIGAR = $CIGAR . $skipped . "N";
			} else {
			    $CIGAR = $CIGAR . $skipped . "D";
			}
			$L = $C2[1] - $C2[0] + 1;
			$CIGAR = $CIGAR . $L . "M";
			$running_length = $running_length + $L;
			$C1[0] = $C2[0];
			$C1[1] = $C2[1];
		    }
		    if($strand eq "+") {
			$b[2] =~ /-(\d+)$/;
			$reverse_end = $1;
			$ilen = $reverse_end - $forward_start;
		    } else {
			$a[2] =~ /-(\d+)$/;
			$forward_end = $1;
			$ilen = $reverse_start - $forward_end;
		    }
		    print $CIGAR . "\t=\t$reverse_start\t$ilen\t$seq_a\t";
		    for($i=0; $i<$alignment_length_a; $i++) {
			print "h";
		    }
		    print "\n";

		    print "$b[0]\t$BITFLAGb\t$a[1]\t$reverse_start\t255\t";
		    @bspans = split(/, /,$b[2]);
		    @C1 = split(/-/,$bspans[0]);
		    $L = $C1[1] - $C1[0] + 1;
		    $CIGAR = $L . "M";
		    for($i=1; $i<@bspans; $i++) {
			@C2 = split(/-/,$bspans[$i]);
			$skipped = $C2[0] - $C1[1] - 1;
			if($skipped >= 15) {
			    $CIGAR = $CIGAR . $skipped . "N";
			} else {
			    $CIGAR = $CIGAR . $skipped . "D";
			}
			$L = $C2[1] - $C2[0] + 1;
			$CIGAR = $CIGAR . $L . "M";
			$C1[0] = $C2[0];
			$C1[1] = $C2[1];
		    }
		    $ilen = -1 * $ilen;
		    print $CIGAR . "\t=\t$reverse_start\t$ilen\t$seq_b\t";
		    for($i=0; $i<$alignment_length_b; $i++) {
			print "h";
		    }
		    print "\n";
		}
		$line1 = <INFILE>;
		chomp($line1);
	    } else { # only forward exists
		$line1 = $line2;
	    }
	}
    } else { # only reverse exists
	$line1 = <INFILE>;
	chomp($line1);
    }
    if($line1 eq '') {
	$flag = 0;
    }
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


# seq.1a  99      chrM:1-16299_strand=+   5920    46      100M    =       6068    248     CACTACCAGTGCTAGCCGCAGGCATTACTATACTACTAACAGACCGCAACCTAAACACAACTTTCTTTGATCCCGCTGGAGGAGGGGACCCAATTCTCTA    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa    XT:A:U  NM:i:0  SM:i:23 AM:i:23 X0:i:1  X1:i:1  XM:i:0  XO:i:0  XG:i:0  MD:Z:100        XA:Z:chr2:1-181748087_strand=+,-22444476,100M,1;

# seq.1b  147     chrM:1-16299_strand=+   6068    46      100M    =       5920    -248    CCTCCCAGGATTTGGAATTATTTCACATGGAGTTACTTACTACTCCGGAAAAAAAGAACCTTTCGGCTATATAGGAATAGTATGAGCAATAATGTCTATT    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa    XT:A:U  NM:i:1  SM:i:23 AM:i:23 X0:i:1  X1:i:1  XM:i:1  XO:i:0  XG:i:0  MD:Z:29T70      XA:Z:chr2:1-181748087_strand=+,+22444328,100M,2;


#seq.5a  chr7    52636876-52636975       CAAAATCCCTGCCTTACCTGAGCCGACTGACCCTGCCACCACCTAGCCTGGCCAGCTGTGTCCTCCAAGCCGAGGTCAATGGCCAGCAGGGCATGCAGCC
#seq.5b  chr7    52636995-52637094       AGTTGAGGAGCATAAGACATAGTGGGAAGGCCCAGAGCCCAGGCCTGTGATCAGGCAGGGCTGTTAAGAAGGCCAAAAAGTCCCTGTCCAGGTCCCAGCC
