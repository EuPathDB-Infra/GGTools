#/usr/bin/perl

if(@ARGV != 8) {
    print STDERR 
"
Usage: make_GU_and_TU.pl -infile NAME -gufile NAME -tufile NAME -type <single|paired>

  Where: -infile is the file output from sort_bowtie.pl
         -gufile is the file to be output containing unique genome alignments
         -gnufile is the file to be output containing non-nique genome alignments
         -type is 'single' for single-end reads and 'paired' for paired-end reads

";
    exit(1);
}

for($i=0; $i<@ARGV; $i++) {
    $optionrecognized = 1;
    if($ARGV[$i] eq "-infile") {
	$i++;
	$infile = $ARGV[$i];
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-gufile") {
	$i++;
	$outfile1 = $ARGV[$i];
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-gnufile") {
	$i++;
	$outfile2 = $ARGV[$i];
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-type") {
	if($ARGV[$i+1] eq "single") {
	    $paired_end = "false";
	    $optionrecognzied = 0;
	}
	if($ARGV[$i+1] eq "paired") {
	    $paired_end = "true";
	    $optionrecognzied = 0;
	}
    }

    if($optionrecognized == 0) {
	die "\nERROR: option $ARGV[$i+1] not recognized\n";
    }
}

$|=1;

$max_fragment_length = 500000;
open(INFILE, $infile) or die "\nERROR: Cannot open -infile '$infile'\n";
$t = `tail -1 $infile`;
$t =~ /seq.(\d+)/;
$num_seqs = $1;
$line = <INFILE>;
chomp($line);
open(OUTFILE1, ">$outfile1") or die "\nERROR: Cannot open -gufile '$outfile1' for writing\n";
open(OUTFILE2, ">$outfile2") or die "\nERROR: Cannot open -tufile '$outfile2' for writing\n";
print "num_seqs = $num_seqs\n";
for($seqnum=1; $seqnum<$num_seqs; $seqnum++) {
    $numa=0;
    $numb=0;
    undef %a_reads;
    undef %b_reads;
    while($line =~ /seq.($seqnum)a/) {
	$seqs_a[$numa] = $line;
	$numa++;
	$line = <INFILE>;
	chomp($line);
    }
    while($line =~ /seq.($seqnum)b/) {
	$seqs_b[$numb] = $line;
	$numb++;
	$line = <INFILE>;
	chomp($line);
    }
    if($numa > 0 || $numb > 0) {
	$num_different_a = 0;
	for($i=0; $i<$numa; $i++) {
	    $line2 = $seqs_a[$i];
	    if(!($line2 =~ /^N+$/)) {
		@a = split(/\t/,$line2);
		$id = $a[0];
		$strand = $a[1];
		$chr = $a[2];
		$chr =~ s/:.*//;
		$start = $a[3]+1;
		$seq = $a[4];
		if($seq =~ /^(N+)/) {
		    $seq =~ s/^(N+)//;
		    $Nprefix = $1;
		    @x = split(//,$Nprefix);
		    $start = $start + @x;
		}
		$seq =~ s/N+$//;
		@x = split(//,$seq);
		$seqlength = @x;
		$end = $start + $seqlength - 1; 
	    }
	    $a_reads{"$id\t$strand\t$chr\t$start\t$end\t$seq"}++;
	    if($a_reads{"$id\t$strand\t$chr\t$start\t$end\t$seq"} == 1) {
		$num_different_a++;
	    }
	}
	$num_different_b = 0;
	for($i=0; $i<$numb; $i++) {
	    $line2 = $seqs_b[$i];
	    if(!($line2 =~ /^N+$/)) {
		@a = split(/\t/,$line2);
		$id = $a[0];
		$strand = $a[1];
		$chr = $a[2];
		$chr =~ s/:.*//;
		$start = $a[3]+1;
		$seq = $a[4];
		if($seq =~ /^(N+)/) {
		    $seq =~ s/^(N+)//;
		    $Nprefix = $1;
		    @x = split(//,$Nprefix);
		    $start = $start + @x;
		}
		$seq =~ s/N+$//;
		@x = split(//,$seq);
		$seqlength = @x;
		$end = $start + $seqlength - 1; 
	    }
	    $b_reads{"$id\t$strand\t$chr\t$start\t$end\t$seq"}++;
	    if($b_reads{"$id\t$strand\t$chr\t$start\t$end\t$seq"} == 1) {
		$num_different_b++;
	    }
	}
    }
# NOTE: the following three if's cover all cases we care about, because if numa > 1 and numb = 0, then that's
# not really ambiguous, blat might resolve it

    if($num_different_a == 1 && $num_different_b == 0) { # unique forward match, no reverse
	foreach $key (keys %a_reads) {
	    $key =~ s/\t\+//;
	    $key =~ s/\t-//;
	    $key =~ /^[^\t]+\t[^\t]+\t([^\t]+\t[^\t]+)\t/;
	    $xx = $1;
	    $yy = $xx;
	    $xx =~ s/\t/-/;
	    $key =~ s/$yy/$xx/;
	    print OUTFILE1 "$key\n";
	}
    }
    if($num_different_a == 0 && $num_different_b == 1) { # unique reverse match, no forward
	foreach $key (keys %b_reads) {
	    $key =~ s/\t\+//;
	    $key =~ s/\t-//;
	    $key =~ /^[^\t]+\t[^\t]+\t([^\t]+\t[^\t]+)\t/;
	    $xx = $1;
	    $yy = $xx;
	    $xx =~ s/\t/-/;
	    $key =~ s/$yy/$xx/;
	    print OUTFILE1 "$key\n";
	}
    }
    if($paired_end eq "false") {
	if($num_different_a > 1) { 
	    foreach $key (keys %a_reads) {
		$key =~ s/\t\+//;
		$key =~ s/\t-//;
		$key =~ /^[^\t]+\t[^\t]+\t([^\t]+\t[^\t]+)\t/;
		$xx = $1;
		$yy = $xx;
		$xx =~ s/\t/-/;
		$key =~ s/$yy/$xx/;
		print OUTFILE2 "$key\n";
	    }
	}
    }
    if(($num_different_a > 0 && $num_different_b > 0) && ($num_different_a * $num_different_b < 1000000)) { 
# forward and reverse matches, must check for consistency, but not if more than 1,000,000 possibilities,
# in that case skip...
	undef %consistent_mappers;
	foreach $akey (keys %a_reads) {
	    foreach $bkey (keys %b_reads) {

		@a = split(/\t/,$akey);
		$aid = $a[0];
		$astrand = $a[1];
		$achr = $a[2];
		$astart = $a[3];
		$aend = $a[4];
		$aseq = $a[5];
		@a = split(/\t/,$bkey);
		$bstrand = $a[1];
		$bchr = $a[2];
		$bstart = $a[3];
		$bend = $a[4];
		$bseq = $a[5];
		if($astrand eq "+" && $bstrand eq "-") {
		    if($achr eq $bchr && $astart <= $bstart && $bstart - $astart < $max_fragment_length) {
			if($bstart > $aend + 1) {
			    $akey =~ s/\t\+//;
			    $akey =~ s/\t-//;
			    $bkey =~ s/\t\+//;
			    $bkey =~ s/\t-//;
			    $akey =~ /^[^\t]+\t[^\t]+\t([^\t]+\t[^\t]+)\t/;
			    $xx = $1;
			    $yy = $xx;
			    $xx =~ s/\t/-/;
			    $akey =~ s/$yy/$xx/;
			    $bkey =~ /^[^\t]+\t[^\t]+\t([^\t]+\t[^\t]+)\t/;
			    $xx = $1;
			    $yy = $xx;
			    $xx =~ s/\t/-/;
			    $bkey =~ s/$yy/$xx/;
			    $consistent_mappers{"$akey\n$bkey\n"}++;
			}
			else {
			    $overlap = $aend - $bstart + 1;
			    @sq = split(//,$bseq);
			    $joined_seq = $aseq;
			    for($i=$overlap; $i<@sq; $i++) {
				$joined_seq = $joined_seq . $sq[$i];
			    }
			    $aid =~ s/a//;
			    $consistent_mappers{"$aid\t$achr\t$astart-$bend\t$joined_seq\n"}++;
			}
		    }
		}
		if($astrand eq "-" && $bstrand eq "+") {
		    if($achr eq $bchr && $bstart <= $astart && $astart - $bstart < $max_fragment_length) {
			if($astart > $bend + 1) {
			    $akey =~ s/\t\+//;
			    $akey =~ s/\t-//;
			    $bkey =~ s/\t\+//;
			    $bkey =~ s/\t-//;
			    $akey =~ /^[^\t]+\t[^\t]+\t([^\t]+\t[^\t]+)\t/;
			    $xx = $1;
			    $yy = $xx;
			    $xx =~ s/\t/-/;
			    $akey =~ s/$yy/$xx/;
			    $bkey =~ /^[^\t]+\t[^\t]+\t([^\t]+\t[^\t]+)\t/;
			    $xx = $1;
			    $yy = $xx;
			    $xx =~ s/\t/-/;
			    $bkey =~ s/$yy/$xx/;
			    $consistent_mappers{"$akey\n$bkey\n"}++;
			}
			else {
			    $overlap = $bend - $astart + 1;
			    @sq = split(//,$bseq);
			    $joined_seq = "";
			    for($i=0; $i<@sq-$overlap; $i++) {
				$joined_seq = $joined_seq . $sq[$i];
			    }
			    $joined_seq = $joined_seq . $aseq;
			    $aid =~ s/a//;
			    $consistent_mappers{"$aid\t$achr\t$bstart-$aend\t$joined_seq\n"}++;
			}
		    }
		}
	    }
	}
	$count = 0;
	foreach $key (keys %consistent_mappers) {
	    $count++;
	    $str = $key;
	}
	if($count == 1) {
	    print OUTFILE1 $str;
	}
	if($count > 1) {
# add something here so that if all consistent mappers agree on some
# exons, then those exons will still get reported, each on its own line
	    foreach $key (keys %consistent_mappers) {
		print OUTFILE2 $key;
	    }
	}
    }
}
close(OUTFILE1);
close(OUTFILE2);
