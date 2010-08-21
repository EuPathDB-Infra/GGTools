#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

$|=1;

if(@ARGV < 7) {
    die "

Usage: make_RUM_junctions_file.pl <rum_unique> <rum_nu> <genome seq> <gene annotations> <all junctions outfile rum-format> <all junctions outfile bed-format> <high quality junctions outfile bed-format> [options]

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
$genome_sequence = $ARGV[2];
$gene_annot = $ARGV[3];
$outfile1 = $ARGV[4];
$outfile2 = $ARGV[5];
$outfile3 = $ARGV[6];

open(OUTFILE1, ">$outfile1") or die "\nError: cannot open file '$outfile1' for writing\n\n";
print OUTFILE1 "intron\tscore\tknown\tcanonical_splice_signals\t\tambiguous\tlong_overlap_unique_reads\tshort_overlap_unique_reads\tlong_overlap_nu_reads\tshort_overlap_nu_reads\n";

open(OUTFILE2, ">$outfile2") or die "\nError: cannot open file '$outfile2' for writing\n\n";
print OUTFILE2 "track\tname=rum_junctions_all\tdescription=\"RUM junctions (all)\" itemRgb=\"On\"\n";

open(OUTFILE3, ">$outfile3") or die "\nError: cannot open file '$outfile3' for writing\n\n";
print OUTFILE3 "track\tname=rum_junctions_hq\tdescription=\"RUM high quality junctions\" itemRgb=\"On\"\n";

# read in known junctions to color them green in the hq track:

open(INFILE, $gene_annot) or die "\nError: cannot open file '$gene_annot' for reading\n\n";
while($line = <INFILE>) {
    @a = split(/\t/, $line);
    $chr = $a[0];
    $a[5] =~ s/\s*,\s*$//;
    $a[6] =~ s/\s*,\s*$//;
    $a[5] =~ s/^\s*,\s*//;
    $a[6] =~ s/^\s*,\s*//;
    @starts = split(/\s*,\s*/,$a[5]);
    @ends = split(/\s*,\s*/,$a[6]);
    for($i=0; $i<@starts-1; $i++) {
	$S = $ends[$i] + 1;
	$E = $starts[$i+1];
	$intron = $chr . ":" . $S . "-" . $E;
	$knownintron{$intron} = 1;
    }
}
close(INFILE);

$faok = "false";
$minintron = 15;

for($i=7; $i<@ARGV; $i++) {
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
    `perl modify_fa_to_have_seq_on_one_line.pl $genome_sequence > $f`;
    open(GENOMESEQ, $f) or die "\nError: cannot open file '$f' for reading\n\n";
} else {
    open(GENOMESEQ, $genome_sequence) or die "\nError: cannot open file '$genome_sequence' for reading\n\n";
}

$FLAG = 0;
while($FLAG == 0) {
    undef %CHR2SEQ;
    undef %allintrons;
    undef @amb;
    undef @badoverlapU;
    undef @goodoverlapU;
    undef @badoverlapNU;
    undef @goodoverlapNU;

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
	    print STDERR "chr=$chr\n";
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
    &printjunctions();
}
close(GENOMESEQ);

sub printjunctions () {
    foreach $intron (keys %allintrons) {
	$amb{$intron} = $amb{$intron} + 0;
	$badoverlapU{$intron} = $badoverlapU{$intron} + 0;
	$goodoverlapU{$intron} = $goodoverlapU{$intron} + 0;
	$badoverlapNU{$intron} = $badoverlapNU{$intron} + 0;
	$goodoverlapNU{$intron} = $goodoverlapNU{$intron} + 0;
	$knownintron{$intron} = $knownintron{$intron} + 0;

# chromosome
# start seg 1: 50 bases upstream from junction start
# start seg 2: 50 bases upstream from junction end
# score: goodoverlap_badoverlap
# 50
# +
# start seg 1: 50 bases upstream from junction start (again)
# start seg 2: 50 bases upstream from junction end (again)
# 
# Color:
#    0,0,128 NAVY (for high quality)
#    255,69,0 RED (for low quality)
# 2
# 50,50
# 0, intron_length + 50
 
	$intron =~ /^(.*):(\d+)-(\d+)$/;
	$chr = $1;
	$start = $2 - 1;
	$end = $3;
	$end2 = $end + 50;
	$start2 = $start - 50;
	$ilen = $end - $start + 50;
	$LEN1 = 50;
	$LEN2 = 50;
	if($start2 < 0) {
	    $adjust = $start2;
	    $start2 = 0;
	    $LEN1 = $LEN1 + $adjust;
	    $ilen = $ilen + $adjust;
	}
	if($goodoverlapU{$intron} > 0 && $goodsplicesignal{$intron} == 1) {
	    $N = $goodoverlapU{$intron} + $goodsplicesignal{$intron};
	    print OUTFILE1 "$intron\t$N\t$knownintron{$intron}\t$goodsplicesignal{$intron}\t$amb{$intron}\t$goodoverlapU{$intron}\t$badoverlapU{$intron}\t$goodoverlapNU{$intron}\t$badoverlapNU{$intron}\n";
	    print OUTFILE2 "$chr\t$start2\t$end2\t$N\t50\t+\t$start2\t$end2\t0,0,128\t2\t$LEN1,$LEN2\t0,$ilen\n";
	    if($knownintron{$intron}==1) {
		print OUTFILE3 "$chr\t$start2\t$end2\t$N\t50\t+\t$start2\t$end2\t0,0,128\t2\t$LEN1,$LEN2\t0,$ilen\n";
	    } else {
		print OUTFILE3 "$chr\t$start2\t$end2\t$N\t$M\t50\t+\t$start2\t$end2\t0,205,102\t2\t$LEN1,$LEN2\t0,$ilen\n";
	    }
	} else {
	    print OUTFILE1 "$intron\t0\t$knownintron{$intron}\t$goodsplicesignal{$intron}\t$amb{$intron}\t$goodoverlapU{$intron}\t$badoverlapU{$intron}\t$goodoverlapNU{$intron}\t$badoverlapNU{$intron}\n";
	    print OUTFILE2 "$chr\t$start2\t$end2\t0\t$M\t50\t+\t$start2\t$end2\t255,69,0\t2\t$LEN1,$LEN2\t0,$ilen\n";
	}
    }
}

sub getjunctions () {
    open(INFILE, $rumU) or die "\nError: cannot open file '$rumU' for reading\n\n";
    while($line = <INFILE>) {
	if(!($line =~ /, /)) {
	    next;
	}
	chomp($line);
	@a = split(/\t/,$line);
	$chr = $a[1];
	if(!(defined $CHR2SEQ{$chr})) {
	    next;
	}
	$seq = $a[4];
	$snt = $a[0];
	$snt =~ s/seq.//;
#	print STDERR "1:seq.$snt\n";
	while($seq =~ /^([^+]*)\+/) {  # removing the insertions
	    $pref = $1;
	    $seq =~ s/^$pref\+[^+]+\+/$pref/;
	}
	$strand = $a[3];
	@SPANS = split(/, /,$a[2]);
	@SEQ = split(/:/, $seq);
	for($i=0; $i<@SPANS-1; $i++) {
	    @c1 = split(/-/,$SPANS[$i]);
	    @c2 = split(/-/,$SPANS[$i+1]);
	    $elen1 = $c1[1] - $c1[0] + 1;
	    $elen2 = $c2[1] - $c2[0] + 1;
	    $ilen = $c2[0] - $c1[1] - 1;
	    $istart = $c1[1]+1;
	    $iend = $c2[0]-1;
	    $intron = $chr . ":" . $istart . "-" . $iend;
	    $altintron = "";
	    if($ilen >= $minintron) {
		$allintrons{$intron} = 1;
		if(!(defined $amb{$intron}) || !($goodsplicesignal{$intron})) {
		    $SEQ[$i] =~ /(.)$/;
		    $leftexon_lastbase = $1;
		    $SEQ[$i+1] =~ /^(.)/;
		    $rightexon_firstbase = $1;
		    $intron_firstbase = substr($CHR2SEQ{$chr}, $istart-1, 1);
		    $intron_lastbase = substr($CHR2SEQ{$chr}, $iend-1, 1);
		    $splice_signal_upstream = substr($CHR2SEQ{$chr}, $istart-1, 2);
		    $splice_signal_downstream = substr($CHR2SEQ{$chr}, $iend-2, 2);
		    if(($splice_signal_upstream eq "GT" && $splice_signal_downstream eq "AG") || ($splice_signal_upstream eq "CT" && $splice_signal_downstream eq "AC")) {
			$goodsplicesignal{$intron} = 1;
		    } else {
			$goodsplicesignal{$intron} = 0;
		    }
		    if($leftexon_lastbase eq $intron_lastbase) {
			$istart_alt = $istart-1;
			$iend_alt = $iend-1;
			$altintron = $chr . ":" . $istart_alt . "-" . $iend_alt;
			$amb{$intron}=1;  # amb for ambiguous
			$amb{$altintron}=1;
			$allintrons{$altintron} = 1;
		    }
		    if($rightexon_firstbase eq $intron_firstbase) {
			$istart_alt = $istart+1;
			$iend_alt = $iend+1;
			$altintron = $chr . ":" . $istart_alt . "-" . $iend_alt;
			$amb{$intron}=1;  # amb for ambiguous
			$amb{$altintron}=1;
			$allintrons{$altintron} = 1;
		    }
		}
		if($elen1 <= 7 || $elen2 <= 7) {
		    $badoverlapU{$intron}++;
		    if($altintron =~ /\S/) {
			$badoverlapU{$altintron}++;			    
		    }
		} else {
		    $goodoverlapU{$intron}++;
		    if($altintron =~ /\S/) {
			$goodoverlapU{$altintron}++;			    
		    }
		}
	    }
	}
    }
    close(INFILE);
    print STDERR "finished Unique\n";
    open(INFILE, $rumNU) or die "\nError: cannot open file '$rumNU' for reading\n\n";
    while($line = <INFILE>) {
	if(!($line =~ /, /)) {
	    next;
	}
	chomp($line);
	@a = split(/\t/,$line);
	if(!(defined $CHR2SEQ{$a[1]})) {
	    next;
	}
	$seq = $a[4];
	while($seq =~ /^([^+]*)\+/) {  # removing the insertions
	    $pref = $1;
	    $seq =~ s/^$pref\+[^+]+\+/$pref/;
	}
	$strand = $a[3];
	$chr = $a[1];
	@SPANS = split(/, /,$a[2]);
	@SEQ = split(/:/, $seq);
	$snt = $a[0];
	$snt =~ s/seq.//;
#	print STDERR "2:seq.$snt\n";
	for($i=0; $i<@SPANS-1; $i++) {
	    @c1 = split(/-/,$SPANS[$i]);
	    @c2 = split(/-/,$SPANS[$i+1]);
	    $elen1 = $c1[1] - $c1[0] + 1;
	    $elen2 = $c2[1] - $c2[0] + 1;
	    $ilen = $c2[0] - $c1[1] - 1;
	    $istart = $c1[1]+1;
	    $iend = $c2[0]-1;
	    $altintron="";
	    if($ilen >= $minintron) {
		$intron = $chr . ":" . $istart . "-" . $iend;
		$allintrons{$intron} = 1;
		if(!(defined $amb{$intron})) {
		    $SEQ[$i] =~ /(.)$/;
		    $leftexon_lastbase = $1;
		    $SEQ[$i+1] =~ /^(.)/;
		    $rightexon_firstbase = $1;
		    $intron_firstbase = substr($CHR2SEQ{$chr}, $istart-1, 1);
		    $intron_lastbase = substr($CHR2SEQ{$chr}, $iend-1, 1);
		    $splice_signal_upstream = substr($CHR2SEQ{$chr}, $istart-1, 2);
		    $splice_signal_downstream = substr($CHR2SEQ{$chr}, $iend-2, 2);
		    if(($splice_signal_upstream eq "GT" && $splice_signal_downstream eq "AG") || ($splice_signal_upstream eq "CT" && $splice_signal_downstream eq "AC")) {
			$goodsplicesignal{$intron} = 1;
		    } else {
			$goodsplicesignal{$intron} = 0;
		    }
		    if($leftexon_lastbase eq $intron_lastbase) {
			$istart_alt = $istart-1;
			$iend_alt = $iend-1;
			$altintron = $chr . ":" . $istart_alt . "-" . $iend_alt;
			$amb{$intron}=1;  # amb for ambiguous
			$amb{$altintron}=1;
			$allintrons{$intron} = 1;
		    }
		    if($rightexon_firstbase eq $intron_firstbase) {
			$istart_alt = $istart+1;
			$iend_alt = $iend+1;
			$altintron = $chr . ":" . $istart_alt . "-" . $iend_alt;
			$amb{$intron}=1;  # amb for ambiguous
			$amb{$altintron}=1;
			$allintrons{$intron} = 1;
		    }
		}
#		    print "elen1 = $elen1\n";
#		    print "elen2 = $elen2\n";
		if($elen1 <=7 || $elen2 <= 7) {
		    $badoverlapNU{$intron}++;
		    if($altintron =~ /\S/) {
			$badoverlapNU{$altintron}++;			    
		    }
		} else {
		    $goodoverlapNU{$intron}++;
		    if($altintron =~ /\S/) {
			$goodoverlapNU{$altintron}++;			    
		    }
		}
	    }
	}
    }
    close(INFILE);
    print STDERR "finished NU\n";
}
