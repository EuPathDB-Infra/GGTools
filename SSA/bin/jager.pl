use strict;
$|=1;

if(@ARGV<4) {
    die "
Usage: perl transcript_inferrer.pl <junctions file> <coverage file> <rum unique> <annot file> [options]

Where:

  <junctions file>  : The high-quality junctions file output by RUM, sorted
                      by chromosome (it should be sorted already when it comes
                      out of RUM).

  <coverage file>   : The coverage file output by RUM, sorted by chromosome
                      (it should be sorted already when it comes out of RUM).

  <rum unique file> : RUM_Unique file output by RUM, the one that is sorted
                       by chromosome/location.

  <annot file>      : Transcript models file, in the format of the RUM gene
                      info file.

Options:

  -minscore : don't use junctions unless they have at least this score
              (default = 1).  Note: this will only be applied to things
              with coverage at the junction of at least 5 times the minscore.

  -maxexon  : Don't infer exons larger than this (default = 2000 bp)

";
}
my $min_intron = 35;

my @junctions_file;
my $junctionsinfile = $ARGV[0];
my $covinfile = $ARGV[1];
my $rumfile = $ARGV[2];
my $annotfile = $ARGV[3];
my %count_coverage_in_span_cache;
my %ave_coverage_in_span_cache;
my $minscore = 1;
my $maxexon = 2000;
my @inferredTranscript;

for(my $i=4; $i<@ARGV; $i++) {
    my $optionrecognized = 0;
    my $double = "false";
    if($ARGV[$i] eq "-minscore") {
	$minscore = $ARGV[$i+1];
	$i++;
	$optionrecognized = 1;
	$double = "true";
    }
    if($ARGV[$i] eq "-maxexon") {
	$maxexon = $ARGV[$i+1];
	$i++;
	$optionrecognized = 1;
	$double = "true";
    }
    if($optionrecognized == 0) {
	if($double eq "true") {
	    my $temp = $ARGV[$i-1];
	    die "\nERROR: option '$temp $ARGV[$i]' not recognized\n";
	} else {
	    die "\nERROR: option '$ARGV[$i]' not recognized\n";
	}
    }
}

open(INFILE, $junctionsinfile) or die "Error: cannot open '$junctionsinfile' for reading\n\n";
my $junctions_ref = &filter_junctions_file();
my %junctions = %{$junctions_ref};
close(INFILE);

#foreach my $chr (keys %junctions) {
#    my $N = @{$junctions{$chr}};
#    for(my $i=0; $i<$N; $i++) {
#	print "junctions{$chr}[$i] = $junctions{$chr}[$i]\n";
#    }
#}


open(INFILE, $annotfile) or die "Error: cannot open '$annotfile' for reading.\n\n";

# read in the transcript models

my %TRANSCRIPT;
my %EXON_temp;
my %INTRON_temp;
my %cnt;
my @A;
my @B;
my %tcnt;
my %ecnt;
my %icnt;

while(my $line = <INFILE>) {
    chomp($line);
    my @a = split(/\t/,$line);

    $a[5] =~ s/\s*,\s*$//;
    $a[6] =~ s/\s*,\s*$//;
    my $chr = $a[0];
    $tcnt{$chr}=$tcnt{$chr}+0;
    $TRANSCRIPT{$chr}[$tcnt{$chr}]{strand} = $a[1];
    $TRANSCRIPT{$chr}[$tcnt{$chr}]{num} = $a[4];
    $TRANSCRIPT{$chr}[$tcnt{$chr}]{start} = $a[2]+1;  # add one to convert to one-based coords
    $TRANSCRIPT{$chr}[$tcnt{$chr}]{end} = $a[3];
    my @s = split(/,/,$a[5]);
    my @e = split(/,/,$a[6]);
    my @c;
    my $transcript_length=0;
    for(my $i=0; $i<@s; $i++) {
	$TRANSCRIPT{$chr}[$tcnt{$chr}]{coords}[2*$i]=$s[$i]+1;  # add one to convert to one-based coords
	$TRANSCRIPT{$chr}[$tcnt{$chr}]{coords}[2*$i+1]=$e[$i];
	my $S = $s[$i]+1;
	my $E = $chr . ":" . $S . "-" . $e[$i];
	$transcript_length = $transcript_length + $e[$i] - $S + 1;
	$EXON_temp{$chr}{$E}{start} = $S;
	$EXON_temp{$chr}{$E}{end} = $e[$i];
	if($i < @s-1) {
	    my $s2 = $e[$i]+1;
	    my $e2 = $s[$i+1];
	    my $E = $chr . ":" . $s2 . "-" . $e2;
	    $INTRON_temp{$chr}{$E}{start} = $s2;
	    $INTRON_temp{$chr}{$E}{end} = $e2;
	}
    }
    $TRANSCRIPT{$chr}[$tcnt{$chr}]{length} = $transcript_length;
    $TRANSCRIPT{$chr}[$tcnt{$chr}]{id} = $a[7];
    $tcnt{$chr}++;
}
close(INFILE);

my %EXON;
my %EXONS;
foreach my $chr (sort {cmpChrs($a,$b)} keys %EXON_temp) {
    $ecnt{$chr} = 0;
    foreach my $exon (sort {$EXON_temp{$chr}{$a}{start} <=> $EXON_temp{$chr}{$b}{start}} keys %{$EXON_temp{$chr}}) {
	$EXON{$chr}[$ecnt{$chr}]{start} = $EXON_temp{$chr}{$exon}{start};
	$EXON{$chr}[$ecnt{$chr}]{end} = $EXON_temp{$chr}{$exon}{end};
	$EXON{$chr}[$ecnt{$chr}]{exon} = $exon;
	my $s = $EXON_temp{$chr}{$exon}{start};
	my $e = $EXON_temp{$chr}{$exon}{end};
	my $exon = "$chr:$s-$e";
	$EXONS{$exon}=1;
	$ecnt{$chr}++;
    }
}
my %INTRON;
foreach my $chr (sort {cmpChrs($a,$b)} keys %INTRON_temp) {
    $icnt{$chr} = 0;
    foreach my $intron (sort {$INTRON_temp{$chr}{$a}{start} <=> $INTRON_temp{$chr}{$b}{start}} keys %{$INTRON_temp{$chr}}) {
	$INTRON{$chr}[$icnt{$chr}]{start} = $INTRON_temp{$chr}{$intron}{start};
	$INTRON{$chr}[$icnt{$chr}]{end} = $INTRON_temp{$chr}{$intron}{end};
	$INTRON{$chr}[$icnt{$chr}]{intron} = $intron;
	$icnt{$chr}++;
    }
}

open(COVFILE, $covinfile) or die "Error: cannot open '$covinfile' for reading\n\n";
my $line = <INFILE>;
chomp($line);
my @a = split(/\t/,$line);
my $chr = $a[0];
my %prev_exon_start;
my $count;
my %junction_score;
my %working_chr;
my @coverage;

# Going to change everything into one-based (inclusive) coordinates, then will change back
# to print.  Did this to keep from going crazy.

my %working_on;

my $cnt = 0;
my %junctionstarts;
my %junctionends;
my @junction_starts;
my @junction_ends;
my @coverage;
my %exonstart2ends;
my %exonend2starts;
my %putative_exons;
my %JUNCTIONS;
my %startloc2exons;
my %endloc2exons;
my %adjacent_and_connected;
my @putative_exon_array;
my %associated;

open(RUMFILE, $rumfile);

foreach my $chr (sort {cmpChrs($a,$b)} keys %junctions) {

    print STDERR "working on chromosome '$chr'\n";

    undef @coverage;
    undef %exonstart2ends;
    undef %exonend2starts;
    undef %junctionstarts;
    undef %junctionends;
    undef @junction_starts;
    undef @junction_ends;
    undef %putative_exons;
    undef %JUNCTIONS;
    undef %startloc2exons;
    undef %endloc2exons;
    undef %adjacent_and_connected;
    undef @putative_exon_array;
    undef %associated;

    my $flag2 = 0;
    # read in the coverage file for one chromosome
    while($flag2 == 0) {
	my $line = <COVFILE>;
	if($line =~ /track type/) {
	    next;
	}
	chomp($line);
	my @a = split(/\t/,$line);
	if($a[0] eq $chr) {
	    for(my $j=$a[1]+1; $j<=$a[2]; $j++) {
		$coverage[$j] = $a[3];
	    }
	} else {
	    # reset the file handle so the last line read will be read again                 
	    my $len = -1 * (1 + length($line));
	    seek(COVFILE, $len, 1);
	    $flag2 = 1;
	}
    }

    my $N = @{$junctions{$chr}};
    for(my $k=0; $k<$N; $k++) {
	$line = $junctions{$chr}[$k];
	@a = split(/\t/,$line);
	my $score = $a[3];
	my $blocksizes = $a[10];
	$blocksizes =~ s/\s*,\s*$//;
	my @B = split(/,/,$blocksizes);
	my $offsets = $a[11];
	$offsets =~ s/\s*,\s*$//;
	my @O = split(/,/,$offsets);
	$junctionstarts{$a[1] + $B[0]} = 1;
	$junctionends{$a[1] + $O[1] + 1} = 1; # added one to change to one-based
	my $x = $a[1] + $B[0];
	my $y = $a[1] + $O[1] + 1;
	my $J = "$chr:$x-$y\n";
	$JUNCTIONS{$J} = $x;
    }

    my $C=0;
    foreach my $s (sort {$a<=>$b} keys %junctionstarts) {
	$junction_starts[$C] = $s;
	$C++;
    }
    $C=0;
    foreach my $s (sort {$a<=>$b} keys %junctionends) {
	$junction_ends[$C] = $s;
	$C++;
    }
    my $start_index = 0;

    for(my $je=0; $je<@junction_ends; $je++) {
	my $flag = 0;
	my $js = $start_index;
	if($js % 1000 == 0) {
	    print STDERR "finished $js\n";
	}

	while($flag == 0) {
	    if($start_index >= @junction_starts) {
		$flag = 1;
		next;
	    }

#	    print "--------\njunction_starts[$js] = $junction_starts[$js]\n";
#	    print "junction_ends[$je] = $junction_ends[$je]\n";

	    if($junction_starts[$js] <= $junction_ends[$je]) {
		$start_index++;
	    }
	    if($junction_starts[$js] > $junction_ends[$je] + $maxexon) {
		$flag = 1;
		next;
	    }
	    if($junction_ends[$je] < $junction_starts[$js] && $junction_starts[$js] <= $junction_ends[$je] + $maxexon) {
		
		my $cov_yes = 0;
		my $annot_yes = 0;
		
		# 1) see if there is coverage across span $junction_ends[$je] to $junction_starts[$js]
		my $N = &count_coverage_in_span($junction_ends[$je], $junction_starts[$js], 1);
		if($N == 0) {
		    $cov_yes = 1;
		}
		
		# 2) see if the span $junction_ends[$je] to $junction_starts[$js] is annotated as an exon
		
		my $exon = "$chr:$junction_ends[$je]-$junction_starts[$js]";
		if($EXONS{$exon} + 0 == 1) {
		    $annot_yes = 1;
		}
		
		# 3) if 2 *or* 3 are 'yes' then call this an exon
		if($cov_yes == 1 || $annot_yes == 1) {
		    my $exon = "$chr:$junction_ends[$je]-$junction_starts[$js]";
		    $putative_exons{$exon}=$junction_ends[$je];
		    if(defined $exonstart2ends{$junction_ends[$je]}) {
			$N = @{$exonstart2ends{$junction_ends[$je]}};
		    } else {
			$N = 0;
		    }
		    $exonstart2ends{$junction_ends[$je]}[$N] = $junction_starts[$js];
		    if(defined $exonend2starts{$junction_starts[$js]}) {
			$N = @{$exonend2starts{$junction_starts[$js]}}+0;
		    } else {
			$N = 0;
		    }
		    $exonend2starts{$junction_starts[$js]}[$N] = $junction_ends[$je];
 		}
	    }		
	    $js++;
	}
    }

    # use %exonstart2ends and %exonend2starts to remove exons that span two exons:
    # in other words remove 2 in cases like the below, where 'x' represents exon:
    # 1) ---xxxxxx--------xxxxxx----
    # 2) ---xxxxxxxxxxxxxxxxxxxx----

    foreach my $exon (keys %putative_exons) {
	$exon =~ /^.*:(\d+)-(\d+)$/;
	my $S = $1;
	my $E = $2;
	my $N1 = @{$exonstart2ends{$S}};
	my $N2 = @{$exonend2starts{$E}};
	for(my $i=0; $i<$N1; $i++) {
	    for(my $j=0; $j<$N2; $j++) {
		if($exonstart2ends{$S}[$i] + $min_intron <= $exonend2starts{$E}[$j]) {
		    delete $putative_exons{$exon};
		}
	    }	    
	}
    }

    $cnt = 0;
    foreach my $exon (sort {$putative_exons{$a}<=>$putative_exons{$b}} keys %putative_exons) {
	$putative_exon_array[$cnt] = $exon;
	$cnt++;
	$exon =~ /^(.*):(\d+)-(\d+)$/;
	my $S = $2;
	my $E = $3;
	my $N = 0;
	if(defined $startloc2exons{$S}) {
	    $N = @{$startloc2exons{$S}};
	} else {
	    $N = 0;
	}
	$startloc2exons{$S}[$N] = $exon;

	if(defined $endloc2exons{$E}) {
	    $N = @{$endloc2exons{$E}};
	} else {
	    $N = 0;
	}
	$endloc2exons{$E}[$N] = $exon;
    }

    foreach my $J (keys %JUNCTIONS) {
	$J =~ /^(.*):(\d+)-(\d+)$/;
	my $S = $2;
	my $E = $3;
	if(defined $endloc2exons{$S} && defined $startloc2exons{$E}) {
	    my $N1 = @{$endloc2exons{$S}};
	    my $N2 = @{$startloc2exons{$E}};
	    for(my $i=0; $i<$N1; $i++) {
		for(my $j=0; $j<$N2; $j++) {
		    # The start of the junction is the end of an exon and the end of the junction
		    # is the start of an exon.
		    $adjacent_and_connected{$endloc2exons{$S}[$i]}[$j] = $startloc2exons{$E}[$j];
		}
	    }
	}
    }

# DEBUG STUFF
#    foreach my $exon (keys %adjacent_and_connected) {
#	my $N = @{$adjacent_and_connected{$exon}};
#	print "------\n ***  $exon\n";
#	for(my $i=0; $i<$N; $i++) {
#	    my $x = $adjacent_and_connected{$exon}[$i];
#	    print "$x\n";
#	}
#    }
# DEBUG STUFF

    # for reach read pair:
    # 1) make list of upstream exons from upstream read
    #     - if read has one span and exon completely contains it, or
    #       if exon contains end span and shares (the right) break point with that span, or
    #       if exon contains internal span and shares both break points with that span.
    # 2) make list of downstream exons the same way, from downstream read
    # 3) associate:
    #    i) exons upstream with exons downstream
    #   ii) exons upstream with exons upstream if separated by at least one exon
    #  iii) exons downstream with exons downstream if separated by at least one exon

    my $exon_index_start = 0;
    my %exons_hitting_read;
    while(1 == 1) {  # one pass through this loop is one read (pair) in RUM_Unique
	$line = <RUMFILE>;
	chomp($line);
	if($line eq '') {
	    last;
	}

	if(!($line =~ /\t$chr\t/)) {
	    # reset the file handle so the last line read will be read again
	    my $len = -1 * (1 + length($line));
	    seek(RUMFILE, $len, 1);
	    last;
	}

	my @a = split(/\t/,$line);
	my $STRAND = $a[3];
	$a[0] =~ /(\d+)/;
	my $seqnum1 = $1;
	$a[2] =~ /^(\d+)-/;
	my $start = $1;
	my $end;
	my $line2;
	my $upstream_spans="";
	my $downstream_spans="";
	my $spans;
	my $merged = "false";
	my $aonly = "false";
	my $bonly = "false";
	my $unmerged_pair = "false";
	if($a[0] =~ /a/) {
	    $line2 = <RUMFILE>;
	    chomp($line2);
	    my @b = split(/\t/,$line2);
	    $b[0] =~ /(\d+)/;
	    my $seqnum2 = $1;
	    if($seqnum1 == $seqnum2 && $b[0] =~ /b/ && $a[0] =~ /a/) {  # it's a proper pair
		if($a[3] eq "+") {
		    $b[2] =~ /-(\d+)$/;
		    $end = $1;
		    $upstream_spans = $a[2];
		    $downstream_spans = $b[2];
		} else {
		    $b[2] =~ /^(\d+)-/;
		    $start = $1;
		    $a[2] =~ /-(\d+)$/;
		    $end = $1;
		    $upstream_spans = $b[2];
		    $downstream_spans = $a[2];
		}
		$unmerged_pair = "true";
	    } else {  # it's an 'a' without a 'b'
		$a[2] =~ /-(\d+)$/;
		$end = $1;
		# reset the file handle so the last line read will be read again
		my $len = -1 * (1 + length($line2));
		seek(RUMFILE, $len, 1);
		$upstream_spans = $a[2];
		$aonly = "true";
	    }
	}
	if($a[0] =~ /b/) { # it's a 'b' without an 'a'
	    $a[2] =~ /-(\d+)$/;
	    $end = $1;	    
	    $upstream_spans = $a[2];
	    $bonly = "true";
	}
	if(!($a[0] =~ /a/) && !($a[0] =~ /b/)) { # it's a merged pair
	    $a[2] =~ /-(\d+)$/;
	    $end = $1;	    
	    $merged = "true";
	    $upstream_spans = $a[2];
	}

	# first remove indels from spans, because we just want introns.
	my @temp = split(/, /, $upstream_spans);
	for(my $i=0; $i<@temp-1; $i++) {
	    $temp[$i] =~ /-(\d+)/;
	    my $x = $1;
	    $temp[$i+1] =~ /(\d+)-/;
	    my $y = $1;
	    if($y - $x < $min_intron) {
		$upstream_spans =~ s/-$x, $y-/-/;
	    }
	}
	@temp = split(/, /, $downstream_spans);
	for(my $i=0; $i<@temp-1; $i++) {
	    $temp[$i] =~ /-(\d+)/;
	    my $x = $1;
	    $temp[$i+1] =~ /(\d+)-/;
	    my $y = $1;
	    if($y - $x < $min_intron) {
		$downstream_spans =~ s/-$x, $y-/-/;
	    }
	}

	# in the following three cases there is no useful paired end info, so just skip
	if($aonly eq "true" && !($upstream_spans =~ /,/)) {
	    next;
	}
	if($bonly eq "true" && !($downstream_spans =~ /,/)) {
	    next;
	}
	if($merged eq "true" && !($spans =~ /,/)) {
	    next;
	}
#	print "-----------\n";
#	print "line = $line\n";
#	print "line2 = $line2\n";
#	print "unmerged_pair = $unmerged_pair\n";
#	print "aonly = $aonly\n";
#	print "bonly = $bonly\n";
#	print "merged = $merged\n";
#	print "upstream_spans = $upstream_spans\n";
#	print "downstream_spans = $downstream_spans\n";


	# make list of upstream exons

	my $exon_index = $exon_index_start - 1;
	my $flag = 0;
	undef %exons_hitting_read;
#	print "upstream search...\n";
	while(2 == 2) { # this loop scrolls through the exons that overlap the upstream spans
	    $exon_index++;
#	    print "putative_exon_array[$exon_index] = $putative_exon_array[$exon_index]\n";
	    $putative_exon_array[$exon_index] =~ /^(.*):(\d+)-(\d+)/;
	    my $chr = $1;
	    my $exonstart = $2;
	    my $exonend = $3;
	    if($flag == 0 && $exonend < $start) {
		$exon_index_start++;
	    }
	    $flag = 1;
	    if($exon_index_start >= @putative_exon_array) {
		last;
	    }
	    if($exonstart > $end) {
		last;
	    }
	    if($exonend < $start) {
		next;
	    }	    
	    # got here means $putative_exon_array[$exon_index] is an exon that overlaps
	    # with ($start, $end) of read(s) span

	    my $flag_s = 0;
	    if(!($upstream_spans =~ /,/)) {
		my @s = split(/-/, $upstream_spans);
		if($s[0] >= $exonstart && $s[1] <= $exonend) {
		    $flag_s = 1;
		}
	    } else {
		my @SPANS = split(/, /, $upstream_spans);
		for(my $span=0; $span<@SPANS; $span++) {
		    my @s = split(/-/, $SPANS[$span]);
		    if($span == 0 && $s[0] >= $exonstart && $s[1] == $exonend) {
			$flag_s = 1;
		    } elsif($span == @SPANS-1 && $s[0] == $exonstart && $s[1] <= $exonend) {
			$flag_s = 1;
		    } elsif($s[0] == $exonstart && $s[1] == $exonend) {
			$flag_s = 1;
		    }
		}
	    }
	    if($flag_s == 1) {
		$exons_hitting_read{$putative_exon_array[$exon_index]}[0]=$exonstart;
		$exons_hitting_read{$putative_exon_array[$exon_index]}[1]=$exonend;
	    }
	}
#	foreach my $exon (keys %exons_hitting_read) {
#	    print "upstream exon: $exon\n";
#	}

	# make list of downstream exons

	$exon_index = $exon_index_start - 1;
	$flag = 0;
#	print "downstream search...\n";
	while(2 == 2) { # this loop scrolls through the exons that overlap the upstream spans
	    $exon_index++;
#	    print "putative_exon_array[$exon_index] = $putative_exon_array[$exon_index]\n";
	    $putative_exon_array[$exon_index] =~ /^(.*):(\d+)-(\d+)/;
	    my $chr = $1;
	    my $exonstart = $2;
	    my $exonend = $3;
	    $flag = 1;
	    if($exon_index_start >= @putative_exon_array) {
		last;
	    }
	    if($exonstart > $end) {
		last;
	    }
	    if($exonend < $start) {
		next;
	    }	    
	    # got here means $putative_exon_array[$exon_index] is an exon that overlaps
	    # with ($start, $end) of read(s) span

	    my $flag_s = 0;
	    if(!($downstream_spans =~ /,/)) {
		my @s = split(/-/, $downstream_spans);
		if($s[0] >= $exonstart && $s[1] <= $exonend) {
		    $flag_s = 1;
		}
	    } else {
		my @SPANS = split(/, /, $downstream_spans);
		for(my $span=0; $span<@SPANS; $span++) {
		    my @s = split(/-/, $SPANS[$span]);
		    if($span == 0 && $s[0] >= $exonstart && $s[1] == $exonend) {
			$flag_s = 1;
		    } elsif($span == @SPANS-1 && $s[0] == $exonstart && $s[1] <= $exonend) {
			$flag_s = 1;
		    } elsif($s[0] == $exonstart && $s[1] == $exonend) {
			$flag_s = 1;
		    }
		}
	    }
	    if($flag_s == 1) {
		$exons_hitting_read{$putative_exon_array[$exon_index]}[0]=$exonstart;
		$exons_hitting_read{$putative_exon_array[$exon_index]}[1]=$exonend;
	    }
	}
#	foreach my $exon (keys %exons_hitting_read) {
#	    print "downstream exon: $exon\n";
#	}

	# Now have list of upstream and downstream exons, these are all things
	# that could be in the same transcript, so now make hash of the associations
	# between all pairs of them.  Note even two upstream or two downstream carry the
	# same info as one up and one down, so record all.

	foreach my $exon1 (keys %exons_hitting_read) {
	    foreach my $exon2 (keys %exons_hitting_read) {	    
		if($exons_hitting_read{$exon1}[1] + $min_intron <= $exons_hitting_read{$exon2}[0]) {
		    $associated{$exon1} = $exon2;
		    print "$exon1 --> $exon2\n";
		}
	    }
	}
    }
}



sub count_coverage_in_span () {
    # This will return the number of bases in the span that
    # have coverage no more than $coverage_cutoff
    my ($start, $end, $coverage_cutoff) = @_;
    my $tmp = $start . ":" . $end . ":" . $coverage_cutoff;
    if(defined $count_coverage_in_span_cache{$tmp}) {
	return $count_coverage_in_span_cache{$tmp};
    }
    my $num_below=0;
    for(my $i=$start; $i<=$end; $i++) {
	if($coverage[$i] < $coverage_cutoff) {
	    $num_below++;
	}
    }
    $count_coverage_in_span_cache{$tmp}=$num_below;
    return $num_below;
}

sub ave_coverage_in_span () {
    # This will return the average depth over bases in the span
    my ($start, $end, $coverage_cutoff) = @_;
    my $tmp = $start . ":" . $end . ":" . $coverage_cutoff;
    my $sum = 0;
    if(defined $ave_coverage_in_span_cache{$tmp}) {
	return $ave_coverage_in_span_cache{$tmp};
    }
    for(my $i=$start; $i<=$end; $i++) {
	$sum = $sum + $coverage[$i];
    }
    my $ave = $sum / ($end - $start + 1);
    $ave_coverage_in_span_cache{$tmp}=$ave;
    return $ave;
}

# Initial filtering to remove low scoring stuff that does not seem
# like it should be part of a transcript:
#
# Score <= 2, unknown, overlaps with something with score >= 20,
# unless it is a perfect alternate splice form without being
# ridiculously long.
# (might want to play with 20 to get whatever works best)

sub filter_junctions_file () {
    my %junction_num;
    my %kept_junctions;
    my $scorefilter = 3;
    my $scorefilter_max = 20;
    my %junctions_file;
    my %starts;
    my %ends;
    while(my $line = <INFILE>) {
	chomp($line);
	my @a = split(/\t/,$line);
	$junction_num{$a[0]}=$junction_num{$a[0]}+0;
	$junctions_file{$a[0]}[$junction_num{$a[0]}][0] = $a[1];
	$junctions_file{$a[0]}[$junction_num{$a[0]}][1] = $a[2];
	$junctions_file{$a[0]}[$junction_num{$a[0]}][2] = $a[4];
	$junctions_file{$a[0]}[$junction_num{$a[0]}][3] = $line;
	if($a[4] >= $scorefilter) {
	    $starts{$a[0]}{$a[1]}=1;
	    $ends{$a[0]}{$a[2]}=1;
	}
	$junction_num{$a[0]}++;
    }
#    foreach my $chr (keys %junction_num) {
#	foreach my $start (sort {$a<=>$b} keys %{$starts{$chr}}) {
#	    print "starts{$chr}{$start}=$starts{$chr}{$start}\n";
#	}
#    }
    foreach my $chr (keys %junction_num) {
	my $start = 0;
	my $kept_counter=0;
	for(my $i=0; $i<$junction_num{$chr}; $i++) {
#	    print "junctions_file{$chr}[$i][2]=$junctions_file{$chr}[$i][2]\n";
#	    print "$junctions_file{$chr}[$i][3]\n";
#	    print "junctions_file{$chr}[$i][0]=$junctions_file{$chr}[$i][0]\n";
#	    print "junctions_file{$chr}[$i][0]=$junctions_file{$chr}[$i][0]\n";
	    if($junctions_file{$chr}[$i][2] > $scorefilter) {
		$kept_junctions{$chr}[$kept_counter] = $junctions_file{$chr}[$i][3];
		$kept_counter++;
	    } elsif($starts{$chr}{$junctions_file{$chr}[$i][0]}+0==1 &&
		    $ends{$chr}{$junctions_file{$chr}[$i][1]}+0==1) {
		$kept_junctions{$chr}[$kept_counter] = $junctions_file{$chr}[$i][3];
		$kept_counter++;
	    } elsif($junctions_file{$chr}[$i][3] =~ /24,116,205/ || $junctions_file{$chr}[$i][3] =~ /16,78,139/) {
		$kept_junctions{$chr}[$kept_counter] = $junctions_file{$chr}[$i][3];
		$kept_counter++;		
	    } else {
		my $flag = 0;
		for(my $j=$start; $j<$junction_num{$chr}; $j++) {
		    if($junctions_file{$chr}[$j][2] >= $scorefilter_max &&
		       $junctions_file{$chr}[$j][0] <= $junctions_file{$chr}[$i][1] &&
		       $junctions_file{$chr}[$j][1] >= $junctions_file{$chr}[$i][0]) {
			   $flag = 1;
		    }
		    if($junctions_file{$chr}[$i][1] < $junctions_file{$chr}[$j][0]) {
			# we've checked everything that can overlap, set $j to jump out of loop
			$j = $junction_num{$chr};
		    }
		}
		if($flag == 0) {
		    $kept_junctions{$chr}[$kept_counter] = $junctions_file{$chr}[$i][3];
		    $kept_counter++;
		}
		while($junctions_file{$chr}[$start][1] < $junctions_file{$chr}[$i][0] && $start < $junction_num{$chr}) {
		    $start++;
		}		    
	    }
	}
    }
    return \%kept_junctions;
}

sub cmpChrs () {
    my $a2_c = lc($b);
    my $b2_c = lc($a);
    if($a2_c =~ /^\d+$/ && !($b2_c =~ /^\d+$/)) {
        return 1;
    }
    if($b2_c =~ /^\d+$/ && !($a2_c =~ /^\d+$/)) {
        return -1;
    }
    if($a2_c =~ /^[ivxym]+$/ && !($b2_c =~ /^[ivxym]+$/)) {
        return 1;
    }
    if($b2_c =~ /^[ivxym]+$/ && !($a2_c =~ /^[ivxym]+$/)) {
        return -1;
    }
    if($a2_c eq 'm' && ($b2_c eq 'y' || $b2_c eq 'x')) {
        return -1;
    }
    if($b2_c eq 'm' && ($a2_c eq 'y' || $a2_c eq 'x')) {
        return 1;
    }
    if($a2_c =~ /^[ivx]+$/ && $b2_c =~ /^[ivx]+$/) {
        $a2_c = "chr" . $a2_c;
        $b2_c = "chr" . $b2_c;
    }
   if($a2_c =~ /$b2_c/) {
	return -1;
    }
    if($b2_c =~ /$a2_c/) {
	return 1;
    }
    # dealing with roman numerals starts here

    if($a2_c =~ /chr([ivx]+)/ && $b2_c =~ /chr([ivx]+)/) {
	$a2_c =~ /chr([ivx]+)/;
	my $a2_roman = $1;
	$b2_c =~ /chr([ivx]+)/;
	my $b2_roman = $1;
	my $a2_arabic = arabic($a2_roman);
    	my $b2_arabic = arabic($b2_roman);
	if($a2_arabic > $b2_arabic) {
	    return -1;
	} 
	if($a2_arabic < $b2_arabic) {
	    return 1;
	}
	if($a2_arabic == $b2_arabic) {
	    my $tempa = $a2_c;
	    my $tempb = $b2_c;
	    $tempa =~ s/chr([ivx]+)//;
	    $tempb =~ s/chr([ivx]+)//;
	    my %temphash;
	    $temphash{$tempa}=1;
	    $temphash{$tempb}=1;
	    foreach my $tempkey (sort {cmpChrs($a,$b)} keys %temphash) {
		if($tempkey eq $tempa) {
		    return 1;
		} else {
		    return -1;
		}
	    }
	}
    }

    if($b2_c =~ /chr([ivx]+)/ && !($a2_c =~ /chr([a-z]+)/) && !($a2_c =~ /chr(\d+)/)) {
	return -1;
    }
    if($a2_c =~ /chr([ivx]+)/ && !($b2_c =~ /chr([a-z]+)/) && !($b2_c =~ /chr(\d+)/)) {
	return 1;
    }
    # roman numerals ends here
    if($a2_c =~ /chr(\d+)$/ && $b2_c =~ /chr.*_/) {
        return 1;
    }
    if($b2_c =~ /chr(\d+)$/ && $a2_c =~ /chr.*_/) {
        return -1;
    }
    if($a2_c =~ /chr([a-z])$/ && $b2_c =~ /chr.*_/) {
        return 1;
    }
    if($b2_c =~ /chr([a-z])$/ && $a2_c =~ /chr.*_/) {
        return -1;
    }
    if($a2_c =~ /chr(\d+)/) {
        my $numa = $1;
        if($b2_c =~ /chr(\d+)/) {
            my $numb = $1;
            if($numa < $numb) {return 1;}
	    if($numa > $numb) {return -1;}
	    if($numa == $numb) {
		my $tempa = $a2_c;
		my $tempb = $b2_c;
		$tempa =~ s/chr\d+//;
		$tempb =~ s/chr\d+//;
		my %temphash;
		$temphash{$tempa}=1;
		$temphash{$tempb}=1;
		foreach my $tempkey (sort {cmpChrs($a,$b)} keys %temphash) {
		    if($tempkey eq $tempa) {
			return 1;
		    } else {
			return -1;
		    }
		}
	    }
        } else {
            return 1;
        }
    }
    if($a2_c =~ /chrx(.*)/ && ($b2_c =~ /chr(y|m)$1/)) {
	return 1;
    }
    if($b2_c =~ /chrx(.*)/ && ($a2_c =~ /chr(y|m)$1/)) {
	return -1;
    }
    if($a2_c =~ /chry(.*)/ && ($b2_c =~ /chrm$1/)) {
	return 1;
    }
    if($b2_c =~ /chry(.*)/ && ($a2_c =~ /chrm$1/)) {
	return -1;
    }
    if($a2_c =~ /chr\d/ && !($b2_c =~ /chr[^\d]/)) {
	return 1;
    }
    if($b2_c =~ /chr\d/ && !($a2_c =~ /chr[^\d]/)) {
	return -1;
    }
    if($a2_c =~ /chr[^xy\d]/ && (($b2_c =~ /chrx/) || ($b2_c =~ /chry/))) {
        return -1;
    }
    if($b2_c =~ /chr[^xy\d]/ && (($a2_c =~ /chrx/) || ($a2_c =~ /chry/))) {
        return 1;
    }
    if($a2_c =~ /chr(\d+)/ && !($b2_c =~ /chr(\d+)/)) {
        return 1;
    }
    if($b2_c =~ /chr(\d+)/ && !($a2_c =~ /chr(\d+)/)) {
        return -1;
    }
    if($a2_c =~ /chr([a-z])/ && !($b2_c =~ /chr(\d+)/) && !($b2_c =~ /chr[a-z]+/)) {
        return 1;
    }
    if($b2_c =~ /chr([a-z])/ && !($a2_c =~ /chr(\d+)/) && !($a2_c =~ /chr[a-z]+/)) {
        return -1;
    }
    if($a2_c =~ /chr([a-z]+)/) {
        my $letter_a = $1;
        if($b2_c =~ /chr([a-z]+)/) {
            my $letter_b = $1;
            if($letter_a lt $letter_b) {return 1;}
	    if($letter_a gt $letter_b) {return -1;}
        } else {
            return -1;
        }
    }
    my $flag_c = 0;
    while($flag_c == 0) {
        $flag_c = 1;
        if($a2_c =~ /^([^\d]*)(\d+)/) {
            my $stem1_c = $1;
            my $num1_c = $2;
            if($b2_c =~ /^([^\d]*)(\d+)/) {
                my $stem2_c = $1;
                my $num2_c = $2;
                if($stem1_c eq $stem2_c && $num1_c < $num2_c) {
                    return 1;
                }
                if($stem1_c eq $stem2_c && $num1_c > $num2_c) {
                    return -1;
                }
                if($stem1_c eq $stem2_c && $num1_c == $num2_c) {
                    $a2_c =~ s/^$stem1_c$num1_c//;
                    $b2_c =~ s/^$stem2_c$num2_c//;
                    $flag_c = 0;
                }
            }
        }
    }
    if($a2_c le $b2_c) {
	return 1;
    }
    if($b2_c le $a2_c) {
	return -1;
    }

    return 1;
}

sub isroman($) {
    my $arg = shift;
    $arg ne '' and
      $arg =~ /^(?: M{0,3})
                (?: D?C{0,3} | C[DM])
                (?: L?X{0,3} | X[LC])
                (?: V?I{0,3} | I[VX])$/ix;
}

sub arabic($) {
    my $arg = shift;
    my %roman2arabic = qw(I 1 V 5 X 10 L 50 C 100 D 500 M 1000);
    my %roman_digit = qw(1 IV 10 XL 100 CD 1000 MMMMMM);
    my @figure = reverse sort keys %roman_digit;
    $roman_digit{$_} = [split(//, $roman_digit{$_}, 2)] foreach @figure;
    isroman $arg or return undef;
    my $last_digit = 1000;
    my $arabic=0;
    foreach (split(//, uc $arg)) {
        my ($digit) = $roman2arabic{$_};
        $arabic -= 2 * $last_digit if $last_digit < $digit;
        $arabic += ($last_digit = $digit);
    }
    $arabic;
}

sub Roman($) {
    my $arg = shift;
    my %roman2arabic = qw(I 1 V 5 X 10 L 50 C 100 D 500 M 1000);
    my %roman_digit = qw(1 IV 10 XL 100 CD 1000 MMMMMM);
    my @figure = reverse sort keys %roman_digit;
    $roman_digit{$_} = [split(//, $roman_digit{$_}, 2)] foreach @figure;
    0 < $arg and $arg < 4000 or return undef;
    my $roman="";
    my $x;
    foreach (@figure) {
        my ($digit, $i, $v) = (int($arg / $_), @{$roman_digit{$_}});
        if (1 <= $digit and $digit <= 3) {
            $roman .= $i x $digit;
        } elsif ($digit == 4) {
            $roman .= "$i$v";
        } elsif ($digit == 5) {
            $roman .= $v;
        } elsif (6 <= $digit and $digit <= 8) {
            $roman .= $v . $i x ($digit - 5);
        } elsif ($digit == 9) {
            $roman .= "$i$x";
        }
        $arg -= $digit * $_;
        $x = $i;
    }
    $roman;
}

sub roman($) {
    lc Roman shift;
}
