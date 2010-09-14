
if(@ARGV < 2) {
    die "
Usage: get_reads.pl <chr> <start> <end> <chrdir> | less -S

Where <chrdir> is a directory with a RUM file that has
been broken into individual chromosomes (this is done for
efficiency).  To do the breaking up for 'RUM_Unique', make
a directory, say 'chrdir' and run the following:
> perl get_reads.pl RUM_Unique -prepare chrdir
After that finishes you can run the program using the
normal four arguments.

Output is piped to 'less -S' so it's easier to read, but
you can also output to a file.  If you name the file with
suffix '.txt' and open it in a web browser it will display
pretty nicely.  You might want to zoom way out to see it
better, - on Firefox control-<minus-sign> shrinks the font
and control-<plus sign> increases it

";
}

$infile = $ARGV[0];
if($ARGV[1] eq "-prepare") {
    $individual_file_dir = $ARGV[2];
    open(INFILE, $ARGV[0]);
    while($line = <INFILE>) {
	chomp($line);
	@a = split(/\t/,$line);
	$chr = $a[1];
	if(!(defined $chromosomes{$chr})) {
	    $chromosomes{$chr}++;
	    open $F1{$chr}, ">" . "$individual_file_dir/" . $chr;
	}
	$FF = $F1{$chr};
	print $FF "$line\n";
    }
    foreach $chr (keys %F1) {
	close($F1{$chr});
    }
    exit();
}

$CHR = $ARGV[0];
$START = $ARGV[1];
$END = $ARGV[2];
$individual_file_dir = $ARGV[3];

#print "$individual_file_dir/$CHR\n";

open(INFILE, "$individual_file_dir/$CHR");
$hit=0;
$upstream_limit=1000000000000;
$downstream_limit=0;
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    $coords = $a[2];
    $strand = $a[3];
    $seq = $a[4];
    $coords =~ /^(\d+)/;
    $start = $1;
    $coords =~ /(\d+)$/;
    $end = $1;
    if($end >= $START && $start <= $END) {
	$hits[$hit][0] = $coords;
	$hits[$hit][1] = $strand;
	$hits[$hit][2] = $seq;
	$hit++;
	while($seq =~ /\+/) {
	    $seq =~ s/^([^+]+)\+[^+]+\+/$1/;
	}
	if($start < $upstream_limit) {
	    $upstream_limit = $start;
	}
	if($end > $downstream_limit) {
	    $downstream_limit = $end;
	}
#	print "$line\n";
    }
}
close(INFILE);

$numhits = $hit;
#print "numhits = $numhits\n";
#print "upstream_limit = $upstream_limit\n";
#print "downstream_limit = $downstream_limit\n";
for($i=0; $i<$numhits; $i++) {
#    print "parsing $i\n";
    @C = split(/, /, $hits[$i][0]);
    $strand = $hits[$i][1];
    @S = split(/:/, $hits[$i][2]);
    for($j=0; $j<@C; $j++) {
	@b = split(/-/, $C[$j]);
	$COORDS[$j][0] = $b[0];
	$COORDS[$j][1] = $b[1];
	@{$SEQ[$j]} = split(//, $S[$j]);
    }
    $coord_span_cnt = 0;
    for($j=$upstream_limit; $j<=$downstream_limit; $j++) {
	if($j >= $COORDS[$coord_span_cnt][0] && $j <= $COORDS[$coord_span_cnt][1]) {
	    $data_exists{$j}++;
	}
	if($j == $COORDS[$coord_span_cnt][1]) {
	    $coord_span_cnt++;
	}
    }
}
$gap = 0;
open(LOG, ">log.txt");
for($j=$upstream_limit; $j<=$downstream_limit; $j++) {
    print LOG "$j\t$data_exists{$j}\n";
    if(!(exists $data_exists{$j}) && $gap == 0) {
	print LOG "here 1\n";
	$gap_start{$j}++;
	$gap = 1;
	$j2 = $j;
    } elsif((exists $data_exists{$j}) && $gap == 1) {
	    print LOG "here 2\n";
	    $gap_end{$j2}=$j;
	    $gap = 0;
    }
}
close(LOG);
foreach $key (sort {$a<=>$b} keys %gap_start) {
    print LOG "$key: $gap_start{$key} : $gap_end{$key}\n";
}

for($i=0; $i<$numhits; $i++) {
    @C = split(/, /, $hits[$i][0]);
    $strand = $hits[$i][1];
    @S = split(/:/, $hits[$i][2]);
    undef @SEQ;
    undef @COORDS;
    for($j=0; $j<@C; $j++) {
	@b = split(/-/, $C[$j]);
	$COORDS[$j][0] = $b[0];
	$COORDS[$j][1] = $b[1];
	@{$SEQ[$j]} = split(//, $S[$j]);
    }
    $coord_span_cnt = 0;
    for($j=$upstream_limit; $j<=$COORDS[@COORDS-1][1]; $j++) {
	if(exists $gap_start{$j}){
#	    print "gap_start{$j} = $gap_start{$j}\n";
#	    print "gap_end{$j} = $gap_end{$j}\n";
	    $gap_length = $gap_end{$j} - $j;
	    $GL = format_large_int($gap_length);
	    print "  <-  $GL bp gap ->  ";
	    $j = $gap_end{$j};
	} else {
	    if($j < $COORDS[$coord_span_cnt][0]) {
		print " ";
	    }
	    if($j >= $COORDS[$coord_span_cnt][0] && $j <= $COORDS[$coord_span_cnt][1]) {
		$k = $j - $COORDS[$coord_span_cnt][0];
		$BASE{$j}{$SEQ[$coord_span_cnt][$k]}++;
		print $SEQ[$coord_span_cnt][$k];
	    }
	    if($j == $COORDS[$coord_span_cnt][1]) {
		$coord_span_cnt++;
	    }
	}
    }
    print "\n";
}

for($j=$upstream_limit; $j<=$downstream_limit; $j++) {
    if(exists $gap_start{$j}){
	$gap_length = $gap_end{$j} - $j;
	$GL = format_large_int($gap_length);
	print "---------------------";
	$j = $gap_end{$j};
    } else {
	print "-";
    }
}
print "--------\n";

for($j=$upstream_limit; $j<=$downstream_limit; $j++) {
    if(exists $gap_start{$j}){
	$gap_length = $gap_end{$j} - $j;
	$GL = format_large_int($gap_length);
	print "  <-  $GL bp gap ->  ";
	$j = $gap_end{$j};
    } else {
	$consensus = "A";
	$A = $BASE{$j}{"A"};
	$max = $A;
	$C = $BASE{$j}{"C"};
	if($C > $max) {
	    $consensus = "C";
	    $max = $C;
	}
	$G = $BASE{$j}{"G"};
	if($G > $max) {
	    $consensus = "G";
	    $max = $G;
	}
	$T = $BASE{$j}{"T"};
	if($T > $max) {
	    $consensus = "T";
	    $max = $T;
	}
	$N = $BASE{$j}{"N"};
	if($N > $max) {
	    $consensus = "N";
	    $max = $N;
	}
	print $consensus;
    }
}
print "\n";
for($j=$upstream_limit; $j<=$downstream_limit; $j++) {
    if(exists $gap_start{$j}){
	$gap_length = $gap_end{$j} - $j;
	$GL = format_large_int($gap_length);
	print "---------------------";
	$j = $gap_end{$j};
    } else {
	print "-";
    }
}
print "--------\n";


sub format_large_int () {
    ($int_f) = @_;
    @a_f = split(//,"$int_f");
    $j_f=0;
    $newint_f = "";
    $n_f = @a_f;
    for($i_f=$n_f-1;$i_f>=0;$i_f--) {
	$j_f++;
	$newint_f = $a_f[$i_f] . $newint_f;
	if($j_f % 3 == 0) {
	    $newint_f = "," . $newint_f;
	}
    }
    $newint_f =~ s/^,//;
    return $newint_f;
}
