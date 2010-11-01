
if(@ARGV < 1) {
    die "
Usage: sam2bed.pl <samfile>

Options: 
    -uniquerecords  : If there are no non-unique mappers.
                      If there are non-unique mappers, the
                      HI and IH tags must be used.

";
}


open(INFILE, $ARGV[0]);

$uniquerecords = "false";
for($i=1; $i<@ARGV; $i++) {
    if($ARGV[$i] eq '-uniquerecords') {
	$uniquerecords = "true";
    }
}

while($line1 = <INFILE>) {
    if($line1 =~ /^@/) {
	next;
    }
    $line2 = <INFILE>;

    chomp($line1);
    chomp($line2);

    if($uniquerecords eq "false") {
	$line1 =~ /IH:i:(\d+)/;
	$num_alignments_this_read = $1;
	if($num_alignments_this_read > 1) {
	    next;
	}
	if(!($line1 =~ /IH:i/ || $line2 =~ /IH:i/)) {
	    next;
	}
    } else {
	$line1 =~ /X0:i:(\d+)/;
	$num_alignments_this_read = $1;
	if($num_alignments_this_read > 1) {
	    next;
	}
    }
    @a = split(/\t/,$line1);
    $a[2] =~ s/:.*//;
    $chr1 = $a[2];
    $cig1 = $a[5];
    $start1 = $a[3];

    @a = split(/\t/,$line2);
    $a[2] =~ s/:.*//;
    $chr2 = $a[2];
    if($chr1 eq "*" && $chr2 eq "*") {
	next;
    }
    $cig2 = $a[5];
    $start2 = $a[3];

    $id = $a[0];

    $spans1 = &cig2spans($cig1, $start1);
    $spans2 = &cig2spans($cig2, $start2);

    if($chr1 eq $chr2) {
	$spans = union($spans1, $spans2);
	@b = split(/, /,$spans);
	for($i=0; $i<@b; $i++) {
	    @c = split(/-/,$b[$i]);
	    $str = "$chr1\t$c[0]\t$c[1]\t+\n";
#	    if($str =~ /hhhhh/) {
#		print "line1 = $line1\n";
#		print "line2 = $line2\n";
#	    }
#	    print "$id\t";
	    print $str;
	}
    } else {
	@b = split(/, /,$spans1);
	for($i=0; $i<@b; $i++) {
	    @c = split(/-/,$b[$i]);
	    $str = "$chr1\t$c[0]\t$c[1]\t+\n";
#	    print "$id\t";
	    print $str;
	}
	@b = split(/, /,$spans2);
	for($i=0; $i<@b; $i++) {
	    @c = split(/-/,$b[$i]);
	    $str = "$chr2\t$c[0]\t$c[1]\t+\n";
#	    print "$id\t";
	    print $str;
	}
    }
}

sub union () {
    ($spans1_u, $spans2_u) = @_;

    undef %chash;
    @a = split(/, /,$spans1_u);
    for($i=0;$i<@a;$i++) {
	@b = split(/-/,$a[$i]);
	for($j=$b[0];$j<=$b[1];$j++) {
	    $chash{$j}++;
	}
    }
    @a = split(/, /,$spans2_u);
    for($i=0;$i<@a;$i++) {
	@b = split(/-/,$a[$i]);
	for($j=$b[0];$j<=$b[1];$j++) {
	    $chash{$j}++;
	}
    }
    $first = 1;
    foreach $pos (sort {$a<=>$b} keys %chash) {
	if($first == 1) {
	    $spans_union = $pos;
	    $first = 0;
	} else {
	    if($pos > $pos_prev + 1) {
		$spans_union = $spans_union . "-$pos_prev, $pos";
	    }
	}
	$pos_prev = $pos;
    }
    $spans_union = $spans_union . "-$pos_prev";
    return $spans_union;
}

sub cig2spans () {
    ($matchstring, $start) = @_;
    $spans_c = "";
    $current_loc = $start;
    $offset = 0;
    while($matchstring =~ /^(\d+)([^\d])/) {
	$num = $1;
	$type = $2;
	$matchstring =~ s/^\d+[^\d]//;
	if($type eq 'M') {
	    $E = $current_loc + $num - 1;
	    if($spans_c =~ /\S/) {
		$spans_c = $spans_c . ", " .  $current_loc . "-" . $E;
	    } else {
		$spans_c = $current_loc . "-" . $E;
	    }
	    $offset = $E - $current_loc + 1;
	    $current_loc = $E;
	}
	if($type eq 'D' || $type eq 'N') {
	    $current_loc = $current_loc + $num + 1;
	}
	if($type eq 'I') {
	    $current_loc++;
	    $offset = $offset  + $num + 1;
	    $offset = $offset + 1;
	}
    }
    $spans2_c = "";
    while($spans2_c ne $spans_c) {
	$spans2_c = $spans_c;
	@b = split(/, /, $spans_c);
	for($i=0; $i<@b-1; $i++) {
	    @c1 = split(/-/, $b[$i]);
	    @c2 = split(/-/, $b[$i+1]);
	    if($c1[1] + 1 >= $c2[0]) {
		$str = "-$c1[1], $c2[0]";
		$spans_c =~ s/$str//;
	    }
	}
    }
    return $spans_c;
}
