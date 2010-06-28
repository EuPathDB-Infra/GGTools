# Written by Gregory R. Grant
# University of Pennsylvania, 2010

if(@ARGV < 1) {

die "
Usage: find_trend.pl <data table> <num perms> <num bins>

Where <data table> has a header row and a first column of ids.

This program searches the profiles (rows of intensities) in the table for
things which are significantly increasing or decreasing, where permutations
of the data obtained by permuting columns are used as a control.  Two
statistics are used, the slope of the least square fit to the data (after
global scaling by the average value across the row), or correlation
with the vector (1, 2, 3, ..., N).

Try <num perms> = 10 and <num bins> = 100

";

}

$| = 1;


$num_bins = $ARGV[2];

open(INFILE, "id2info");
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    $a[1] =~ s/---\s*//g;
    $id2info{$a[0]} = $a[1];
}
close(INFILE);
open(INFILE, $ARGV[0]);
$numperms = $ARGV[1];

$line = <INFILE>;
$row = 0;
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    if(!($id2info{$a[0]} =~ /control/ || $id2info{$a[0]} =~ /AFFX/ || $id2info{$a[0]} =~ /spike/ || $id2info{$a[0]} =~ /ncrna/ || $id2info{$a[0]} =~ /pseudogene/)) {
	$n = @a;
	$ave = 0;
	for($i=1; $i<$n; $i++) {
	    $ave = $ave + $a[$i];
	}
	$ave = $ave /($n-1);
	$names[$row] = $a[0];
	for($i=1; $i<$n; $i++) {
	    $data[$row][$i-1] = $a[$i] / $ave;
	}
	$row++;
    }
}
close(INFILE);

$numrows = $row;
$vect_length = $n-1;
for($i=0; $i<$vect_length; $i++) {
    $X[$i] = $i;
    $Z[$i] = $vect_length - $i;
}

for($i=0; $i<$vect_length; $i++) {
    $Y[$i] = $data[0][$i];
}
$slope = LS_SLOPE(\@X, \@Y);
$max_slope = $slope;
$min_slope = $slope;
for($row = 0; $row<$numrows; $row++) {
    $str0 = "";
    for($i=0; $i<$vect_length; $i++) {
	$Y[$i] = $data[$row][$i];
	$str0 = $str0 . "\t$data[$row][$i]";
    }
    $datastring[$row] = $str0;
    $slope = LS_SLOPE(\@X, \@Y);
    if($slope > $max_slope) {
	$max_slope = $slope;
    }
    if($slope < $min_slope) {
	$min_slope = $slope;
    }

    $rcorr = rank_correlation(\@Z, \@Y);
    if($rcorr > $max_rcorr) {
	$max_rcorr = $rcorr;
    }
    if($rcorr < $min_rcorr) {
	$min_rcorr = $rcorr;
    }
}
print "max_slope = $max_slope\n";
print "min_slope = $min_slope\n\n";
print "max_rcorr = $max_rcorr\n";
print "min_rcorr = $min_rcorr\n\n";

for($i=0; $i<=$num_bins; $i++) {
    $cutoff_up_slope[$i] = $max_slope / $num_bins * $i;
    $cutoff_down_slope[$i] = $min_slope / $num_bins * $i;
}
for($i=0; $i<=$num_bins; $i++) {
    $cutoff_up_rcorr[$i] = $max_rcorr / $num_bins * $i;
    $cutoff_down_rcorr[$i] = $min_rcorr / $num_bins * $i;
}

#for($i=0; $i<=$num_bins; $i++) {
#    print "cutoff_down_slope[$i] = $cutoff_down_slope[$i]\n";
#}
#for($i=0; $i<=$num_bins; $i++) {
#    print "cutoff_down_rcorr[$i] = $cutoff_down_rcorr[$i]\n";
#}

for($i=0; $i<=$num_bins; $i++) {
    $count_up_slope[$i]=0;
    $count_down_slope[$i]=0;
}
for($i=0; $i<=$num_bins; $i++) {
    $count_up_rcorr[$i]=0;
    $count_down_rcorr[$i]=0;
}

for($row = 0; $row<$numrows; $row++) {
    for($i=0; $i<$vect_length; $i++) {
	$Y[$i] = $data[$row][$i];
    }
    $slope = LS_SLOPE(\@X, \@Y);

    if($slope >= 0) {
	$trend[$row] = "up";
	for($i=0; $i<=$num_bins; $i++) {
	    if($slope <= $cutoff_up_slope[$i]) {
		$count_up_slope[$i]++;
		$bin_up_slope[$row] = $i;
		$bin_down_slope[$row] = 0;
		$i = $num_bins+1;
	    }
	}
    }
    if($slope > $cutoff_up_slope[$num_bins]) { # in case loss of precision makes the ineq do something funny
	$count_up_slope[$num_bins]++;
	$bin_up_slope[$row] = $num_bins;
	$bin_down_slope[$row] = 0;
    }    
    if($slope < 0) {
	$trend[$row] = "down";
	for($i=0; $i<=$num_bins; $i++) {
	    if($slope >= $cutoff_down_slope[$i]) {
		$count_down_slope[$i]++;
		$bin_down_slope[$row] = $i;
		$bin_up_slope[$row] = 0;
		$i = $num_bins+1;
	    }
	}
    }
    if($slope < $cutoff_down_slope[$num_bins]) { # in case loss of precision makes the ineq do something funny
	$count_down_slope[$num_bins]++;
	$bin_down_slope[$row] = $num_bins;
	$bin_up_slope[$row] = 0;
    }

    $rcorr = rank_correlation(\@Z, \@Y);
    if($rcorr >= 0) {
	for($i=0; $i<=$num_bins; $i++) {
	    if($rcorr <= $cutoff_up_rcorr[$i]) {
		$count_up_rcorr[$i]++;
		$bin_up_rcorr[$row] = $i;
		$bin_down_rcorr[$row] = 0;
		$i = $num_bins+1;
	    }
	}
    }
    if($rcorr > $cutoff_up_rcorr[$num_bins]) { # in case loss of precision makes the ineq do something funny
	$count_up_rcorr[$num_bins]++;
	$bin_up_rcorr[$row] = $num_bins;
	$bin_down_rcorr[$row] = 0;
    }    
    if($rcorr < 0) {
	for($i=0; $i<=$num_bins; $i++) {
	    if($rcorr >= $cutoff_down_rcorr[$i]) {
		$count_down_rcorr[$i]++;
		$bin_down_rcorr[$row] = $i;
		$bin_up_rcorr[$row] = 0;
		$i = $num_bins+1;
	    }
	}
    }
    if($rcorr < $cutoff_down_rcorr[$num_bins]) { # in case loss of precision makes the ineq do something funny
	$count_down_rcorr[$num_bins]++;
	$bin_down_rcorr[$row] = $num_bins;
	$bin_up_rcorr[$row] = 0;
    }

}

$cumulative_count_down_slope_unpermunted[$num_bins] = $count_down_slope[$num_bins];
$cumulative_count_up_slope_unpermunted[$num_bins] = $count_up_slope[$num_bins];
for($i=$num_bins-1; $i>=0; $i--) {
    $cumulative_count_down_slope_unpermunted[$i] = $cumulative_count_down_slope_unpermunted[$i+1] + $count_down_slope[$i];
    $cumulative_count_up_slope_unpermunted[$i] = $cumulative_count_up_slope_unpermunted[$i+1] + $count_up_slope[$i];
}
#for($i=0; $i<=$num_bins; $i++) {
#    print "cumulative_count_down_slope_unpermunted[$i] = $cumulative_count_down_slope_unpermunted[$i]\n";
#}
#print "\n";
#for($i=0; $i<=$num_bins; $i++) {
#    print "cumulative_count_up_slope_unpermunted[$i] = $cumulative_count_up_slope_unpermunted[$i]\n";
#}
#print "\n";

$cumulative_count_down_rcorr_unpermunted[$num_bins] = $count_down_rcorr[$num_bins];
$cumulative_count_up_rcorr_unpermunted[$num_bins] = $count_up_rcorr[$num_bins];
for($i=$num_bins-1; $i>=0; $i--) {
    $cumulative_count_down_rcorr_unpermunted[$i] = $cumulative_count_down_rcorr_unpermunted[$i+1] + $count_down_rcorr[$i];
    $cumulative_count_up_rcorr_unpermunted[$i] = $cumulative_count_up_rcorr_unpermunted[$i+1] + $count_up_rcorr[$i];
}
#for($i=0; $i<=$num_bins; $i++) {
#    print "cumulative_count_down_rcorr_unpermunted[$i] = $cumulative_count_down_rcorr_unpermunted[$i]\n";
#}
#print "\n";
#for($i=0; $i<=$num_bins; $i++) {
#    print "cumulative_count_up_rcorr_unpermunted[$i] = $cumulative_count_up_rcorr_unpermunted[$i]\n";
#}
#print "\n";

for($perm=0; $perm<$numperms; $perm++) {
#    print "\nperm = $perm\n";
    for($i=0; $i<$vect_length; $i++) {
	$indexes[$i] = $i;
    }
    for($i=0; $i<=$num_bins; $i++) {
	$count_up_slope[$i]=0;
	$count_down_slope[$i]=0;
    }
    for($i=0; $i<=$num_bins; $i++) {
	$count_up_rcorr[$i]=0;
	$count_down_rcorr[$i]=0;
    }
    if($perm > 0) {  # causes the unpermuted to be the first perm
	for($j=0; $j<$vect_length * 100; $j++) {  # this should be enough swappings to randomize it
	    $index1 = int(rand($vect_length));
	    $index2 = $index1;
	    while($index1 == $index2) {
		$index2 = int(rand($vect_length));
	    }
	    $temp = $indexes[$index1];
	    $indexes[$index1] = $indexes[$index2];
	    $indexes[$index2] = $temp;
	}
    }
    for($row=0; $row<$numrows; $row++) {
	for($i=0; $i<$vect_length; $i++) {
	    $Y[$i] = $data[$row][$indexes[$i]]; # @Y gets the randomized row
	}
	$slope = LS_SLOPE(\@X, \@Y);
	if($slope >= 0) {
	    for($i=0; $i<=$num_bins; $i++) {
		if($slope <= $cutoff_up_slope[$i]) {
		    $count_up_slope[$i]++;
		    $i = $num_bins+1;
		}
	    }
	}
	if($slope > $cutoff_up_slope[$num_bins]) { # in case loss of precision makes the ineq do something funny
	    $count_up_slope[$num_bins]++;
	}    
	if($slope < 0) {
	    for($i=0; $i<=$num_bins; $i++) {
		if($slope >= $cutoff_down_slope[$i]) {
		    $count_down_slope[$i]++;
		    $i = $num_bins+1;
		}
	    }
	}
	if($slope < $cutoff_down_slope[$num_bins]) { # in case loss of precision makes the ineq do something funny
	    $count_down_slope[$num_bins]++;
	}

	$rcorr = rank_correlation(\@Z, \@Y);
	if($rcorr >= 0) {
	    for($i=0; $i<=$num_bins; $i++) {
		if($rcorr <= $cutoff_up_rcorr[$i]) {
		    $count_up_rcorr[$i]++;
		    $i = $num_bins+1;
		}
	    }
	}
	if($rcorr > $cutoff_up_rcorr[$num_bins]) { # in case loss of precision makes the ineq do something funny
	    $count_up_rcorr[$num_bins]++;
	}    
	if($rcorr < 0) {
	    for($i=0; $i<=$num_bins; $i++) {
		if($rcorr >= $cutoff_down_rcorr[$i]) {
		    $count_down_rcorr[$i]++;
		    $i = $num_bins+1;
		}
	    }
	}
	if($rcorr < $cutoff_down_rcorr[$num_bins]) { # in case loss of precision makes the ineq do something funny
	    $count_down_rcorr[$num_bins]++;
	}

    }
    $cumulative_count_down_slope[$num_bins] = $count_down_slope[$num_bins];
    $cumulative_count_up_slope[$num_bins] = $count_up_slope[$num_bins];
    for($i=$num_bins-1; $i>=0; $i--) {
	$cumulative_count_down_slope[$i] = $cumulative_count_down_slope[$i+1] + $count_down_slope[$i];
	$cumulative_count_up_slope[$i] = $cumulative_count_up_slope[$i+1] + $count_up_slope[$i];
    }
    for($i=0; $i<=$num_bins; $i++) {
	$cumulative_count_down_slope_ave[$i] = $cumulative_count_down_slope_ave[$i] + $cumulative_count_down_slope[$i];
	$cumulative_count_up_slope_ave[$i] = $cumulative_count_up_slope_ave[$i] + $cumulative_count_up_slope[$i];
    }

    $cumulative_count_down_rcorr[$num_bins] = $count_down_rcorr[$num_bins];
    $cumulative_count_up_rcorr[$num_bins] = $count_up_rcorr[$num_bins];
    for($i=$num_bins-1; $i>=0; $i--) {
	$cumulative_count_down_rcorr[$i] = $cumulative_count_down_rcorr[$i+1] + $count_down_rcorr[$i];
	$cumulative_count_up_rcorr[$i] = $cumulative_count_up_rcorr[$i+1] + $count_up_rcorr[$i];
    }
    for($i=0; $i<=$num_bins; $i++) {
	$cumulative_count_down_rcorr_ave[$i] = $cumulative_count_down_rcorr_ave[$i] + $cumulative_count_down_rcorr[$i];
	$cumulative_count_up_rcorr_ave[$i] = $cumulative_count_up_rcorr_ave[$i] + $cumulative_count_up_rcorr[$i];
    }

}
for($i=0; $i<=$num_bins; $i++) {
    $cumulative_count_down_slope_ave[$i] = $cumulative_count_down_slope_ave[$i] / $numperms;
    $cumulative_count_up_slope_ave[$i] = $cumulative_count_up_slope_ave[$i] / $numperms;
}
for($i=0; $i<=$num_bins; $i++) {
    $cumulative_count_down_rcorr_ave[$i] = $cumulative_count_down_rcorr_ave[$i] / $numperms;
    $cumulative_count_up_rcorr_ave[$i] = $cumulative_count_up_rcorr_ave[$i] / $numperms;
}

#for($i=0; $i<=$num_bins; $i++) {
#    print "cumulative_count_down_slope_ave[$i] = $cumulative_count_down_slope_ave[$i]\n";
#}
#print "\n";
#for($i=0; $i<=$num_bins; $i++) {
#    print "cumulative_count_up_slope_ave[$i] = $cumulative_count_up_slope_ave[$i]\n";
#}
#print "\n";

#for($i=0; $i<=$num_bins; $i++) {
#    print "cumulative_count_down_rcorr_ave[$i] = $cumulative_count_down_rcorr_ave[$i]\n";
#}
#print "\n";
#for($i=0; $i<=$num_bins; $i++) {
#    print "cumulative_count_up_rcorr_ave[$i] = $cumulative_count_up_rcorr_ave[$i]\n";
#}
#print "\n";

$FDR_down_slope[0] = 1;
$FDR_up_slope[0] = 1;
$FDR_down_rcorr[0] = 1;
$FDR_up_rcorr[0] = 1;

for($i=1; $i<=$num_bins; $i++) {
    if($cumulative_count_down_slope_unpermunted[$i] > 0) {
	$FDR_down_slope[$i] = $cumulative_count_down_slope_ave[$i] / $cumulative_count_down_slope_unpermunted[$i];
    }
    else {
	$FDR_down_slope[$i] = 1;
    }
    if($cumulative_count_up_slope_unpermunted[$i] > 0) {
	$FDR_up_slope[$i] = $cumulative_count_up_slope_ave[$i] / $cumulative_count_up_slope_unpermunted[$i];
    }
    else {
	$FDR_up_slope[$i] = 1;
    }
    if($FDR_up_slope[$i] > 1) {
	$FDR_up_slope[$i] = 1;
    }
    if($FDR_down_slope[$i] > 1) {
	$FDR_down_slope[$i] = 1;
    }
}
for($i=1; $i<=$num_bins; $i++) {
    if($cumulative_count_down_rcorr_unpermunted[$i] > 0) {
	$FDR_down_rcorr[$i] = $cumulative_count_down_rcorr_ave[$i] / $cumulative_count_down_rcorr_unpermunted[$i];
    }
    else {
	$FDR_down_rcorr[$i] = 1;
    }
    if($cumulative_count_up_rcorr_unpermunted[$i] > 0) {
	$FDR_up_rcorr[$i] = $cumulative_count_up_rcorr_ave[$i] / $cumulative_count_up_rcorr_unpermunted[$i];
    }
    else {
	$FDR_up_rcorr[$i] = 1;
    }
    if($FDR_up_rcorr[$i] > 1) {
	$FDR_up_rcorr[$i] = 1;
    }
    if($FDR_down_rcorr[$i] > 1) {
	$FDR_down_rcorr[$i] = 1;
    }
}

$FDR_down_slope_monotonized[0] = $FDR_down_slope[0];
$m = $FDR_down_slope[0];
for($i=1; $i<=$num_bins; $i++) {
    if($FDR_down_slope[$i] >= $m) {
	$FDR_down_slope_monotonized[$i] = $m;
    }
    else {
	$FDR_down_slope_monotonized[$i] = $FDR_down_slope[$i];
	$m = $FDR_down_slope[$i];
    }
}
$FDR_up_slope_monotonized[0] = $FDR_up_slope[0];
$m = $FDR_up_slope[0];
for($i=1; $i<=$num_bins; $i++) {
    if($FDR_up_slope[$i] >= $m) {
	$FDR_up_slope_monotonized[$i] = $m;
    }
    else {
	$FDR_up_slope_monotonized[$i] = $FDR_up_slope[$i];
	$m = $FDR_up_slope[$i];
    }
}

$FDR_down_rcorr_monotonized[0] = $FDR_down_rcorr[0];
$m = $FDR_down_rcorr[0];
for($i=1; $i<=$num_bins; $i++) {
    if($FDR_down_rcorr[$i] >= $m) {
	$FDR_down_rcorr_monotonized[$i] = $m;
    }
    else {
	$FDR_down_rcorr_monotonized[$i] = $FDR_down_rcorr[$i];
	$m = $FDR_down_rcorr[$i];
    }
}
$FDR_up_rcorr_monotonized[0] = $FDR_up_rcorr[0];
$m = $FDR_up_rcorr[0];
for($i=1; $i<=$num_bins; $i++) {
    if($FDR_up_rcorr[$i] >= $m) {
	$FDR_up_rcorr_monotonized[$i] = $m;
    }
    else {
	$FDR_up_rcorr_monotonized[$i] = $FDR_up_rcorr[$i];
	$m = $FDR_up_rcorr[$i];
    }
}

#for($i=0; $i<=$num_bins; $i++) {
#    print("FDR_down_slope[$i] = $FDR_down_slope[$i]\n");
#}
#for($i=0; $i<=$num_bins; $i++) {
#    print("FDR_up_slope[$i] = $FDR_up_slope[$i]\n");
#}

#for($i=0; $i<=$num_bins; $i++) {
#    print("FDR_down_rcorr[$i] = $FDR_down_rcorr[$i]\n");
#}
#for($i=0; $i<=$num_bins; $i++) {
#    print("FDR_up_rcorr[$i] = $FDR_up_rcorr[$i]\n");
#}

#for($i=0; $i<=$num_bins; $i++) {
#    print("FDR_down_slope_monotonized[$i] = $FDR_down_slope_monotonized[$i]\n");
#}
#for($i=0; $i<=$num_bins; $i++) {
#    print("FDR_up_slope_monotonized[$i] = $FDR_up_slope_monotonized[$i]\n");
#}

#for($i=0; $i<=$num_bins; $i++) {
#    print("FDR_down_rcorr_monotonized[$i] = $FDR_down_rcorr_monotonized[$i]\n");
#}
#for($i=0; $i<=$num_bins; $i++) {
#    print("FDR_up_rcorr_monotonized[$i] = $FDR_up_rcorr_monotonized[$i]\n");
#}

#print "ID\tFDR\tdirection\tFDR slope up\tFDR slope down\tFDR rcorr up\tFDR rcorr down\tinfo\tdata\n";
#print "ID\tb_up_slope\tFDR_up_slope\tb_down_slope\tFDR_down_slope\tb_up_rcorr\tFDR_up_rcorr\tb_down_corr\tFDR_down_rcorr\n";
print "ID\tFDR\tdirection\ttype\tinfo\tdata\n";
for($row=0; $row<$numrows; $row++) {
    $b_up_slope = $bin_up_slope[$row];
    $b_down_slope = $bin_down_slope[$row];
    $b_up_rcorr = $bin_up_rcorr[$row];
    $b_down_rcorr = $bin_down_rcorr[$row];
    $F_up_slope = $FDR_up_slope_monotonized[$b_up_slope];
    $min_f = $F_up_slope;
    $min_at = "slope";
    $F_down_slope = $FDR_down_slope_monotonized[$b_down_slope];
    if($F_down_slope < $min_f) {
	$min_f = $F_down_slope;
    }
    $F_up_rcorr = $FDR_up_rcorr_monotonized[$b_up_rcorr];
    if($F_up_rcorr < $min_f) {
	$min_f = $F_up_rcorr;
	$min_at = "trend";
    }
    $F_down_rcorr = $FDR_down_rcorr_monotonized[$b_down_rcorr];
    if($F_down_rcorr < $min_f) {
	$min_f = $F_down_rcorr;
	$min_at = "trend";
    }
    $ds = $datastring[$row];
    $info = $id2info{$names[$row]};
    
#    print "$names[$row]\t$min_f\t$trend[$row]\t$min_at\t$F_up_slope\t$F_down_slope\t$F_up_rcorr\t$F_down_rcorr\t$info$ds\n";
    print "$names[$row]\t$min_f\t$trend[$row]\t$min_at\t$info$ds\n";
}

sub LS_SLOPE {
    my ($X_vect_ref, $Y_vect_ref) = @_;
    my @X_vect = @{$X_vect_ref};
    my @Y_vect = @{$Y_vect_ref};

    $X_vect_length = @X_vect;
    $Y_vect_length = @Y_vect;
    
    if($X_vect_length != $Y_vect_length) {
	print("ERROR: tried to run the LS subroutine with vectors of differing lengths\n");
	exit("");
    }
    
    $y=0;
    $x=0;
    $yy=0;
    $xx=0;
    $xy=0;
    
    for($i=0; $i<$X_vect_length; $i++) {
	$x=$x+$X_vect[$i];
	$y=$y+$Y_vect[$i];
	$xx=$xx+($X_vect[$i])*($X_vect[$i]);
	$yy=$yy+($Y_vect[$i])*($Y_vect[$i]);
	$xy=$xy+($X_vect[$i])*($Y_vect[$i]);
    }
    
    $slope =($X_vect_length*$xy-$x*$y)/($X_vect_length*$xx-$x*$x);  # slope

    return $slope;
}

sub rank_correlation {
    ($vector1_ref, $vector2_ref) = @_;
    @vector1 = @{$vector1_ref};
    @vector2 = @{$vector2_ref};
    $length_vector1=@vector1;
    $length_vector2=@vector2;
    @ranks1 = ranks(@vector1);
    @ranks2 = ranks(@vector2);
    for($i=0; $i<$length_vector1; $i++) {
	$ranks1[$i] = $length_vector1 - $ranks1[$i] + 1;
    }
    $r = int(correlation(\@ranks1, \@ranks2) * 1000) / 1000;

    return $r;
}


sub ranks {
    @positions = sort {$_[$b]<=>$_[$a]} (0 .. $#_);
    @ranks = sort {$positions[$a]<=>$positions[$b]} (0 .. $#_);
    map {$_+1} @ranks;
}

sub correlation {
    ($vector1, $vector2) = @_;
    $length_vector1=@{$vector1};
    $length_vector2=@{$vector2};
    $i;
    $j;
    $m;
    $mean1;
    $mean2;
    @values1;
    @values2;
    $length_values;
    $S;
    $t;
    $sd1;
    $sd2;

    if($length_vector1 != $length_vector2) {
	die "Error: tried to correlate vectors of different lengths\n";
    }
    $length_vector = $length_vector1;
    $j=0;
    $m=0;
    $eq = 0;
    for($i=0;$i<$length_vector;$i++) {
	if(((defined $vector1->[$i]) && ($vector1->[$i] != "")) && !($vector1->[$i] =~ /[a-zA-Z]/)) {
	    if(((defined $vector2->[$i]) && ($vector2->[$i] != "")) && !($vector2->[$i] =~ /[a-zA-Z]/)) {
		$values1[$j] = $vector1->[$i];
		$m1=$m1+$vector1->[$i];
		$values2[$j] = $vector2->[$i];
		$m2=$m2+$vector2->[$i];
		$eq = $eq + ($values1[$j] - $values2[$j])**2;
		$j++;
	    }
	}
    }

    if($eq == 0) {
	return 1;
    }

    $length_values = $j;


    $x2 = 0;
    $y2 = 0;
    $x=0;
    $y=0;
    $xy=0;
    for($i=0; $i<$length_values; $i++) {
	$x = $x + $values1[$i];
	$x2 = $x2 + $values1[$i]**2;
	$y = $y + $values2[$i];
	$y2 = $y2 + $values2[$i]**2;
	$xy = $xy + $values1[$i] * $values2[$i];
    }

    $ssx = $x2 - $x*$x/$length_values;
    $ssy = $y2 - $y*$y/$length_values;
    $ssxy = $xy - $x*$y/$length_values;

    $r = $ssxy / sqrt($ssx * $ssy);

    return $r;
}    
