#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

if(@ARGV < 2) {
    die "
Usage compare_covs.pl <cov1> <cov2> [options]

Where: <cov1> and <cov2> are two coverage files.

Options: -open : if the intervals in the coverage files do not contain the right endpoint.

";
}

$open = "false";
for($i=2; $i<@ARGV; $i++) {
    $optionrecognized = 0;
    if($ARGV[$i] eq "-open") {
	$open = "true";
	$optionrecognized = 1;
    }
    if($optionrecognized == 0) {
	print "\n\nERROR: Option $ARGV[$i] not recognized.\n\n";
	exit();
    }
}

if($open eq "true") {
    $adjust = 1;
}
else {
    $adjust = 0;
}


$total1 = 0;
$numpositions = 0;
$numpositions1 = 0;
$numpositions2 = 0;
$intersection = 0;

open(INFILE, $ARGV[0]);
while($line = <INFILE>) {
    chomp($line);
    if($line =~ /^(\S+)\t(\d+)\t(\d+)\t(\d+)/) {
	$chr = $1;
	$total1 = $total1 + ($3 - $2 + 1 - $adjust) * $4;
	$chrs{$chr}++;
	$allchrs{$chr}++;
	if($chrs{$chr} == 1) {
	    $c = $chr;
	    $c =~ s/[^a-zA-Z0-9_]//g;
	    open $F1{$chr}, ">" . $c . ".1";
	    $FF = $F1{$chr};
	    print $FF "$line\n";
	}
	else {
	    $FF = $F1{$chr};
	    print $FF "$line\n";
	}
    }
}
close(INFILE);
foreach $chr (keys %F1) {
    close $F1{$chr};
}

if($total1 < 10000) {
    print "total file 1: $total1\n";
} else {
    $f = format_large_int($total1);
    print "total file 1: $f\n";
}

undef %chrs;
open(INFILE, $ARGV[1]);
$total2 =0;
while($line = <INFILE>) {
    chomp($line);
    if($line =~ /^(\S+)\t(\d+)\t(\d+)\t(\d+)/) {
	$chr = $1;
	$total2 = $total2 + ($3 - $2 + 1 - $adjust) * $4;
	$chrs{$chr}++;
	$allchrs{$chr}++;
	if($chrs{$chr} == 1) {
	    $c = $chr;
	    $c =~ s/[^a-zA-Z0-9_]//g;
	    open $F2{$chr}, ">" . $c . ".2";
	    $FF = $F2{$chr};
	    print $FF "$line\n";
	}
	else {
	    $FF = $F2{$chr};
	    print $FF "$line\n";
	}
    }
}
close(INFILE);
foreach $chr (keys %F2) {
    close $F2{$chr};
}
$f = format_large_int($total2);
print "total file 2: $f\n";

$diff = 0;
foreach $CHR (keys %allchrs) {
    undef %cov1;
    $c = $CHR;
    $c =~ s/[^a-zA-Z0-9_]//g;
    $c1 = $c . ".1";
    $c2 = $c . ".2";
    open(INFILE1, $c1);
    while($line = <INFILE1>) {
	chomp($line);
	if($line =~ /^(\S+)\t(\d+)\t(\d+)\t(\d+)/) {
	    $chr = $1;
	    $start = $2;
	    $end = $3;
	    $cov = $4;
	    if($chr eq $CHR) {
		for($i=$start; $i<=$end-$adjust; $i++) {
		    $cov1{$i}=$cov;
		    $numpositions1++;
		}
	    }
	}
    }
    close(INFILE1);
    
    open(INFILE2, $c2);
    while($line = <INFILE2>) {
	chomp($line);
	if($line =~ /^(\S+)\t(\d+)\t(\d+)\t(\d+)/) {
	    $chr = $1;
	    $start = $2;
	    $end = $3;
	    $cov = $4;
	    if($chr eq $CHR) {
		for($i=$start; $i<=$end-$adjust; $i++) {
		    if(exists $cov1{$i}) {
			$absval = absvalue($cov1{$i}, $cov);
			$diff = $diff + $absval;
			if($cov1{$i} >= 1) {
			    $CV1 = $CV1 + $cov/$cov1{$i};
			    $NCV1++;
			}
			if($cov1{$i} >= 5) {
			    $CV5 = $CV5 + $cov/$cov1{$i};
			    $NCV5++;
			}
			if($cov1{$i} >= 10) {
			    $CV10 = $CV10 + $cov/$cov1{$i};
			    $NCV10++;
			}
			delete $cov1{$i};
			$intersection++;
		    }
		    else {
			$diff = $diff + $cov;		
		    }
		    $numpositions2++;
		    $numpositions++;
		}
	    }
	}
    }

    close(INFILE2);
    foreach $pos (keys %cov1) {
	$diff = $diff + $cov1{$pos};
	$numpositions++;
	if($cov1{$pos} >= 1) {
	    $NCV1++;
	}
	if($cov1{$pos} >= 5) {
	    $NCV5++;
	}
	if($cov1{$pos} >= 10) {
	    $NCV10++;
	}
    }
    if(-e $c1) {
	unlink($c1);
    }
    if(-e $c2) {
	unlink($c2);
    }
}

$f = format_large_int($diff);
print "Total Difference: $f\n";

$percent_diff = $diff / ($total1 + $total2);  # $total1 + $total2 is the max possible difference

print "percent diff: $percent_diff\n";

$ave_diff = $diff / $numpositions;

if($NCV1 > 0) {
    $ave_CV1 = $CV1 / $NCV1;
} else {
    $ave_CV1 = 1;
}
if($NCV5 > 0) {
    $ave_CV5 = $CV5 / $NCV5;
} else {
    $ave_CV5 = 1;
}
if($NCV10 > 0) {
    $ave_CV10 = $CV10 / $NCV10;
} else {
    $ave_CV10 = 1;
}

#print "ave diff: $ave_diff\n";
print "ave ratio cov2 to cov1 for all locations where cov1 >= 1: $ave_CV1\n";
print "ave ratio cov2 to cov1 for all locations where cov1 >= 5: $ave_CV5\n";
print "ave ratio cov2 to cov1 for all locations where cov1 >= 10: $ave_CV10\n";
$ave_depth1 = $total1 / $numpositions1;
$ave_depth2 = $total2 / $numpositions2;
$f = format_large_int($numpositions);
print "combined support: $f\n";
$f = format_large_int($numpositions1);
print "support1 = $f\n";
$f = format_large_int($numpositions2);
print "support2 = $f\n";
$f = format_large_int($intersection);
print "intersection = $f\n";
#print "ave depth 1: $ave_depth1\n";
#print "ave depth 2: $ave_depth2\n";

sub absvalue () {
    ($x,$y) = @_;
    if($x >= $y) {
	$z = $x - $y;
	return $z;
    }
    else {
	$z = $y - $x;
	return $z;
    }
}

sub format_large_int () {
    ($int) = @_;
    @a = split(//,"$int");
    $j=0;
    $newint = "";
    $n = @a;
    for($i=$n-1;$i>=0;$i--) {
	$j++;
	$newint = $a[$i] . $newint;
	if($j % 3 == 0) {
	    $newint = "," . $newint;
	}
    }
    $newint =~ s/^,//;
    return $newint;
}
