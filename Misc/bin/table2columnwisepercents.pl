#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

if(@ARGV<1 || $ARGV[0] eq "help") {
    die "
Usage: table2columnwisepercents.pl <table> [options]

<table> is assumed to have header of column names and first column of row names.

This script replaces each value by it's percentile in it's column.

Note: this script is not very smart about ties, but for a large table of
continuous data that shouldn't be a problem.

Options:  -sort  n : sort on column n
          -sort2 n : sort on column n to break ties when -sort specified.
          -sort3 n : sort on column n to break ties when -sort and -sort2 specified.

By default sorting is in increasing order.  Make the column for -sort negative to sort
in decreasing ordered (make -sort negative only, not sort2, sort3, etc...).

";
}


$sort = "false";
$sort2 = "false";
$sort3 = "false";
$sortdecreasing = "false";
for($i=1; $i<@ARGV; $i++) {
    $optionrecognized = 0;
    if($ARGV[$i] eq "-sort") {
	$sort = "true";
	$i++;
	$COLUMN4SORT = $ARGV[$i];
	if($COLUMN4SORT =~ /^-\d+$/) {
	    $sortdecreasing = "true";
	    $COLUMN4SORT =~ s/^-//;
	}
	if(!($COLUMN4SORT =~ /^\d+$/) || $COLUMN4SORT == 0) {
	    print "\nERROR: $COLUMN4SORT is not a valid column number to sort on.\n\n";
	    exit();
	}
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-sort2") {
	$sort2 = "true";
	$i++;
	$COLUMN4SORT2 = $ARGV[$i];
	if(!($COLUMN4SORT2 =~ /^\d+$/) || $COLUMN4SORT2 == 0) {
	    print "\nERROR: $COLUMN4SORT2 is not a valid column number to sort on.\n\n";
	    exit();
	}
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-sort3") {
	$sort3 = "true";
	$i++;
	$COLUMN4SORT3 = $ARGV[$i];
	if(!($COLUMN4SORT3 =~ /^\d+$/) || $COLUMN4SORT3 == 0) {
	    print "\nERROR: $COLUMN4SORT3 is not a valid column number to sort on.\n\n";
	    exit();
	}
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-sort4") {
	$sort4 = "true";
	$i++;
	$COLUMN4SORT4 = $ARGV[$i];
	if(!($COLUMN4SORT4 =~ /^\d+$/) || $COLUMN4SORT4 == 0) {
	    print "\nERROR: $COLUMN4SORT4 is not a valid column number to sort on.\n\n";
	    exit();
	}
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-sort5") {
	$sort5 = "true";
	$i++;
	$COLUMN4SORT5 = $ARGV[$i];
	if(!($COLUMN4SORT5 =~ /^\d+$/) || $COLUMN4SORT5 == 0) {
	    print "\nERROR: $COLUMN4SORT5 is not a valid column number to sort on.\n\n";
	    exit();
	}
	$optionrecognized = 1;
    }
    if($optionrecognized == 0) {
	print "\nERROR: option $ARGV[$i] not recognized\n";
	exit();
    }
}

open(INFILE, $ARGV[0]);
$line = <INFILE>;
print $line;
$row=0;
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    $numcols = @a;      # this count includes the id column
    $ids[$row] = $a[0];
    for($col=1;$col<$numcols;$col++) {
	$value[$col]{$row} = $a[$col];
    }
    $row++;
}
$numrows = $row;
for($col=1;$col<$numcols;$col++) {
    $c=1;
    foreach $row (sort {$value[$col]{$a}<=>$value[$col]{$b}} keys %{$value[$col]}) {
	if($value[$col]{$row} > 0) {
	    $percents[$row][$col] = int($c/$numrows * 1000) / 1000;
	}
	else {
	    $percents[$row][$col] = 0;
	}
	$c++;
    }
}
if($sort eq "false") {
    for($row=0; $row<$numrows; $row++) {
	print $ids[$row];
	for($col=1; $col<$numcols; $col++) {
	    print "\t$percents[$row][$col]";
	}
	print "\n";
    }
} else {
    if($sort2 eq "false") {
	for($row=0; $row<$numrows; $row++) {
	    $str = $ids[$row];
	    for($col=1; $col<$numcols; $col++) {
		$str = $str . "\t$percents[$row][$col]";
	    }
	    $str = $str . "\n";
	    $thing{$str} = $percents[$row][$COLUMN4SORT];
	}
	if($sortdecreasing eq "false") {
	    foreach $key (sort {$thing{$a} <=> $thing{$b}} keys %thing) {
		print $key;
	    }
	} else {
	    foreach $key (sort {$thing{$b} <=> $thing{$a}} keys %thing) {
		print $key;
	    }
	}
    } else {
	if($sort3 eq "false") {
	    for($row=0; $row<$numrows; $row++) {
		$str = $ids[$row];
		for($col=1; $col<$numcols; $col++) {
		    $str = $str . "\t$percents[$row][$col]";
		}
		$str = $str . "\n";
		$thing{$str} = $percents[$row][$COLUMN4SORT];
		$thing2{$str} = $percents[$row][$COLUMN4SORT2];
	    }
	    if($sortdecreasing eq "false") {
		foreach $key (sort {$thing{$a} <=> $thing{$b} || $thing2{$a} <=> $thing2{$b}} keys %thing) {
		    print $key;
		}
	    } else {
		foreach $key (sort {$thing{$b} <=> $thing{$a} || $thing2{$b} <=> $thing2{$a}} keys %thing) {
		    print $key;
		}
	    }
	} else {
	    for($row=0; $row<$numrows; $row++) {
		$str = $ids[$row];
		for($col=1; $col<$numcols; $col++) {
		    $str = $str . "\t$percents[$row][$col]";
		}
		$str = $str . "\n";
		$thing{$str} = $percents[$row][$COLUMN4SORT];
		$thing2{$str} = $percents[$row][$COLUMN4SORT2];
		$thing3{$str} = $percents[$row][$COLUMN4SORT3];
	    }
	    if($sortdecreasing eq "false") {
		foreach $key (sort {$thing{$a} <=> $thing{$b} || $thing2{$a} <=> $thing2{$b} || $thing3{$a} <=> $thing3{$b}} keys %thing) {
		    print $key;
		}
	    } else {
		foreach $key (sort {$thing{$b} <=> $thing{$a} || $thing2{$b} <=> $thing2{$a} || $thing3{$b} <=> $thing3{$a}} keys %thing) {
		    print $key;
		    
		}
	    }
	}
    }
}
