#!/usr/bin/perl

# Written by Gregory R. Grant
# Unviersity of Pennsylvania, 2010

if(@ARGV<1) {
    die "
Usage: table2columnwisepercents.pl <table> [options]

Where <table> is a table of rows and coluns of numeric values.
Table is expected to have a first row of column names and a first column
of row names.

Note: this script is not very smart about ties, but for a large table of
continuous data that shouldn't be a problem.

Options:

    -sort n  : sort on column n

";
    
}

# Table assumbed to have header and 1st column of row ids

$sort = "false";
for($i=1; $i<@ARGV; $i++) {
    $optionrecognized = 0;
    if($ARGV[$i] eq "-sort") {
	$sort = "true";
	$i++;
	$COLUMN4SORT = $ARGV[$i];
	if(!($COLUMN4SORT =~ /^\d+$/) || $COLUMN4SORT == 0) {
	    print "\nERROR: $COLUMN4SORT is not a valid column number to sort on.\n\n";
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
    for($row=0; $row<$numrows; $row++) {
	$str = $ids[$row];
	for($col=1; $col<$numcols; $col++) {
	    $str = $str . "\t$percents[$row][$col]";
	}
	$str = $str . "\n";
	$thing{$str} = $percents[$row][$COLUMN4SORT];
    }
    foreach $key (sort {$thing{$a} <=> $thing{$b}} keys %thing) {
	print $key;
    }
}
