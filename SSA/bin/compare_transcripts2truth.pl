
if(@ARGV<2) {
    die "
Usage: compare_transcripts2truth.pl <truth> <inferred>

Where:  <truth> is the bedfile of the true transcripts
        <inferred> is the bedfile of the inferred transcripts

";
}

$truth = $ARGV[0];
$inferred = $ARGV[1];

print STDERR "reading the truth file\n";
open(INFILE, $truth);
$line = <INFILE>;
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    $chr = $a[0];
    $START = $a[1];
    $blocksizes = $a[10];
    $blocksizes =~ s/\s*,\s*$//;
    @B = split(/,/,$blocksizes);
    $offsets = $a[11];
    $offsets =~ s/\s*,\s*$//;
    @O = split(/,/,$offsets);
    $starts = "";
    $ends = "";
    for($i=0;$i<@O;$i++) {
	$start = $START + $O[$i];
	$starts = $starts . "$start,";
	$end = $start + $B[$i];
	$ends = $ends . "$end,";
    }
    $starts =~ s/\s*,\s*$//;
    $ends =~ s/\s*,\s*$//;
    $structure = $chr . ":" . $starts . ":" . $ends;
    $truth{$structure}++;
#    print "structure = $structure\n";
    if($starts =~ /,/) {
	$starts =~ s/^\d+,//;
	$ends =~ s/,\d+$//;
	$structure_without_ends = $chr . ":" . $starts . ":" . $ends;
#	print "structure_without_ends=$structure_without_ends\n";
	$truth_without_ends{$structure_without_ends}++;
    }
}
close(INFILE);

$num_true=0;
$num_false=0;
open(INFILE, $inferred);
$line = <INFILE>;
print STDERR "reading the inferred file\n";
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    $chr = $a[0];
    $START = $a[1];
    $blocksizes = $a[10];
    $blocksizes =~ s/\s*,\s*$//;
    @B = split(/,/,$blocksizes);
    $offsets = $a[11];
    $offsets =~ s/\s*,\s*$//;
    @O = split(/,/,$offsets);
    $starts = "";
    $ends = "";
    for($i=0;$i<@O;$i++) {
	$start = $START + $O[$i];
	$starts = $starts . "$start,";
	$end = $start + $B[$i];
	$ends = $ends . "$end,";
    }
    $starts =~ s/\s*,\s*$//;
    $ends =~ s/\s*,\s*$//;

    @S = split(/,/, $starts);
    @E = split(/,/, $ends);
    $start_string = $S[0] . ",";
    $end_string = "";
    $N = @S;
    for($i=1; $i<$N; $i++) {
        $intronlength = $S[$i] - $E[$i-1];
        $realstart = $E[$i-1] + 1;
        $realend = $S[$i];
        $length = $realend - $realstart + 1;
        if($length >= 10) {
            $start_string = $start_string . $S[$i] . ",";
            $end_string = $end_string . $E[$i-1] . ",";
        }
        else {
#	    print "$line\n";
        }
    }
    $end_string = $end_string . $E[$N-1] . ",";
    $start_string =~ s/,\s*$//;
    $end_string =~ s/,\s*$//;

    $structure = $chr . ":" . $start_string . ":" . $end_string;
#    $structure = $chr . ":" . $starts . ":" . $ends;
#    print "structure = $structure\n";
    $inferred{$structure}++;
    if(defined $truth{$structure}) {
	$num_true++;
	print "$structure\n";
    } else {
	$num_false++;
	if($starts =~ /,/) {
	    $starts =~ s/^\d+,//;
	    $ends =~ s/,\d+$//;
	    $structure_without_ends = $chr . ":" . $starts . ":" . $ends;
#	    print "structure_without_ends=$structure_without_ends\n";
	    $inferred_without_ends{$structure_without_ends}++;
	    if(defined $truth_without_ends{$structure_without_ends}) {
		$num_true_without_ends++;
	    }
	}
    }
}
close(INFILE);
$total = $num_true + $num_false;
$percent = int($num_true / $total * 10000) / 100;;
print "True Positive Rate: $percent ($num_true)\n";
$percent = int($num_false / $total * 10000) / 100;;
print "False Positive Rate: $percent ($num_false)\n";

$num_true_extra = $num_true + $num_true_without_ends;
$percent = int($num_true_extra / $total * 10000) / 100;;
print "True Positive Rate Without Ends: $percent ($num_true_extra)\n";


