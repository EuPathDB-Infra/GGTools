#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennslyvania, 2010

$|=1;
if(@ARGV<5) {
    die "
Usage: compare2simulated.pl <working_dir> <simulated_readsfile_name> <simulated_bedfile_name> <rum config> <name> [options]

<name> is some general name to identify this comparision.

The files <simulated_readsfile_name> and <simulated_bedfile_name> are
expected to be in <working_dir>, specify them without full path.

The bed file should have the sequence numbers in the first column.

This script takes the fasta and bed files output from the simulator,
runs RUM and generates a comparison.

Options: 
   -minidentity n : run blat with minidenitity equal to n, which must be an integer and 0<=n<=100

";
}

$minidentity = 93;
for($i=5; $i<@ARGV; $i++) {
    $optionrecognized = 0;
    if($ARGV[$i] eq "-minidentity") {
	$minidentity = $ARGV[$i+1];
	$i++;
	if(!($minidentity =~ /^\d+$/ && $minidentity <= 100)) {
	    die "\nERROR: minidentity must be an integer between zero and 100.\nYou have given '$minidentity'.\n\n";
	}
	$optionrecognized = 1;
    }
    if($optionrecognized == 0) {
	die "\nError: option '$ARGV[$i]' not recognized.\n";
    }
}

$working_dir = $ARGV[0];
$simulated_readsfile = $ARGV[1];
$rum_cofig = $ARGV[3];
`perl RUM_runner.pl $rum_cofig $working_dir/$simulated_readsfile $working_dir 1 $ARGV[4] -minidentity $minidentity`;
$simulated_bedfile = $ARGV[2];
$name = $ARGV[4];

open(OUTFILE, ">$working_dir/temp_truth_all.bed");
open(INFILE, "$working_dir/$simulated_bedfile");
while($line = <INFILE>) {
    $line =~ s/^\S+\t//;
    print OUTFILE $line;
}
close(OUTFILE);
close(INFILE);
`java -Xmx2000m M2C $working_dir/temp_truth_all.bed $working_dir/temp_truth_all.cov $working_dir/temp_truth_all.log -ucsc -chunks 2`;
 
open(INFILE, "$working_dir/RUM_Unique");
while($line = <INFILE>) {
    $line =~ /seq.(\d+)/;
    $uniquemappers_seqnums{$1}++;
}
close(INFILE);


open(INFILE, "$working_dir/$simulated_bedfile");
open(OUTFILE, ">$working_dir/temp_truth_unique-only.bed");
while($line = <INFILE>) {
    $line =~ /seq.(\d+)/;
    if(defined $uniquemappers_seqnums{$1}) {
	$line =~ s/^\S+\t//;
	print OUTFILE $line;
    }
}

`java -Xmx2000m M2C $working_dir/temp_truth_unique-only.bed $working_dir/temp_truth_unique-only.cov $working_dir/temp_truth_unique-only.log -ucsc -chunks 2`;

$RUM_cov = "RUM_" . $name . ".cov";

$filename = "$working_dir/RUM_comparison_" . $name . ".txt";
open(OUTFILE, ">$filename");

$comparator_output = `perl scripts/compare_covs.pl $working_dir/temp_truth_all.cov $working_dir/$RUM_cov -open`;
print OUTFILE "Comparison to all reads:\n";
print OUTFILE $comparator_output;
print OUTFILE "\n";

$comparator_output = `perl scripts/compare_covs.pl $working_dir/temp_truth_unique-only.cov $working_dir/$RUM_cov -open`;
print OUTFILE "\nComparison only on uniquely mapping reads:\n";
print OUTFILE $comparator_output;
print OUTFILE "\n";

close(OUTFILE);
