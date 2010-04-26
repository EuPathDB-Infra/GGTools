#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennslyvania, 2010

$|=1;
if(@ARGV<3) {
    die "
Usage: compare2simulated.pl <working_dir> <simulated_readsfile_name> <simulated_bedfile_name> <name>

<name> is some general name to identify this comparision.

The files <simulated_readsfile_name> and <simulated_bedfile_name> are
expected to be in <working_dir>, specify them without full path.

The bed file should have the sequence numbers in the first column.

This script takes the fasta and bed files output from the simulator,
runs RUM and generates a comparison.

";
}
$working_dir = $ARGV[0];
$simulated_readsfile = $ARGV[1];
`perl RUM_runner.pl rum.config_mm9 $working_dir/$simulated_readsfile $working_dir 1 $ARGV[3]`;
$simulated_bedfile = $ARGV[2];
$name = $ARGV[3];

open(OUTFILE, ">$working_dir/temp_truth_all.bed");
open(INFILE, "$working_dir/$simulated_bedfile");
while($line = <INFILE>) {
    $line =~ s/^\S+\t//;
    print OUTFILE $line;
}
close(OUTFILE);
close(INFILE);
`java -Xmx2000m M2C $working_dir/temp_truth_all.bed $working_dir/temp_truth_all.cov -ucsc -chunks 2`;
 
open(INFILE, "$working_dir/RUM_NU");
while($line = <INFILE>) {
    $line =~ /seq.(\d+)/;
    $NUseqnums{$1}++;
}
close(INFILE);

open(INFILE, "$working_dir/$simulated_bedfile");
open(OUTFILE, ">$working_dir/temp_truth_unique-only.bed");
while($line = <INFILE>) {
    $line =~ /seq.(\d+)/;
    if(!(defined $NUseqnums{$1})) {
	$line =~ s/^\S+\t//;
	print OUTFILE $line;
    }
}

`java -Xmx2000m M2C $working_dir/temp_truth_unique-only.bed $working_dir/temp_truth_unique-only.cov -ucsc -chunks 2`;

$RUM_cov = "RUM_" . $name . ".cov";

$comparator_output = `perl scripts/compare_covs.pl $working_dir/temp_truth_unique-only.cov $working_dir/$RUM_cov -open`;

$filename = "RUM_comparison_" . $name . ".txt";
open(OUTFILE, $filename);
print OUTFILE $comparator_output;
print OUTFILE "\n";
close(OUTFILE);
