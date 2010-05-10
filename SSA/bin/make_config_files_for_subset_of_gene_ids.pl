#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennslyvania, 2010

if(@ARGV < 4) {
    die "
Usage:  make_config_files_for_subset_of_gene_ids.pl <stem> <ids> <configdir> <outdir>

  * <stem> is the suffix that will qualify these config files.

  * <ids> is a space separated list of ids or the name of a file of ids.

  * <configdir> is the directory where the master config files are.

  * <outdir> is the directory where the output files are to be written.

";
}

$stem = $ARGV[0];
if(-e $ARGV[1]) {
    open(INFILE, $ARGV[1]);
    $cnt=0;
    while($line = <INFILE>) {
	chomp($line);
	$ids{$line}++;
	$cnt++;
    }
}
else {
    for($cnt=1;$cnt<@ARGV;$cnt++) {
	$ids{$ARGV[$cnt]}++;
    }
}

$configdir = $ARGV[2];
if(-e $configdir) {
    die "Error: directory with config files '$configdir' does not seem to exist.\n\n";
}
$configdir =~ s!/$!!;
$configdir = $configdir . "/";
$outdir = $ARGV[3];
if(-e $outdir) {
    die "Error: output directory '$outdir' does not seem to exist.\n\n";
}
$outdir =~ s!/$!!;
$outdir = $outdir . "/";

$simulator_config_geneinfo_file = $configdir . "simulator_config_geneinfo";

open(INFILE, $simulator_config_geneinfo_file);
$filename = $outdir . "simulator_config_geneinfo_" . $stem;
open(OUTFILE, ">$filename");
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    @b = split(/::::/,$a[7]);
    $flag = 0;
    for($i=0; $i<@b; $i++) {
	$b[$i] =~ s/\(.*\)//;
	if($ids{$b[$i]} + 0 > 0) {
	    if($flag == 0) {
		print OUTFILE "$line\n";
		$flag = 1;
	    }
	}
    }
}
close(INFILE);
close(OUTFILE);

$simulator_config_featurequantifications_file = $configdir . "simulator_config_featurequantifications";

open(INFILE, $simulator_config_featurequantifications_file);
$filename = $outdir . "simulator_config_featurequantifications_" . $stem;
open(OUTFILE, ">$filename");
print OUTFILE "--------------------------------------------------------------------\n";
$line = <INFILE>;
until($line  =~ /---------------------/) {
    $line = <INFILE>;
}
while(1 == 1) {
    $line = <INFILE>;
    chomp($line);
    if($line eq '') {
	last;
    }
    $idline = $line;
    $idline =~ s/\s*(\+|-)\s*$//;
    @b = split(/::::/,$idline);
    $flag = 0;
    for($i=0; $i<@b; $i++) {
	$b[$i] =~ s/\(.*\)//;
	if($ids{$b[$i]} + 0 > 0) {
	    if($flag == 0) {
		print OUTFILE "$line\n";
		$flag = 1;
		until($line =~ /---------------------/ || $line eq '') {
		    $line = <INFILE>;
		    if($line eq '') {
			last;
		    }
		    print OUTFILE $line;
		    if($line =~ /intron/) {
			@a = split(/\t/,$line);
			$introns{$a[1]}++;
		    }
		}
	    }
	}
    }
}
close(INFILE);
close(OUTFILE);

$simulator_config_intronseq_file = $configdir . "simulator_config_intronseq";

open(INFILE, $simulator_config_intronseq_file);
$filename = $outdir . "simulator_config_intronseq_" . $stem;
open(OUTFILE, ">$filename");
$line = <INFILE>;
chomp($line);
$flag = 0;
while($flag == 0) {
    if($line =~ /^>/) {
	$id = $line;
	$id =~ s/^>//;
	if($introns{$id}+0>0) {
	    print OUTFILE "$line\n"; 
	    $line = <INFILE>;
	    chomp($line);
	    if($line eq '') {
		$flag = 1;
	    }
	    until($line =~ /^>/ || $flag == 1) {
		print OUTFILE "$line\n"; 
		$line = <INFILE>;
		chomp($line);
		if($line eq '') {
		    $flag = 1;
		}
	    }
	}
	else {
	    $line = <INFILE>;
	    chomp($line);
	    until($line =~ /^>/ || $flag == 1) {
		$line = <INFILE>;
		chomp($line);
		if($line eq '') {
		    $flag = 1;
		}
	    }
	}
    }
}
close(INFILE);
close(OUTFILE);

$simulator_config_geneseq_file = $configdir . "simulator_config_geneseq";

open(INFILE, $simulator_config_geneseq_file);
$filename = $outdir . "simulator_config_geneseq_" . $stem;
open(OUTFILE, ">$filename");
$line = <INFILE>;
chomp($line);
$flag = 0;
while($flag == 0) {
    if($line =~ /^>/) {
	$id = $line;
	$id =~ s/:[^:]+:[^:]+$//;
	$id =~ s/^>//;
	if($ids{$id}+0>0) {
	    print OUTFILE "$line\n"; 
	    $line = <INFILE>;
	    chomp($line);
	    if($line eq '') {
		$flag = 1;
	    }
	    until($line =~ /^>/ || $flag == 1) {
		print OUTFILE "$line\n"; 
		$line = <INFILE>;
		chomp($line);
		if($line eq '') {
		    $flag = 1;
		}
	    }
	}
	else {
	    $line = <INFILE>;
	    chomp($line);
	    until($line =~ /^>/ || $flag == 1) {
		$line = <INFILE>;
		chomp($line);
		if($line eq '') {
		    $flag = 1;
		}
	    }
	}
    }
}
close(INFILE);
close(OUTFILE);

