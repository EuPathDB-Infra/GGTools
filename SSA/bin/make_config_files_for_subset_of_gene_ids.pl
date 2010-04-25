#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennslyvania, 2010

if(@ARGV < 1) {
    die "
Usage:  make_config_files_for_subset_of_gene_ids.pl <stem> <ids>

  * <stem> is the suffix that will qualify these config files.

  * <ids> is a space separated list of ids or the name of a file of ids.

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

open(INFILE, "simulator_config_geneinfo");
$filename = "simulator_config_geneinfo_" . $stem;
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

open(INFILE, "simulator_config_featurequantifications");
$filename = "simulator_config_featurequantifications_" . $stem;
open(OUTFILE, ">$filename");
print OUTFILE "--------------------------------------------------------------------\n";
while($line = <INFILE>) {
    if($line =~ /---------------------/) {
	$line = <INFILE>;
	chomp($line);
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
}
close(INFILE);
close(OUTFILE);

open(INFILE, "simulator_config_intronseq");
$filename = "simulator_config_intronseq_" . $stem;
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

open(INFILE, "simulator_config_geneseq");
$filename = "simulator_config_geneseq_" . $stem;
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

