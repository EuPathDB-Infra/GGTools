#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennslyvania, 2010
$|=1;

if(@ARGV < 4) {
    die "
Usage: compare_2_vs_1.pl <quant_file1> <quant_file2> <quant_file3> <min_depth> [options]

Where:  <quant_file1> and <quant_file2> are feature quantification files for condition 1
        and <quant_file3> is a feature quantification file for condition 2

        <min_depth> is the minimum average depth-of-coverage that each gene must have
        in at least one of the two conditions in order for the gene to be considered.

This program returns features (genes, exons, introns) ranked by their likelihood
of being differenially expressed, as determined by the T-statistic.  Features returned
are ranked by T-stat. This is a pretty thin T-stat given only two reps in one condition
and one in the other, but that is the minimal necessary amount to get a valid definition
of the T-stat.  Note an adjustement factor is added to the denominator of the T-stat to
avoid blowing up near zero.

By default this program returns genes.  Use the options to output at the exon or intron level.

Options: -exons   : rank exons by fold change
         -introns : rank introns by fold change
         -annot f : f is a file of gene annotation with col 2 = ucsc id and col 3 = refseq id
                    col 1 is the annotation that is used in the output.

Note: -exons and -introns can both be set, but if either is set it will not report gene-level
results, to get genes and exons/introns you have to run twice.

";
}
$min_depth = $ARGV[3];
if(!($min_depth =~ /^\d+$/)) {
    die "\nError: <min_depth> must be a positive integer, not '$min_depth'\n\n";
}
$exons="false";
$introns="false";
for($i=4; $i<@ARGV; $i++) {
    $optionrecognized = 0;
    if($ARGV[$i] eq "-exons") {
	$exons = "true";
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-introns") {
	$introns = "true";
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-annot") {
 	$annotfile = $ARGV[$i+1];
	$i++;
	$optionrecognized = 1;
    }
    if($optionrecognized == 0) {
	print STDERR "\nError: option \'$ARGV[$i]\' is not recognized\n\n";
	exit(0);
    }
}

open(INFILE, $annotfile) or die "Error: Cannot open file '$annotfile' for reading.\n\n";
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    $ucsc{$a[1]} = $a[0];
    $refseq{$a[2]} = $a[0];
}
close(INFILE);

open(INFILE, $ARGV[0]);
$line = <INFILE>;
$line = <INFILE>;
$line = <INFILE>; # this is the ----------- ... line
$genecounter = 0;
$gene_ave = 0;
$exon_ave = 0;
$intron_ave = 0;
while($line = <INFILE>) { # this is the gene id line
    chomp($line);
    $geneid[$genecounter] = $line;
    $line = <INFILE>; # this is the header line that each gene has
    $line = <INFILE>; # this is the line for the full gene (sum of all exons)
    chomp($line);
    $line =~ s/ //g;
    @a = split(/\t/,$line);
    $gene_loc[$genecounter] = $a[1];
    $gene_ave_count[$genecounter] = $a[3];
    $gene_ave_norm[$genecounter] = $a[4];
    $gene_ave = $gene_ave + $a[4];
    $line = <INFILE>; # this is the first exon line
    chomp($line);
    until(($line =~ /----------/) || !($line =~ /\S/)) {
	$line =~ s/ //g;
	@a = split(/\t/,$line);
	if($a[0] =~ /exon/) {
	    $exon_ave_count{$a[1]} = $a[3];
	    $exon_ave_norm{$a[1]} = $a[4];
	}
	if($a[0] =~ /intron/) {
	    $intron_ave_count{$a[1]} = $a[3];
	    $intron_ave_norm{$a[1]} = $a[4];
	}
	$line = <INFILE>;
	chomp($line);
    }
    $genecounter++;
}
close(INFILE);
$gene_ave = $gene_ave / $genecounter;

open(INFILE, $ARGV[1]);
$line = <INFILE>;
$line = <INFILE>;
$line = <INFILE>; # this is the ----------- ... line
$genecounter = 0;
$gene_ave2 = 0;
while($line = <INFILE>) { # this is the gene id line
    chomp($line);
    $geneid2[$genecounter] = $line;
    $line = <INFILE>; # this is the header line that each gene has
    $line = <INFILE>; # this is the line for the full gene (sum of all exons)
    chomp($line);
    $line =~ s/ //g;
    @a = split(/\t/,$line);
    $gene_loc2[$genecounter] = $a[1];
    $gene_ave_count2[$genecounter] = $a[3];
    $gene_ave_norm2[$genecounter] = $a[4];
    $gene_ave2 = $gene_ave2 + $a[4];
    $line = <INFILE>; # this is the first exon line
    chomp($line);
    until(($line =~ /----------/) || !($line =~ /\S/)) {
	$line =~ s/ //g;
	@a = split(/\t/,$line);
	if($a[0] =~ /exon/) {
	    $exon_ave_count2{$a[1]} = $a[3];
	    $exon_ave_norm2{$a[1]} = $a[4];
	}
	if($a[0] =~ /intron/) {
	    $intron_ave_count2{$a[1]} = $a[3];
	    $intron_ave_norm2{$a[1]} = $a[4];
	}
	$line = <INFILE>;
	chomp($line);
    }
    $genecounter++;
}
close(INFILE);
$gene_ave2 = $gene_ave2 / $genecounter;

open(INFILE, $ARGV[2]);
$line = <INFILE>;
$line = <INFILE>;
$line = <INFILE>; # this is the ----------- ... line
$genecounter = 0;
$gene_ave3 = 0;
while($line = <INFILE>) { # this is the gene id line
    chomp($line);
    $geneid3[$genecounter] = $line;
    $line = <INFILE>; # this is the header line that each gene has
    $line = <INFILE>; # this is the line for the full gene (sum of all exons)
    chomp($line);
    $line =~ s/ //g;
    @a = split(/\t/,$line);
    $gene_loc3[$genecounter] = $a[1];
    $gene_ave_count3[$genecounter] = $a[3];
    $gene_ave_norm3[$genecounter] = $a[4];
    $gene_ave3 = $gene_ave3 + $a[4];
    $line = <INFILE>; # this is the first exon line
    chomp($line);
    until(($line =~ /----------/) || !($line =~ /\S/)) {
	$line =~ s/ //g;
	@a = split(/\t/,$line);
	if($a[0] =~ /exon/) {
	    $exon_ave_count2{$a[1]} = $a[3];
	    $exon_ave_norm3{$a[1]} = $a[4];
	}
	if($a[0] =~ /intron/) {
	    $intron_ave_count3{$a[1]} = $a[3];
	    $intron_ave_norm3{$a[1]} = $a[4];
	}
	$line = <INFILE>;
	chomp($line);
    }
    $genecounter++;
}
close(INFILE);
$gene_ave3 = $gene_ave3 / $genecounter;

$ratio_adjustment_factor = ($gene_ave + $gene_ave2 + $gene_ave3) / 30;

if($exons eq "false" && $introns eq "false") {
    print "GENE_ID\tLocation\tintensity1\tintensity2\tadj-ratio\n";
    for($i=0; $i<$genecounter; $i++) {
	$F = $gene_ave_norm[$i]-$gene_ave_norm2[$i];
	if($F < 0) {
	    $F = -1 * $F;
	}
	$gene_ratio[$i] = ( $gene_ave_norm3[$i] - ($gene_ave_norm[$i] + $gene_ave_norm2[$i])/2 ) / ($F + $ratio_adjustment_factor);
	if($gene_ratio[$i] < 0) {
	    $gene_ratio[$i] = -1 * $gene_ratio[$i];
	}
	if(($gene_ave_count[$i]+$gene_ave_count2[$i])/2 >= $min_depth || $gene_ave_count3[$i] >= $min_depth) {
	    $geneid[$i] =~ s/\s+(\+|-)\s*$//;
	    if($geneid[$i] =~ /([^:]*)\([^)]*ucsc[^)]*\)/) {
		$u = $ucsc{$1};
		if($u =~ /\S/) {
		    $name = $u;
		}
	    }
	    if($geneid[$i] =~ /([^:]*)\([^)]*refseq[^)]*\)/) {
		$u = $refseq{$1};
		if($u =~ /\S/) {
		    $name = $u;
		}
	    }
	    if($name =~ /\S/) {
		$gene_hash{"$name\t$gene_loc[$i]\t$gene_ave_norm[$i]\t$gene_ave_norm2[$i]\t$gene_ave_norm3[$i]"} = $gene_ratio[$i];
	    }
	    else {
		$gene_hash{"$geneid[$i]\t$gene_loc[$i]\t$gene_ave_norm[$i]\t$gene_ave_norm2[$i]\t$gene_ave_norm3[$i]"} = $gene_ratio[$i];
	    }
	}
    }
    foreach $key (sort {$gene_hash{$b}<=>$gene_hash{$a}} keys %gene_hash) {
	$gene_ratio = $gene_hash{$key};
	print "$key\t$gene_ratio\n";
    }
}
if($exons eq "true") {
    foreach $exon (keys %exon_ave_norm) {
	$F = $exon_ave_norm{$exon}-$exon_ave_norm2{$exon};
	if($F < 0) {
	    $F = -1 * $F;
	}
	$exon_ratio{$exon} = ( $exon_ave_norm3{$exon} - ($exon_ave_norm{$exon} + $exon_ave_norm2{$exon})/2 ) / ($F+$ratio_adjustment_factor);
	if($exon_ratio{$exon} < 0) {
	    $exon_ratio{$exon} = -1 * $exon_ratio{$exon};
	}
	if(($exon_ave_count{$exon}+$exon_ave_count2{$exon})/2 >= $min_depth || $exon_ave_count3{$exon} >= $min_depth) {
	    $feature_hash{"$exon\tEXON\t$exon_ave_norm{$exon}\t$exon_ave_norm2{$exon}\t$exon_ave_norm3{$exon}"} = $exon_ratio{$exon};
	}	
    }
}

if($introns eq "true") {
    foreach $intron (keys %intron_ave_norm) {
	$F = $intron_ave_norm{$intron}-$intron_ave_norm2{$intron};
	if($F < 0) {
	    $F = -1 * $F;
	}
	$intron_ratio{$intron} = ( $intron_ave_norm3{$intron} - ($intron_ave_norm{$intron} + $intron_ave_norm2{$intron})/2 ) / ($F+$ratio_adjustment_factor);
	if($intron_ratio{$intron} < 0) {
	    $intron_ratio{$intron} = -1 * $intron_ratio{$intron};
	}
	if(($intron_ave_count{$intron}+$intron_ave_count2{$intron})/2 >= $min_depth || $intron_ave_count3{$intron} >= $min_depth) {
	    $feature_hash{"$intron\tINTRON\t$intron_ave_norm{$intron}\t$intron_ave_norm2{$intron}\t$intron_ave_norm3{$intron}"} = $intron_ratio{$intron};
	}	
    }
}
if($exons eq "true" || $introns eq "true") {
    print "Location\ttype\tcond1.1\tcond1.2\tcond1.3\tstat\n";
    foreach $key (sort {$feature_hash{$b}<=>$feature_hash{$a}} keys %feature_hash) {
	$feature_ratio = $feature_hash{$key};
	print "$key\t$feature_ratio\n";
    }
}
