#/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

if(@ARGV<1 || $ARGV[0] eq "/help/") {
    print "\nUsage: featurequant2geneprofiles.pl <feature_quantification_files> [options]\n\n";
    print "Profiles are output for all genes/exons/introns, by default.  To change this use the options below.\n\n";
    print "<feature_quantification_files> is a space separated list of feature quantification files.\n\n";
    print "options:\n";
    print "     -genes    : output values for genes only\n";
    print "     -exons    : output values for exons only\n";
    print "     -introns  : output values for introns only\n";
    print "     -features : output values for exons and introns only\n";
    print "     -names=\"name1,,,name2,,, ... ,,,nameN\" : create a header line with these names.\n";
    print "     -simple   : if in -exon or -intron mode, then print out as two column table without qualifiers.\n";
    print "     -printheader   : print a header line.\n";
    print "\n";
    exit();
}

for($i=0; $i<@ARGV; $i++) {
    if(!($ARGV[$i] =~ /^-/)) {
	$numfiles = $i + 1;
    }
}
if($numfiles < 1) {
    print "ERROR: no valid files given.\n";
    exit();
}

$genesonly = "false";
$exonsonly = "false";
$intronsonly = "false";
$featuresonly = "false";
$all = "true";
$printheader = "false";
$simple = "false";
for($i=$numfiles; $i<@ARGV; $i++) {
    $optionrecognized = 0;
    if($ARGV[$i] eq "-genes") {
	$genesonly = "true";
	$all = "false";
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-exons") {
	$exonsonly = "true";
	$all = "false";
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-introns") {
	$intronsonly = "true";
	$all = "false";
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-features") {
	$featuresonly = "true";
	$all = "false";
	$optionrecognized = 1;
    }
    if($ARGV[$i] eq "-simple") {
	$simple = "true";
	$optionrecognized = 1;
    }
    if($ARGV[$i] =~ /names?=(.*)/) {
	@names = split(/,,,/,$1);
	if(@names != $numfiles) {
	    print "Error: number of names must equal the number of files.\n";
	    exit();
	}
	$printheader = "true";
	$optionrecognized = 1;
    }    
    if($ARGV[$i] eq "-printheader") {
	$printheader = "true";
	$optionrecognized = 1;
    }
    if($optionrecognized == 0) {
	print "\nERROR: option $ARGV[$i] not recognized\n";
	exit();
    }
}

for($i=0; $i<$numfiles; $i++) {
    $CNT = 0;
    if(!(-e $ARGV[$i])) {
	print "\nERROR: file $ARGV[$i] not found.\n\n";
	exit();
    }
    open(INFILE, $ARGV[$i]);
    while($line = <INFILE>) {
	chomp($line);
	if($line =~ /---------------------------/) {
	    $line = <INFILE>;
	    chomp($line);
	    $line =~ s/\t(\+|-)//;
	    $geneid = $line;
	    $ALL[$CNT][$i][0] = $geneid;
	}
	if($line =~ /gene/) {
	    if($genesonly eq "true" || $all eq "true") {
		@a = split(/\t/,$line);
		$profile{$geneid}[$i] = $a[4];
		$genelocation{$geneid}[$i] = $a[1];
	    }
	    $a[0] =~ s/^\s+//;
	    $a[0] =~ s/\s+$//;
	    $ALL[$CNT][$i][1] = $a[1];
	    $ALL[$CNT][$i][2] = $a[4];
	    $CNT++;
	}
	if($line =~ /exon/) {
	    if($exonsonly eq "true" || $all eq "true" || $featuresonly eq "true") {
		@a = split(/\t/,$line);
		$exon{$a[1]}[$i] = $a[4];
	    }
	    $a[0] =~ s/^\s+//;
	    $a[0] =~ s/\s+$//;
	    $ALL[$CNT][$i][0] = $a[0];
	    $ALL[$CNT][$i][1] = $a[1];
	    $ALL[$CNT][$i][2] = $a[4];
	    $CNT++;
	}
	if($line =~ /intron/) {
	    if($intronsonly eq "true" || $all eq "true" || $featuresonly eq "true") {
		@a = split(/\t/,$line);
		$intron{$a[1]}[$i] = $a[4];
	    }
	    $a[0] =~ s/^\s+//;
	    $a[0] =~ s/\s+$//;
	    $ALL[$CNT][$i][0] = $a[0];
	    $ALL[$CNT][$i][1] = $a[1];
	    $ALL[$CNT][$i][2] = $a[4];
	    $CNT++;
	}
    }
}

if($printheader eq "true") {
    if($simple eq "false") {
	print "\t";
    }
    print "name";
    for($i=0; $i<$numfiles; $i++) {
	if($names[$i] =~ /\S/) {
	    print "\t$names[$i]";
	} else {
	    $j = $i+1;
	    print "\tfile_$j";
	}
    }
    print "\n";
}

if($all eq "true") {
    for($cnt=0; $cnt<$CNT; $cnt++) {
	print "$ALL[$cnt][0][0]\t$ALL[$cnt][0][1]";
	for($i=0; $i<$numfiles;$i++) {
	    print "\t$ALL[$cnt][$i][2]";
	}
	print "\n";
    }
}

if($genesonly eq "true") {
    foreach $geneid (keys %genelocation) {
	print "$geneid";
	if($locations eq "true") {
	    print "\t$genelocation{$geneid}";
	}
	for($i=0; $i<$numfiles; $i++) {
	    print "\t$profile{$geneid}[$i]";
	}
	print "\n";
    }
}

if($exonsonly eq "true" || $featuresonly eq "true") {
    foreach $exonid (sort {$exon{$a}<=>$exon{$b}} keys %exon) {
	if($simple eq "false") {
	    print "$exonid\tEXON";
	} else {
	    print "$exonid";
	}
	for($i=0; $i<$numfiles; $i++) {
	    print "\t$exon{$exonid}[$i]";
	}
	print "\n";
    }
}

if($intronsonly eq "true" || $featuresonly eq "true") {
    foreach $intronid (sort {$intron{$a}<=>$intron{$b}} keys %intron) {
	if($simple eq "false") {
	    print "$intronid\tINTRON";
	} else {
	    print "$intronid";
	}
	for($i=0; $i<$numfiles; $i++) {
	    print "\t$intron{$intronid}[$i]";
	}
	print "\n";
    }
}


# --------------------------------------------------------------------
# PFD0028w        +
#     Type        Location                Count   Ave_Cnt Ave_Nrm Length
#     gene        Pf3D7_04:67853-67936    0       0       0       84
#   exon 1        Pf3D7_04:67853-67936    0       0       0       84
# --------------------------------------------------------------------
