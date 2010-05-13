#!/usr/bin/perl

# Written by Gregory R Grant
# University of Pennsylvania, 2010

if(@ARGV < 2) {
    die "
Usage: excise_regions_from_cov.pl <cov file> <regions file> [options]

Where:  <cov file> is a coverage file
        <regions file> is a file of regions to excise, each row should look like chr:start-end

Options: 
        -pad n:  excise also anything within n bases of the regions, default = 0

";
}

$pad = 0;
for($i=2; $i<@ARGV; $i++) {
    $optionrecognized = 0;
    if($ARGV[$i] eq "-pad") {
	$i++;
	$pad = $ARGV[$i];
	$optionrecognized = 1;
    }

    if($optionrecognized == 0) {
	die "\nERROR: option '$ARGV[$i-1] $ARGV[$i]' not recognized\n";
    }
}

open(COVFILE, $ARGV[0]);
open(REGIONS, $ARGV[1]);

$cnt = 0;
while($region = <REGIONS>) {
    chomp($regions);
    $region =~ /^(.*):(\d+)-(\d+)/;
    $chr[$cnt] = $1;
    $start[$cnt] = $2;
    $end[$cnt] = $3;
    $cnt++;
}
close(REGIONS);
while($line = <COVFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    $flag = 0;
    for($i=0; $i<$cnt; $i++) {
	if($a[1] <= ($end[$i]+$pad) && $a[2] >= ($start[$i]-$pad) && $a[0] eq $chr[$i]) {
	    $flag = 1;
	}
    }
    if($flag == 0) {
	print "$line\n";
    }
}
close(COVFILE);
