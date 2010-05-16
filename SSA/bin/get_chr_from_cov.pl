#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

$|=1;

if(@ARGV < 2) {
    die "
Usage: get_chr_from_cov.pl <cov file> <chr> [options]

Output is to standard out, redirect to a file using \">\"

Options: -name x  : change the name of the track to x. 

If you want to create a file for more than one chromosome, first do each
separately and give each a different name using -name, and then concatenate
the files together using:
> cat file1 file2 .... filen > merged_file

";
}

for($i=0; $i<@ARGV; $i++) {
    $optionrecognized=0;
    if($ARGV[$i] eq "-name") {
	$name = $ARGV[$i+1];
	$i++;
	$optionrecognized=1;
    }
}

open(INFILE, $ARGV[0]);
# line of input (wig file) looks like this:
# chr9    3342879 3342987 1
$chr = $ARGV[1];
$line = <INFILE>;
if($name =~ /\S/) {
    $line =~ s/(name=\"[^\"]*\")/name=\"$name\"/;
    $line =~ s/(description=\"[^\"]*\")/description=\"$name\"/;
}
print $line;
while($line = <INFILE>) {
    if($line =~ /^$chr\t/) {
	print $line;
    }
}
