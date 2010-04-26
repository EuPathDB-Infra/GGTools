#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

# assuming a line of $ARGV[0] looks like this:
# seq.13392702    chr3    40576630-40576665, 40583257-40583328

if(!($ARGV[0] =~ /\S/)) {
    print STDERR "\nUsage: perl make_bed.pl infile [outfile]\n\nThis makes a bed file with one span per line (strand col = + for all rows), from a file\nthat has rows of tab delimited entries that look like this:\n\nseq.13551265    chr13   57531714-57531787, 57537471-57537504\n\nor this\n\nchr13   57531714-57531787, 57537471-57537504\n\nthe first argument is input file, the second is the output file, if no\nsecond argument then will output to std out.\n";
}
open(INFILE, $ARGV[0]);
$outmode = "screen";
if($ARGV[1] =~ /\S/) {
    open(OUTFILE, ">$ARGV[1]");
    $outmode = "file";
}

while($line = <INFILE>) {
    chomp($line);
    $line =~ s/seq.\d+.?\t//;
    @a = split(/\t/,$line);
    @b = split(/, /,$a[1]);
    $N = @b;
    for($i=0; $i<$N; $i++) {
        @c = split(/-/,$b[$i]);
        if($outmode eq "file") {
            print OUTFILE "$a[0]\t$c[0]\t$c[1]\t+\n";
        }
        else {
            print "$a[0]\t$c[0]\t$c[1]\t+\n";
        }
    }
}
