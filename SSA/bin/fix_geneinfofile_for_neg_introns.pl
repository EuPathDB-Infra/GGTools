#!/usr/bin/perl

# Written by Gregory R. Grant
# Universiry of Pennsylvania, 2010

if(@ARGV < 1) {
    die "
Usage: fix_geneinfofile_for_neg_introns.pl <gene info file>

This script takes a UCSC gene annotation file and outputs a file that removes
introns of zero or negative length.  You'd think there shouldn't be such introns
but for some annotation sets there are.

The annotation file has to be downloaded with the following fields:
1) name
2) chrom
3) strand
4) txStart
5) txEnd
6) cdsStart
7) cdsEnd
8) exonCount
9) exonStarts
10) exonEnds

This script is part of the pipeline of scripts used to create RUM indexes.
For more information see the library file: 'how2setup_genome-indexes_forPipeline.txt'.

";
}

open(INFILE, $ARGV[0]);
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/, $line);
    $chr = $a[0];
    $starts = $a[5];
    $ends = $a[6];
    $starts =~ s/,\s*$//;
    $ends =~ s/,\s*$//;
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
#       print "length = $length\n";
        if($length > 0) {
            $start_string = $start_string . $S[$i] . ",";
            $end_string = $end_string . $E[$i-1] . ",";
        }
        else {
#           print "$line\n";
            $a[4]--;
        }
    }
    $end_string = $end_string . $E[$N-1] . ",";;
    print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$start_string\t$end_string\t$a[7]\n";
}


# chr1    -       4481008 4486494 5       4481008,4483180,4483852,4485216,4486371,        4482749,4483547,4483944,4486023,4486494,        uc007aez.1(mm9_ucsc_known_genes.txt)::::NM_011441(mm9_refseq_genes.txt)

# chr2    +       67011890        67024269        3       67011890,67022499,67023182,     67012086,67023182,67024269,     uc008jxi.1(mm9_ucsc_known_genes.txt)
