#/usr/bin/perl

if(@ARGV != 1) {
    die "
Usage: sort_bowtie.pl <file>

Where <file> is the name of the file that has been output by Bowtie and the
read names are expected to be of the form seq.Na and seq.Nb where N is a
positive integer.

This script sorts the rows so that rows with seq.Na and seq.Nb are before
rows with seq.Ma and seqMb whenever N<M.  Furthermore for fixed N seq.Na is
before seq.Nb.

";

}

$filename = $ARGV[0];
open(INFILE, $filename);
$tempfilename = $filename . "_temp1." . $chunk;
open(OUTFILE, ">$tempfilename");
while($line = <INFILE>) {
    chomp($line);
    $line =~ s/^([^\t]+)\t//;
    $name = $1;
    $name =~ s/seq.//;
    $name =~ /(\d+)(a|b)/;
    print OUTFILE "$1\t$2\t$line\n";
}
close(OUTFILE);
close(INFILE);
$tempfilename2 = $filename . "_temp2." . $chunk;
$x = `sort -T . -n $tempfilename > $tempfilename2`;
$x = `rm $tempfilename`;
open(INFILE, $tempfilename2);
$sortedfilename = $type . "_sorted." . $chunk;
open(OUTFILE, ">$sortedfilename");
while($line = <INFILE>) {
    chomp($line);
    $line =~ s/^(\d+)\t(.)//;
    print OUTFILE "seq.$1$2$line\n";
}
close(OUTFILE);
close(INFILE);
$x = `rm $tempfilename2`;
