if(@ARGV<1) {
    die "
Usage: make_sam_header.pl <chrom info>

Where <chrom info> is the tab delimited file with two columns: chr, size

This file can be obtained from UCSC by going to the table browser
and choosing group: \"All Tables\", table: \"chrominfo\", output format:
\"selected fields from primary and related tables\" and on the next
page choose 'chrom' and 'size'.

";
}

open(INFILE, $ARGV[0]);
while($line = <INFILE>) {
    chomp($line);
    if(!($line =~ /^#/)) {
	$line =~ /(.*)\s+(\d+)/;
	$chr = $1;
	$size = $2;
	print "\@SQ\tSN:$chr\tLN:$size\n";
    }
}
close(INFILE);
