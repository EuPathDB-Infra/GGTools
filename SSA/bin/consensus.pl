#!/usr/bin/perl

$| = 1;

if(@ARGV < 1) {
    die "
Usage: consensus.pl <RUM file> <chr> <start> <end> <genome>

Where: <RUM file> is a RUM_Unique or RUM_NU file
       <chr> <start> <end> specify the span you want the consensus of
       <genome> is a fasta file of the genome, with sequence on one line

";
}

open(INFILE, $ARGV[0]);

$CHR = $ARGV[1];
$START = $ARGV[2];
$END = $ARGV[3];
open(GENOMESEQ, $ARGV[4]);

print STDERR "Going to find $CHR in the genome, this could take a few minutes...\n\n";
$FLAG = 0;
while($FLAG == 0) {
    $line = <GENOMESEQ>;
    if($line eq '') {
	$FLAG = 1;
    } else {
	chomp($line);
	$line =~ />(.*):1-(\d+)_strand=./;
	$chr = $1;
	$ref_seq = <GENOMESEQ>;
	chomp($ref_seq);
	print STDERR "chr $chr\n";
	if($chr eq $CHR) {
	    $CHRSEQ = $ref_seq;
	    $FLAG = 1;
	    last;
	}
    }
}
close(GENOMESEQ);

$SEQ = substr($CHRSEQ, $START-1, $END-$START+1);
@REFSEQ = split(//,$SEQ);
for($i=0; $i<@REFSEQ; $i++) {
    $REF{$i+$START} = $REFSEQ[$i];
}

while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    $chr = $a[1];
    if($chr eq $CHR) {
	@spans = split(/, /, $a[2]);
	for($i=0; $i<@spans; $i++) {
	    @b = split(/-/,$spans[$i]);
	    if($b[1] >= $START && $b[0] <= $END) {
#		print "$line\n";
		$a[3] =~ s/\+.*\+//;
		@seqs = split(/:/, $a[3]);
		$seq = $seqs[$i];
		@S = split(//,$seq);
		for($j=$b[0]; $j<$b[1]; $j++) {
		    $data{$j}{$S[$j-$b[0]]}++;
		}
	    }
	}
    }
}
close(INFILE);

for($i=$START; $i<=$END; $i++) {
    if(defined $data{$i}) {
	print "----\nlocation: $i ($REF{$i})\n";
	$max = 0;
	$min = 1000000000;
	foreach $key (keys %{$data{$i}}) {
	    if($data{$i}{$key} > $max) {
		$max = $data{$i}{$key};
		$maxbase = $key;
	    }
	    if($data{$i}{$key} < $min) {
		$min = $data{$i}{$key};
		$minbase = $key;
	    }
	    print "$key: $data{$i}{$key}\n";
	}
	if($maxbase ne $REF{$i}) {
	    print "SNP CANDIDATE\n";
	}
	if($maxbase ne $minbase && $max/$min <= 5 && $min >= 2) {
	    print "ALLELE SPECIFIC CANDIDATE\n";
	}
    }
}
