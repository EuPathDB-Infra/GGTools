
if(@ARGV<2) {
    die "
Usage: fix_sam.pl <sam file> <last seq num>

";
}

$samfile = $ARGV[0];
$samfile_idsfixed = $samfile . "_idsfixed";

open(INFILE, $samfile);
open(OUTFILE, ">$samfile_idsfixed");
while($line = <INFILE>) {
    if(!($line =~ /^@/)) {
	$line =~ s/\[.*\]//; # this fixes splicemap
	$line =~ s/^\d+~//; # this and the next two fix mapsplice
	$line =~ s!/1!a!;
	$line =~ s!/2!b!;
	print OUTFILE $line;
    }
}
close(INFILE);
close(OUTFILE);

$samfile_sorted = $samfile . "_sorted";
`perl sort_where_lines_start_seq.numa_or_seq.numb.pl $samfile_idsfixed $samfile_sorted`;

open(INFILE, $samfile_sorted);
$lastseqnum = $ARGV[1];

$flag = 0;
$line = <INFILE>;
chomp($line);
$acnt = 0;
$bcnt = 0;
$seqnum = 0;
while($flag == 0) {
    $line =~ /^seq.(\d+)(.)/;
    $seqnum_prev = $seqnum;
    $seqnum = $1;
    if($seqnum > $seqnum_prev + 1) {
	for($i=$seqnum_prev+1; $i<$seqnum; $i++) {
	    print "seq.$i";
	    print "a\t141\t*\t0\t255\t*\t*\t0\t0\tNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
	    print "seq.$i";
	    print "b\t141\t*\t0\t255\t*\t*\t0\t0\tNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
	}
    }

    $seqnuma = $seqnum . "a";
    $seqnumb = $seqnum . "b";
    $type = $2;
    if($type eq "a") {
	$a[$acnt] = $line;
	$acnt++;
    }
    if($type eq "b") {
	$b[$bcnt] = $line;
	$bcnt++;
    }
    $line = <INFILE>;
    chomp($line);
    while($line =~ /^seq.$seqnuma/ || $line =~ /^seq.$seqnumb/) {
	$line =~ /^seq.\d+(.)/;
	$type = $1;
	if($type eq "a") {
	    $a[$acnt] = $line;
	    $acnt++;
	}
	if($type eq "b") {
	    $b[$bcnt] = $line;
	    $bcnt++;
	}
	$line = <INFILE>;
	chomp($line);
    }
    if($acnt <= $bcnt) {
	for($i=0; $i<$acnt; $i++) {
	    $j=$i+1;
	    print "$a[$i]\tIH:i:$bcnt\tHI:i:$j\n$b[$i]\tIH:i:$bcnt\tHI:i:$j\n";
	}
	for($i=$acnt; $i<$bcnt; $i++) {
	    $j=$i+1;
	    print "seq.$seqnum";
	    print "a\t141\t*\t0\t255\t*\t*\t0\t0\tNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tIH:i:$bcnt\tHI:i:$j\n";
	    print "$b[$i]\tIH:i:$bcnt\tHI:i:$j\n";
	}
    }
    if($acnt > $bcnt) {
	for($i=0; $i<$bcnt; $i++) {
	    $j=$i+1;
	    print "$a[$i]\tIH:i:$acnt\tHI:i:$j\n$b[$i]\tIH:i:$acnt\tHI:i:$j\n";
	}
	for($i=$bcnt; $i<$acnt; $i++) {
	    $j=$i+1;
	    print "$a[$i]\tIH:i:$acnt\tHI:i:$j\n";
	    print "seq.$seqnum";
	    print "b\t141\t*\t0\t255\t*\t*\t0\t0\tNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tIH:i:$acnt\tHI:i:$j\n";
	}
    }
    undef @a;
    undef @b;
    $acnt = 0;
    $bcnt = 0;
    if($line eq '') {
	$flag = 1;
    }
}
if($seqnum < $lastseqnum) {
    for($i=$seqnum+1; $i<=$lastseqnum; $i++) {
	print "seq.$i";
	print "a\t141\t*\t0\t255\t*\t*\t0\t0\tNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
	print "seq.$i";
	print "b\t141\t*\t0\t255\t*\t*\t0\t0\tNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    }
}
