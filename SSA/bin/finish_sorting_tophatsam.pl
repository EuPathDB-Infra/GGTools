#!/usr/bin/perl

$nomap_record = "\t69\t*\t0\t255\t*\t*\t0\t0\tXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

open(INFILE, $ARGV[0]);
$flag = 0;
$line = <INFILE>;
chomp($line);
$line =~ /seq.(\d+)./;
$seqnum = $1;
$aline_cnt = 0;
$bline_cnt = 0;
while($flag == 0) {

    undef @a_lines;
    undef @b_lines;
    until(!($line =~ /seq.$seqnum/)) {
	$line =~ /seq.\d+(.)/;
	$type = $1;
	if($type eq 'a') {
	    $a_lines[$aline_cnt] = $line;
	    $aline_cnt++;
	}
	if($type eq 'b') {
	    $b_lines[$bline_cnt] = $line;
	    $bline_cnt++;
	}
	$line = <INFILE>;
	chomp($line);
	if($line eq '') {
	    $flag = 1;
	}
    }
    if($aline_cnt > $bline_cnt) {
	$max_cnt = $aline_cnt;
    } else {
	$max_cnt = $bline_cnt;
    }
    for($i=0; $i<$max_cnt; $i++) {
	if($i < @a_lines) {
	    print "$a_lines[$i]\n";
	} else {
	    print "seq.$seqnum";
	    print "a$nomap_record\n"
	}
	if($i < @b_lines) {
	    print "$b_lines[$i]\n";
	} else {
	    print "seq.$seqnum";
	    print "b$nomap_record\n"
	}
    }
    
    $line =~ /seq.(\d+)(.)/;
    $seqnum = $1;
    $type = $2;
    $aline_cnt = 0;
    $bline_cnt = 0;
}
