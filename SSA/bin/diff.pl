
if(@ARGV<2) {
    die "
Usage: perl diff.pl <file1> <file2> [options]

Reports the first line on which file1 and file2 differ, if any.

Options: -all  : show all lines that are different.

";
}

open(INFILE1, $ARGV[0]);
open(INFILE2, $ARGV[1]);

$all = "false";
for($i=2; $i<@ARGV; $i++) {
    if($ARGV[$i] eq '-all') {
	$all = "true";
    }
}

$cnt = 1;
while($line1 = <INFILE1>) {
    if($cnt % 100000 == 0) {
	print STDERR "$cnt\n";
    }
    $line2 = <INFILE2>;
    if($line2 eq '') {
	print "Reached the end of '$ARGV[2]'\n";
	last;
    }
    if(!($line1 eq $line2)) {
	print "-----------------\nline $cnt:\n\n";
	print "$ARGV[0]:\n$line1\n";
	print "$ARGV[1]:\n$line2\n";
	if($all eq 'false') {
	    print "-----------------\n";
	    exit();
	}
    }
    $cnt++;
}
if($line1 eq '') {
    $line2 = <INFILE2>;
    if($line2 =~ /\S/) {
	print "Reached the end of '$ARGV[1]'\n";
    } else {
	print "Reached the end of both files\n";
    }
}
close(INFILE1);
close(INFILE2);
print "-----------------\n";
