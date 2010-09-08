# chrI    +       334     649     1       334,    649,    YAL069W(sgd)::::YAL069W(ensemble)
use Roman;

open(INFILE, $ARGV[0]);
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    $chr = $a[0];
    $start = $a[2];
    $end = $a[3];
    $hash{$chr}{$line}[0] = $start;
    $hash{$chr}{$line}[1] = $end;
}
close(INFILE);

foreach $chr (sort cmpChrs keys %hash) {
    foreach $line (sort {$hash{$chr}{$a}[0]<=>$hash{$chr}{$b}[0] || ($hash{$chr}{$a}[0]==$hash{$chr}{$b}[0] && $hash{$chr}{$a}[1]<=>$hash{$chr}{$b}[1])} keys %{$hash{$chr}}) {
	chomp($line);
	if($line =~ /\S/) {
	    print $line;
	    print "\n";
	}
    }
}

sub cmpChrs () {
    $a2_c = lc($b);
    $b2_c = lc($a);
    if($a2_c =~ /chr(\d+)$/ && !($b2_c =~ /chr(\d+)$/)) {
        return 1;
    }
    if($b2_c =~ /chr(\d+)$/ && !($a2_c =~ /chr(\d+)$/)) {
        return -1;
    }
    if($a2_c =~ /chr([a-z])$/ && !($b2_c =~ /chr(\d+)$/) && !($b2_c =~ /chr[a-z]+$/)) {
        return 1;
    }
    if($b2_c =~ /chr([a-z])$/ && !($a2_c =~ /chr(\d+)$/) && !($a2_c =~ /chr[a-z]+$/)) {
        return -1;
    }
    if($a2_c =~ /chr([ivx]+)/ && $b2_c =~ /chr([ivx]+)/) {
	$a2_c =~ /chr([ivx]+)/;
	$a2_roman = $1;
	$b2_c =~ /chr([ivx]+)/;
	$b2_roman = $1;
	$a2_arabic = arabic($a2_roman);
    	$b2_arabic = arabic($b2_roman);
	if($a2_arabic >= $b2_arabic) {
	    return -1;
	} else {
	    return 1;
	}
    }
    if($b2_c =~ /chr([ivx]+)/ && $a2_c =~ /chr([ivx]+)/) {
	$b2_c =~ /chr([ivx]+)/;
	$b2_roman = $1;
	$a2_c =~ /chr([ivx]+)/;
	$a2_roman = $1;
	$b2_arabic = arabic($b2_roman);
    	$a2_arabic = arabic($a2_roman);
	if($b2_arabic >= $a2_arabic) {
	    return 1;
	} else {
	    return -1;
	}
    }
    if($b2_c =~ /chr([ivx]+)/ && !($a2_c =~ /chr([a-z]+)/) && !($a2_c =~ /chr(\d+)/)) {
	return -1;
    }
    if($a2_c =~ /chr([ivx]+)/ && !($b2_c =~ /chr([a-z]+)/) && !($b2_c =~ /chr(\d+)/)) {
	return 1;
    }
    if($b2_c =~ /chr([ivx]+)/ && ($a2_c =~ /chrm/)) {
	return -1;
    }
    if($a2_c =~ /chr([ivx]+)/ && ($b2_c =~ /chrm/)) {
	return 1;
    }
    if($a2_c =~ /chr[^xy\d]/ && (($b2_c =~ /chrx/) || ($b2_c =~ /chry/))) {
        return -1;
    }
    if($b2_c =~ /chr[^xy\d]/ && (($a2_c =~ /chrx/) || ($a2_c =~ /chry/))) {
        return 1;
    }

    if($a2_c =~ /chr(\d+)/) {
        $numa = $1;
        if($b2_c =~ /chr(\d+)/) {
            $numb = $1;
            if($numa <= $numb) {return 1;} else {return -1;}
        } else {
            return 1;
        }
    }
    if($a2_c =~ /chr([a-z]+)/) {
        $letter_a = $1;
        if($b2_c =~ /chr([a-z]+)/) {
            $letter_b = $1;
            if($letter_a le $letter_b) {return 1;} else {return -1;}
        } else {
            return -1;
        }
    }
    $flag_c = 0;
    while($flag_c == 0) {
        $flag_c = 1;
        if($a2_c =~ /^([^\d]*)(\d+)/) {
            $stem1_c = $1;
            $num1_c = $2;
            if($b2_c =~ /^([^\d]*)(\d+)/) {
                $stem2_c = $1;
                $num2_c = $2;
                if($stem1_c eq $stem2_c && $num1_c < $num2_c) {
                    return 1;
                }
                if($stem1_c eq $stem2_c && $num1_c > $num2_c) {
                    return -1;
                }
                if($stem1_c eq $stem2_c && $num1_c == $num2_c) {
                    $a2_c =~ s/^$stem1_c$num1_c//;
                    $b2_c =~ s/^$stem2_c$num2_c//;
                    $flag_c = 0;
                }
            }
        }
    }

    return 1;
}
