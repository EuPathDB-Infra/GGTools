open(INFILE, "UCSC_mouse_knowngene_id_mapping");
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    $ucsc{$a[1]} = $a[0];
    $refseq{$a[2]} = $a[0];
}
close(INFILE);

open(INFILE, $ARGV[0]);
$line = <INFILE>;
$line = <INFILE>;
$line = <INFILE>; # this is the ----------- ... line
$genecounter = 0;
while($line = <INFILE>) { # this is the gene id line
    chomp($line);
    $geneid[$genecounter] = $line;
    $line = <INFILE>; # this is the header line that each gene has
    $line = <INFILE>; # this is the line for the full gene (sum of all exons)
    chomp($line);
    $line =~ s/ //g;
    @a = split(/\t/,$line);
    $gene_loc[$genecounter] = $a[1];
    $gene_ave_count[$genecounter] = $a[3];
    $gene_ave_nrm[$genecounter] = $a[4];
    $line = <INFILE>; # this is the first exon line
    chomp($line);
    until(($line =~ /----------/) || !($line =~ /\S/)) {
	$line =~ s/ //g;
	@a = split(/\t/,$line);
	if($a[0] =~ /exon/) {
	    $a[0] =~ /(\d+)$/;
	    $exonnumber = $1 - 1;
	    $exon_loc[$genecounter][$exonnumber] = $a[1];
	    $exon_ave_count[$genecounter][$exonnumber] = $a[3];
	    $exon_ave_norm[$genecounter][$exonnumber] = $a[4];
	    $num_exons[$genecounter]++;
	}
	if($a[0] =~ /intron/) {
	    $a[0] =~ /(\d+)$/;
	    $intronnumber = $1 - 1;
	    $intron_loc[$genecounter][$intronnumber] = $a[1];
	    $intron_ave_count[$genecounter][$intronnumber] = $a[3];
	    $intron_ave_norm[$genecounter][$intronnumber] = $a[4];
	    $num_introns[$genecounter]++;
	}
	$line = <INFILE>;
	chomp($line);
    }
    $genecounter++;
}
close(INFILE);

open(INFILE, $ARGV[1]);
$line = <INFILE>;
$line = <INFILE>;
$line = <INFILE>; # this is the ----------- ... line
$genecounter = 0;
while($line = <INFILE>) { # this is the gene id line
    chomp($line);
    $geneid2[$genecounter] = $line;
    $line = <INFILE>; # this is the header line that each gene has
    $line = <INFILE>; # this is the line for the full gene (sum of all exons)
    chomp($line);
    $line =~ s/ //g;
    @a = split(/\t/,$line);
    $gene_loc2[$genecounter] = $a[1];
    $gene_ave_count2[$genecounter] = $a[3];
    $gene_ave_nrm2[$genecounter] = $a[4];
    $line = <INFILE>; # this is the first exon line
    chomp($line);
    until(($line =~ /----------/) || !($line =~ /\S/)) {
	$line =~ s/ //g;
	@a = split(/\t/,$line);
	if($a[0] =~ /exon/) {
	    $a[0] =~ /(\d+)$/;
	    $exonnumber = $1 - 1;
	    $exon_loc2[$genecounter][$exonnumber] = $a[1];
	    $exon_ave_count2[$genecounter][$exonnumber] = $a[3];
	    $exon_ave_norm2[$genecounter][$exonnumber] = $a[4];
	    $num_exons2[$genecounter]++;
	}
	if($a[0] =~ /intron/) {
	    $a[0] =~ /(\d+)$/;
	    $intronnumber = $1 - 1;
	    $intron_loc2[$genecounter][$intronnumber] = $a[1];
	    $intron_ave_count2[$genecounter][$intronnumber] = $a[3];
	    $intron_ave_norm2[$genecounter][$intronnumber] = $a[4];
	    $num_introns2[$genecounter]++;
	}
	$line = <INFILE>;
	chomp($line);
    }
    $genecounter++;
}
close(INFILE);

print "GENE_ID\tLocation\tWT_intensity\tMUT_intensity\tratio\n";
for($i=0; $i<$genecounter; $i++) {
    if($gene_ave_nrm2[$i] > 0) {
	$gene_ratio[$i] = $gene_ave_nrm[$i] / $gene_ave_nrm2[$i];
    }
    else {
	$gene_ratio[$i] = 10000;
    }
    if($gene_ratio[$i] < 1) {
	if($gene_ratio[$i] == 0) {
	    $gene_ratio[$i] = 10000;
	}
	else {
	    $gene_ratio[$i] = 1/$gene_ratio[$i];
	}
    }
    if($gene_ave_count[$i] >= 5 || $gene_ave_count2[$i] >= 5) {
	$geneid[$i] =~ s/\s+(\+|-)\s*$//;
	if($geneid[$i] =~ /([^:]*)\(ucscgenes\)/) {
	    $u = $ucsc{$1};
	    if($u =~ /\S/) {
		$name = $u;
	    }
	}
	if($geneid[$i] =~ /([^:]*)\(refseq\)/) {
	    $u = $refseq{$1};
	    if($u =~ /\S/) {
		$name = $u;
	    }
	}
	if($name =~ /\S/) {
	    $gene_hash{"$name\t$gene_loc[$i]\t$gene_ave_nrm[$i]\t$gene_ave_nrm2[$i]"} = $gene_ratio[$i];
	}
	else {
	    $gene_hash{"$geneid[$i]\t$gene_loc[$i]\t$gene_ave_nrm[$i]\t$gene_ave_nrm2[$i]"} = $gene_ratio[$i];
	}
    }
#    print "$geneid[$i]\t$gene_loc[$i]\t$gene_ratio[$i]\n";
#    print "geneid[$i] = $geneid[$i]\t$geneid2[$i]\n";
#    print "gene_loc[$i] = $gene_loc[$i]\t$gene_loc2[$i]\n";
#    print "gene_ave_count[$i] = $gene_ave_count[$i]\t$gene_ave_count2[$i]\n";
#    print "gene_ave_nrm[$i] = $gene_ave_nrm[$i]\t$gene_ave_nrm2[$i]\n";
    for($j=0; $j<$num_exons[$i]; $j++) {
	if($exon_ave_norm2[$i][$j] > 0) {
	    $exon_ratio[$i][$j] = $exon_ave_norm[$i][$j] / $exon_ave_norm2[$i][$j];
	}
	else {
	    $exon_ratio[$i][$j] = 10000;
	}
	$k = $j+1;
#	print "exon $k\t$exon_loc[$i][$j]\t$exon_ratio[$i][$j]\n";
#	print "exon_loc[$i][$j] = $exon_loc[$i][$j]\t$exon_loc2[$i][$j]\n";
#	print "exon_ave_count[$i][$j] = $exon_ave_count[$i][$j]\t$exon_ave_count2[$i][$j]\n";
#	print "exon_ave_norm[$i][$j] = $exon_ave_norm[$i][$j]\t$exon_ave_norm2[$i][$j]\n";
    }
    for($j=0; $j<$num_introns[$i]; $j++) {
	if($intron_ave_norm2[$i][$j] > 0) {
	    $intron_ratio[$i][$j] = $intron_ave_norm[$i][$j] / $intron_ave_norm2[$i][$j];
	}
	else {
	    $gene_ratio[$i] = 10000;
	}
	$k = $j+1;
#	print "intron $k\t$intron_loc[$i][$j]\t$intron_ratio[$i][$j]\n";
#	print "intron_loc[$i][$j] = $intron_loc[$i][$j]\t$intron_loc2[$i][$j]\n";
#	print "intron_ave_count[$i][$j] = $intron_ave_count[$i][$j]\t$intron_ave_count2[$i][$j]\n";
#	print "intron_ave_norm[$i][$j] = $intron_ave_norm[$i][$j]\t$intron_ave_norm2[$i][$j]\n";
    }
}
#print "-------------------------\n";

foreach $key (sort {$gene_hash{$b}<=>$gene_hash{$a}} keys %gene_hash) {
    $gene_ratio = $gene_hash{$key};
    print "$key\t$gene_ratio\n";
}
