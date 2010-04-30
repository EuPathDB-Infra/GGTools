$| = 1;

if(@ARGV < 1) {
    print "\nusage: reads_simulator.pl <num_reads> <name> [options]\n\n";
    print "<name> str is used in the names of the output files.\nUse only alphanumeric, underscores and dashes.\n";

    print "\nThis program outputs a fasta file of reads and a bed file representing the truth of where these reads map.\n";
    print "   - bed file is in one-based coords and contains both endpoints of each span.\n";
    print "Also output are: a file of indels, a file of substitutions, and a file of novel splice forms.\n";
    print "\n   options:\n";
    print "      -numgenes n     : Choose n>=1 genes at random from a master pool of gene models (default n = 100000).\n";
    print "      -error x        : Set the error rate for generating a wrong base\n                        anywhere in the read to 0<=x<=1 (default x = 0.005).\n";
    print "      -indelfreq x    : Set indel rate to 0<=x<1 (default x = 0.0005).\n";
    print "      -nalt n         : Set the number of novel splice forms per gene to n>1 (default n = 2).\n";
    print "      -palt x         : Set the percentage of signal coming from novel splice\n                        forms to 0<=x<=1 (default x = 0.2).\n";
    print "      -readlength n   : Set readlength to n>0 bases (default n = 100).\n";
    print "      -sn             : Add sequence number to the first column of the bed file.\n";
    print "      -filenamestem x : The stem x will be added to the four default filenames for the four\n                        required files.  Use this if you have made your own config files (see\n                        below about config files and custom config files).\n";
    print "                        E.g. \"simulator_config_geneinfo\" becomes \"simulator_config_geneinfo_x\".\n";
    print "      -subfreq x      : Set substitution rate to 0<=x<1 (default x = 0.001).\n";
    print "      -tlen n         : Set the length of the low quality tail to n bases (default n = 10).\n";
    print "      -tpercent x     : Set the percent of tails that are low quality to 0<=x<=1\n                        (default x = 0).\n";
    print "      -tqual x        : Set quality of the low quality tail to 0<=x<=1 (default x = 0.8).\n";
    print "\n";
    print "This program depends on four files:\n";
    print "  1) simulator_config_geneinfo\n";
    print "  2) simulator_config_geneseq\n";
    print "  3) simulator_config_intronseq\n";
    print "  4) simulator_config_featurequantifications\n\n";
    print "To create such files for a subset of genes use the script:\n";
    print "   - make_config_files_for_subset_of_gene_ids.pl\n";
    print "Run it with no parameters for the usage\n";
    print "To use those config files with this program use the option -filenamestem\n";
    print "\n";
    exit(0);
}

$name = $ARGV[1];

if($name =~ /^(-|_)/) {
    print STDERR "\nError: invalid name '$name' - name cannot start with a dash or an underscore.\n\n";
    exit(0);
}
if($name =~ /[^a-zA-Z0-9_\-\.]/) {
    print STDERR "\nError: name must be only letters, numbers, dashes, underscores and periods.\n\n";
    exit(0);
}
if(!($name =~ /[a-zA-Z0-9_\-\.]/)) {
    print STDERR "\nERROR: invalid name '$name', cannot be empty and must use only letters,\n       numbers, dashes, underscores and periods.\n\n";
    exit(0);
}

$READLENGTH = 100;
$substitutionfrequency = .001;
$indelfrequency = .0005;
$base_error = .005;
$low_qual_tail_length = 10;
$percent_of_tails_that_are_low_qual = 0;
$quality_of_low_qual_tail = .8;
$percent_alt_spliceforms = .2;
$num_alt_splice_forms_per_gene = 2;
$NUMGENES = 100000;
$stem = "";

use Math::Random qw(:all);
$num_reads = $ARGV[0];
if(!($num_reads =~ /^\d+$/) || ($num_reads <= 0)) {
    print STDERR "\nError: number of reads must be an integer, not '$ARGV[0]'.\n\n";
    exit(0);
}
$sum_of_gene_intensities = 0;
$sum_of_intron_intensities = 0;
$genecnt=0;
$exoncount_total=0;
$introncount_total=0;

$seq_num_in_bedfile = "false";
for($i=2; $i<@ARGV; $i++) {
    $option_recognized = 0;
    if($ARGV[$i] eq "-sn") {
	$seq_num_in_bedfile = "true";
	$option_recognized = 1;
    }
    if($ARGV[$i] eq "-filenamestem") {
	$i++;
	$stem = $ARGV[$i];
	$option_recognized = 1;
	if($stem =~ /^(-|_)/) {
	    print STDERR "\nError: -filenamestem cannot start with a dash or an underscore.\n\n";
	    exit(0);
	}
	if($stem =~ /[^a-zA-Z0-9._-]/) {
	    print STDERR "\nError: -filenamestem must be only letters, numbers, dashes, periods and underscores.\n\n";
	    exit(0);
	}

	$simulator_config_geneinfo = $simulator_config_geneinfo . "_$stem";
	$simulator_config_featurequantifications = $simulator_config_featurequantifications . "_$stem";
	$simulator_config_geneseq = $simulator_config_geneseq . "_$stem";
	$simulator_config_intronseq = $simulator_config_intronseq . "_$stem";
    }
    if($ARGV[$i] eq "-tqual") {
	$i++;
	$quality_of_low_qual_tail = $ARGV[$i];
	$option_recognized = 1;
	if(!($quality_of_low_qual_tail =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -tqual has to be non-negative and no more than one.\n\n";
	    exit(0);
	}
	$quality_of_low_qual_tail = $quality_of_low_qual_tail + 0;
	if($quality_of_low_qual_tail > 1 || $quality_of_low_qual_tail < 0) {
	    print STDERR "\nError: -tqual has to be non-negative and no more than one.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-tpercent") {
	$i++;
	$percent_of_tails_that_are_low_qual = $ARGV[$i];
	$option_recognized = 1;
	if(!($percent_of_tails_that_are_low_qual =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -tpercent has to be non-negative and no more than one.\n\n";
	    exit(0);
	}
	$percent_of_tails_that_are_low_qual = $percent_of_tails_that_are_low_qual + 0;
	if($percent_of_tails_that_are_low_qual > 1 || $percent_of_tails_that_are_low_qual < 0) {
	    print STDERR "\nError: -tpercent has to be non-negative and no more than one.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-tlen") {
	$i++;
	$low_qual_tail_length = $ARGV[$i];
	$option_recognized = 1;
	if(!($low_qual_tail_length =~ /^\d+$/)) {
	    print STDERR "\nError: -teln must be a positive integer.\n\n";
	    exit(0);
	}
	if($low_qual_tail_length < 1) {
	    print STDERR "\nError: -teln must be a positive integer.\nIf you want there to be no tail error don't set -tlen\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-error") {
	$i++;
	$base_error = $ARGV[$i];
	$option_recognized = 1;
	if(!($base_error =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -error has to be non-negative and less than one\n\n";
	    exit(0);
	}
	$base_error = $base_error + 0;
	if($base_error >= 1 || $base_error < 0) {
	    print STDERR "\nError: -error has to be non-negative and less than one\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-nalt") {
	$i++;
	$num_alt_splice_forms_per_gene = $ARGV[$i];
	$option_recognized = 1;
	if(!($num_alt_splice_forms_per_gene =~ /^\d+$/)) {
	    print STDERR "\nError: -nalt has to be a postive integer.\nIf you want there to be no alternate splice forms don't set -nalt, set -palt 0\n\n";
	    exit(0);
	}
	if($num_alt_splice_forms_per_gene < 1) {
	    print STDERR "\nError: -nalt has to be a postive integer.\nIf you want there to be no alternate splice forms don't set -nalt, set -palt 0\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-numgenes") {
	$i++;
	$NUMGENES = $ARGV[$i];
	$option_recognized = 1;
	if(!($NUMGENES =~ /^\d+$/)) {
	    print STDERR "\nError: -numgenes has to be a postive integer.\n\n";
	    exit(0);
	}
	if($NUMGENES < 1) {
	    print STDERR "\nError: -numgenes has to be a postive integer.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-palt") {
	$i++;
	$percent_alt_spliceforms = $ARGV[$i];
	$option_recognized = 1;
	if(!($percent_alt_spliceforms =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -palt has to be non-negative and no more than one.\n\n";
	    exit(0);
	}
	$percent_alt_spliceforms = $percent_alt_spliceforms + 0;
	if($percent_alt_spliceforms > 1 || $percent_alt_spliceforms < 0) {
	    print STDERR "\nError: -palt has to be non-negative and no more than one.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-indelfreq") {
	$i++;
	$indelfrequency = $ARGV[$i];
	$option_recognized = 1;
	if(!($indelfrequency =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -indelfreq has to be strictly between 0 and 1.\n\n";
	    exit(0);
	}
	$indelfrequency = $indelfrequency + 0;
	if($indelfrequency >= 1 || $indelfrequency < 0) {
	    print STDERR "\nError: -indelfreq has to be strictly between 0 and 1.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-subfreq") {
	$i++;
	$substitutionfrequency = $ARGV[$i];
	$option_recognized = 1;
	if(!($substitutionfrequency =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -subfreq has to be less than one non-negative.\n\n";
	    exit(0);
	}
	$substitutionfrequency = $substitutionfrequency + 0;
	if($substitutionfrequency >= 1 || $substitutionfrequency < 0) {
	    print STDERR "\nError: -subfreq has to be less than one non-negative.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-readlength") {
	$i++;
	$READLENGTH = $ARGV[$i];
	$option_recognized = 1;
	if(!($READLENGTH =~ /^\d+$/)) {
	    print STDERR "\nError: -readlength has to be a postive integer.\n\n";
	    exit(0);
	}
	$READLENGTH = $READLENGTH + 0;
	if($READLENGTH < 1) {
	    print STDERR "\nError: -readlength has to be a postive integer.\n\n";
	    exit(0);
	}
    }
    if($option_recognized == 0) {
	print STDERR "\nError: option $ARGV[$i] not recognized.\n\n";
	exit(0);
    }
}

if($READLENGTH <= $low_qual_tail_length && $percent_of_tails_that_are_low_qual > 0) {
    print STDERR "\nERROR: low quality tail length must be less than the readlength.\n\n";
    exit(0);
}
if(!($stem =~ /\S/)) {
    # Here construct random set of $NUMGENES genes by calling make_config...pl script on master-list files

    $simulator_config_geneinfo = "simulator_config_geneinfo_temp";
    $simulator_config_featurequantifications = "simulator_config_featurequantifications_temnp";
    $simulator_config_geneseq = "simulator_config_geneseq_temp";
    $simulator_config_intronseq = "simulator_config_intronseq_temp";
}

$bedfilename = "simulated_reads_$name" . ".bed";
$fafilename = "simulated_reads_$name" . ".fa";
$substitutionsfilename = "simulated_reads_substitutions_$name" . ".txt";
$indelsfilename = "simulated_reads_indels_$name" . ".txt";

open(SIMBEDOUT, ">$bedfilename");
open(SIMFAOUT, ">$fafilename");
open(SIMSUBSOUT, ">$substitutionsfilename");
open(SIMINDELSOUT, ">$indelsfilename");

open(INFILE, $simulator_config_geneinfo);

while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    @b = split(/::::/,$a[7]);
    for($i=0; $i<@b; $i++) {
	$b[$i] =~ s/\(.*\)//;
	$a[5] =~ s/,\s+$//;
	$a[6] =~ s/,\s+$//;
	$starts{$b[$i]} = $a[5];
	$ends{$b[$i]} = $a[6];
	$chr{$b[$i]} = $a[0];
	$strand{$b[$i]} = $a[1];
    }
}
close(INFILE);

open(INFILE, $simulator_config_featurequantifications);

while($line = <INFILE>) {
    chomp($line);
    if($line =~ /(exon)/) {
	@a = split(/\t/,$line);
	@a2 = @{$exon2gene{$a[1]}};
	$N = @a2 + 0;
	$exon2gene{$a[1]}[$N] = $geneid;
	$gene2exon{$geneid}[$exoncnt] = $a[1];
	$exoncnt++;
    }
    if($line =~ /(intron)/) {
	@a = split(/\t/,$line);
	@a2 = @{$intron2gene{$a[1]}};
	$N = @a2 + 0;
	$intron2gene{$a[1]}[$N] = $geneid;
	$gene2intron{$geneid}[$introncnt] = $a[1];
	$intron2intensity{$a[1]} = $a[2];
	$introncnt++;
	$gene2introncnt{$geneid}{$a[1]} = $introncnt;
	##print SIMFAOUT "introncnt = $introncnt\n";
	#print SIMFAOUT "intron = $a[1]\n";
    }
    if($line =~ /--------/) {
	$line = <INFILE>;
	chomp($line);
	$line =~ s/\(.*//;
	$geneid = $line;
	$exoncnt=0;
	$introncnt=0;
	$line = <INFILE>;
	$line = <INFILE>;
	chomp($line);
	@a = split(/\t/,$line);
	$gene_intensity[$genecnt] = $a[2];
	$sum_of_gene_intensities = $sum_of_gene_intensities + $gene_intensity[$genecnt];
	$genes[$genecnt] = $geneid;
	$genecnt++;
    }
}
close(INFILE);
$genecnt--;
$numgenes = $genecnt;
print "$numgenes genes total\n";
print "sum_of_gene_intensities = $sum_of_gene_intensities\n";

# making alternate (unknown) splice forms
for($i=0; $i<$numgenes; $i++) {
    # $i is the original gene id
    # $genecnt is the same gene as $i but will be modified to be an alternate splice form
    $geneid = $genes[$i];
    @a = @{$gene2exon{$geneid}};
    $exoncnt = @a;
    $original_gene_intensity = $gene_intensity[$i];
    # the following adjusts the original intensity because now it will be split between
    # the normal and alt spice form
    $gene_intensity[$i] = $original_gene_intensity * (1-$percent_alt_spliceforms);
    $genes2[$i] = $geneid;
    for($j=0; $j<$num_alt_splice_forms_per_gene; $j++) {
	$geneid_x = $geneid . "_$j";
	$genes2[$genecnt+$j*$numgenes] = $geneid_x;
	$gene_intensity[$genecnt+$j*$numgenes] = $original_gene_intensity * $percent_alt_spliceforms/$num_alt_splice_forms_per_gene;
	$exoncnt_x = 0;
	$str = $starts{$geneid};
	$str =~ s/^\s*,\s*//;
	$str =~ s/\s*,\s*$//;
	@S = split(/,/,$str);
	$str = $ends{$geneid};
	$str =~ s/^\s*,\s*//;
	$str =~ s/\s*,\s*$//;
	@E = split(/,/,$str);
	$FLAG1 = 0;
	$FLAG2 = 0;
	while($FLAG1 == 0 || $FLAG2 == 0) {
	    $FLAG1 = 0;
	    $FLAG2 = 0;
	    $starts_new = "";
	    $ends_new = "";
	    $exoncnt_x = 0;
	    delete($gene2exon{$geneid_x});
	    for($k=0; $k<$exoncnt; $k++) {
		$flip = int(rand(2));
		if($flip == 1) {
		    if($strand{$geneid} eq "+") {
			$starts_new = $starts_new . "$S[$k],";
			$ends_new = $ends_new . "$E[$k],";
		    }
		    if($strand{$geneid} eq "-") {
			$starts_new = "$S[$numexons-$k-1]," . $starts_new;
			$ends_new = "$E[$numexons-$k-1]," . $ends_new;
		    }
		    $gene2exon{$geneid_x}[$exoncnt_x] = $gene2exon{$geneid}[$k];
		    $exoncnt_x++;
		    $FLAG1 = 1; # This makes sure we kept at least one exon
		}
		else {
		    $FLAG2 = 1; # This makes sure we skipped at least one exon
		}
	    }
	    if($exoncnt == 1) {
		$FLAG2 = 1;
	    }
	}
	$starts{$geneid_x} = $starts_new;
	$ends{$geneid_x} = $ends_new;
	$strand{$geneid_x} = $strand{$geneid};
	$chr{$geneid_x} = $chr{$geneid};
    }
    $genecnt++;
}
@genes = @genes2;
$genecnt = @genes;

$alttranscriptsfilename = "transcripts_$name" . ".txt";
open(ALTTRANSCRIPTS, ">$alttranscriptsfilename");
for($i=0; $i<@genes; $i++) {
    $geneid = $genes[$i];
    @a = @{$gene2exon{$geneid}};
    $exoncnt = @a;
    $str1 = $starts{$geneid};
    $str2 = $ends{$geneid};
    $str3 = $strand{$geneid};
    $str4 = $chr{$geneid};
    print ALTTRANSCRIPTS "---------\ngenes[$i] = $geneid\n";
    print ALTTRANSCRIPTS "starts = $str1\n";
    print ALTTRANSCRIPTS "ends = $str2\n";
    print ALTTRANSCRIPTS "strand = $str3\n";
    print ALTTRANSCRIPTS "chr = $str4\n";
    for($j=0; $j<$exoncnt; $j++) {
	print ALTTRANSCRIPTS "$a[$j]\n";
    }
}
close(ALTTRANSCRIPTS);

if($sum_of_gene_intensities <= 0) {
    print STDERR "\nERROR: None of the genes have expression level above zero,\ncheck your feature quantification config file.\n\n";
    exit(0);
}

for($i=0; $i<$genecnt; $i++) {
    $gene_density[$i] = $gene_intensity[$i] / $sum_of_gene_intensities;
    $gene_distribution[$i] = $gene_density[$i];
}
print "gene_distribution[0] = $gene_distribution[0]\n";
for($i=1; $i<$genecnt; $i++) {
    $gene_distribution[$i] = $gene_distribution[$i] + $gene_distribution[$i-1];
    print "gene_distribution[$i] = $gene_distribution[$i]\n";
}

$introncount_total=0;
$sum_of_intron_intensities = 0;
foreach $intron (keys %intron2gene) {
    $introns[$introncount_total] = $intron;
    $intron_intensity[$introncount_total] = $intron2intensity{$intron};
    $sum_of_intron_intensities = $sum_of_intron_intensities + $intron_intensity[$introncount_total];
    $introncount_total++;
}

print "sum_of_intron_intensities = $sum_of_intron_intensities\n";
# The following formula is based on the fact that introns are not treated as isolated but are padded with
# the rest of the exons.
$intron_freq = (2 * $sum_of_intron_intensities) / ($sum_of_gene_intensities + $sum_of_intron_intensities);
$iftemp = $sum_of_intron_intensities / ($sum_of_gene_intensities + $sum_of_intron_intensities);

print "intron frequency: $iftemp\n";
print "padded intron freq = $intron_freq\n";

for($i=0; $i<$introncount_total; $i++) {
    $intron_density[$i] = $intron_intensity[$i] / $sum_of_intron_intensities;
    $intron_distribution[$i] = $intron_density[$i];
}
for($i=1; $i<$introncount_total; $i++) {
    $intron_distribution[$i] = $intron_distribution[$i] + $intron_distribution[$i-1];
}

#for($i=0; $i<$introncount_total; $i++) {
#    print "intron_distribution[$i] = $intron_distribution[$i]\t$introns[$i]\t$intron_intensity[$i]\n";
#}


open(INFILE, $simulator_config_geneseq);

 # this file has seqs on one line and minus strand seqs revserse complemented
 # header line looks like this: >NM_175370:chr1:58714964-58752833_-
$CNT=0;
print "\nReading gene sequences\n";
while($line = <INFILE>) {
    $CNT++;
    if($CNT % 10000 == 0) {
	print "$CNT\n";
    }
    chomp($line);
    $line =~ /:(.*):\d+-\d+/;
    $chr = $1;
    $line =~ s/:.*//;
    $line =~ s/>//;
    $gene_id = $line;
    $SEQ = <INFILE>;
    chomp($SEQ);
    $seq{$gene_id} = $SEQ;
    @S = split(/,/, $starts{$gene_id});
    @E = split(/,/, $ends{$gene_id});
    $offset = 0;
    for($exon=0; $exon<@S; $exon++) {
	$S[$exon]++;
	$len = $E[$exon] - $S[$exon] + 1;
	$EXONSEQ = substr($SEQ,$offset,$len);
	$exonname = $chr . ":" . $S[$exon] . "-" . $E[$exon];
	$exonseq{$exonname} = $EXONSEQ;
	$offset = $offset + $len;
    }
}
close(INFILE);

open(INFILE, $simulator_config_intronseq);

print "\nReading intron sequences\n";
$flag = 0;
$line = <INFILE>;
$CNT=0;
while($flag == 0) {
    $CNT++;
    if($CNT % 10000 == 0) {
	print "$CNT\n";
    }
    $line =~ /^>(.*)/;
    $intron = $1;
    $seq = "";
    $line = <INFILE>;
    until($line =~ /^>/ || $flag == 1) {
	chomp($line);
	if($line !~ /\S/) {
	    $flag = 1;
	}
	else {
	    if($line !~ /^>/) {
		$seq = $seq . $line;
	    }
	}
	$line = <INFILE>;
    }
    $intronseq{$intron} = $seq;
}
close(INFILE);

# The following puts substitutions and indels into each exon
foreach $exon (keys %exon2gene) {
    $try_cnt = 0;
    $exon =~ /^(.*):(\d+)-(\d+)/;
    $chr = $1;
    $start = $2;
    $end = $3;
    $length = $end - $start + 1;
    $num_substitutions = random_binomial(1, $length, $substitutionfrequency);
    undef %substitutions_locs;
    undef @indels_temp;
    $SEQ = $exonseq{$exon};
    for($i=0; $i<$num_substitutions; $i++) {
	$LOC = int(rand($length)) + 1;
	while($substitutions_locs{$LOC}+0>0) {
	    $LOC = int(rand($length)) + 1;
	}
	$substitutions_locs{$LOC}++;
	$substitutions{$exon}[$i] = $LOC;
	$orig = substr($SEQ,$LOC-1,1);
	$B = getrandombase();
	while($B eq $orig) {
	    $B = getrandombase();
	}
	$Z = substr($SEQ,$LOC-1,1,$B);
	$exonseq{$exon} = $SEQ;
	$C = $start + $LOC - 1;
	print SIMSUBSOUT "$exon\t$C\t$Z->$B\n";
    }

    $num_indels = random_binomial(1, $length, $indelfrequency);
    $flag = 0;
    while($flag == 0) {  # the following gets the indel locations and makes sure
                         # indels are at least two bases apart and are in different
                         # places from the substitutions
	$flag = 1;
	undef %indels_locs_temp;
	undef %indels_locs;
	for($i=0; $i<$num_indels; $i++) {
	    $LOC = int(rand($length)+1);
	    while(($substitutions_locs{$LOC}+0>0) && ($indels_locs_temp{$LOC}+0>0)) {
		$LOC = int(rand($length)+1);
	    }
	    $indels_locs_temp{$LOC}++;
	}
	foreach $LOC (keys %indels_locs_temp) {
	    if($indels_locs_temp{$LOC} > 0) {
		$indels_locs{$LOC}++;
	    }
	}
	foreach $LOC1 (keys %indels_locs) {
	    foreach $LOC2 (keys %indels_locs) {
		$X = $LOC1 - $LOC2;
		if($X < 2 && $X > -2 && $X != 0) {
		    $flag = 0;
		}
	    }	    
	}	
	$indelcounter=0;
	foreach $LOC (sort {$b<=>$a} keys %indels_locs) {
	    $indellength = int(random_exponential(1, 1));
	    while($indellength < 1) {
		$indellength = int(random_exponential(1, 1));
	    }
	    $flip = int(rand(2));
	    if($flip == 1) {                                       # INSERTION
		$insert = "";
		for($i=0; $i<$indellength; $i++) {
		    $insert = $insert . getrandombase();
		}
		$indels_temp[$indelcounter][0] = $LOC;
		$indels_temp[$indelcounter][1] = $indellength;
		$indels_temp[$indelcounter][2] = $insert;
		$indelcounter++;
	    }
	    else {                                                 # DELETION
		# have to make sure we don't delete a substitution or part of an insertion or overlap
                # two deletions, and make sure it doesn't overrun the end of the sequence
		foreach $LOC2 (keys %substitutions_locs) {
		    if($LOC2 <= $LOC + $indellength && $LOC2 > $LOC) {
			$flag = 0;
		    }
		}
		foreach $LOC2 (keys %indels_locs) {
		    if($LOC2 <= $LOC + $indellength && $LOC2 > $LOC) {
			$flag = 0;
		    }
		    if($LOC + $indellength >= $length) {
			$flag = 0;
		    }		    
		}
		$indels_temp[$indelcounter][0] = $LOC;
		$indels_temp[$indelcounter][1] = -1 * $indellength;
		$indelcounter++;
	    }
	}
	if($flag == 1) {
	    for($j=0;$j<$indelcounter;$j++) {
		$LOC = $indels_temp[$j][0];
		$indellength = $indels_temp[$j][1];
		$insert = $indels_temp[$j][2];
		$indels{$exon}[$j][0] = $LOC;  # This is the location within the exon, where the
                                               # first base in the exon is location 1.
		$indels{$exon}[$j][1] = $indellength;
		$indels{$exon}[$j][2] = $insert;
		if($indellength > 0) {
		    $Z = substr($SEQ,$LOC,0,$insert);
		    $exonseq{$exon} = $SEQ;
		    print SIMINDELSOUT "$exon\t$LOC\t$indellength\t$insert\n";
		}
		if($indellength < 0) {
		    $l = -1 * $indellength;
		    $Z = substr($SEQ,$LOC,$l,"");
		    $exonseq{$exon} = $SEQ;
		    print SIMINDELSOUT "$exon\t$LOC\t$indellength\t$Z\n";
		}
	    }
	} else {
	    $try_cnt++;
	    if($try_cnt > 500) {
		$num_indels = int($num_indels/2);
		$try_cnt=0;
	    }
	}	
    }
}
close(SIMINDELSOUT);
close(SIMSUBSOUT);

# Now need to assemble the exons back into genes
print "\nAssembling exons\n";
foreach $geneid (keys %gene2exon) {
    @a = @{$gene2exon{$geneid}};
    $numexons = @a;
    $SEQ = "";
    $indelcounter=0; # counts the number of indels per *gene*
    $offset = 0;
    for($j=0; $j<$numexons; $j++) {
	if($strand{$geneid} eq "+") {
	    $i = $j;
	} else {
	    $i = $numexons - $j - 1;
	}
	$exon = $a[$i];
	$exon =~ /^[^:]+:(\d+)-(\d+)/;
	$exonstart = $1;
	$exonend = $2;
	@a2 = @{$indels{$exon}};
        # need to make hash that associates *genes* to indels in gene specific coords (start of gene = 1)
	for($k=0;$k<@a2;$k++) {
	    $location_in_gene = $offset + $a2[$k][0];
	    $gene2indel_temp{$geneid}[$indelcounter][0] = $location_in_gene;
	    $gene2indel_temp{$geneid}[$indelcounter][1] = $a2[$k][1];  # the length of the indel
	    $indelcounter++;
	}
	$offset = $offset + $exonend - $exonstart + 1;
	$SEQ = $SEQ . $exonseq{$a[$i]};
    }
    $seq{$geneid} = $SEQ;
    $FLAG = 0;
    @a = @{$gene2indel_temp{$geneid}};
    while($FLAG == 0) {
	$FLAG = 1;
	for($i=0; $i<@a-1; $i++) {
	    if($a[$i][0] > $a[$i+1][0]) {
		$temp = $a[$i][0];
		$a[$i][0] = $a[$i+1][0];
		$a[$i+1][0] = $temp;
		$temp = $a[$i][1];
		$a[$i][1] = $a[$i+1][1];
		$a[$i+1][1] = $temp;
		$FLAG = 0;
	    }
	}
    }
    for($i=0; $i<@a; $i++) {
	$gene2indel{$geneid}[$i][0] = $a[$i][0];
	$gene2indel{$geneid}[$i][1] = $a[$i][1];
    }
}

print "\nGenerating reads\n";
$CNT=1;
while( 1 == 1) {
    if($CNT % 50000 == 0) {
	print "$CNT reads done\n";
    }
    $R = rand(1);
    if($R < $intron_freq) {
	# pick an intron at random
	$R = rand(1);
	$c = 0;
	while($intron_distribution[$c] < $R && $c<(@intron_distribution-1)) {
	    $c++;
	}
	$INTRON2 = $introns[$c];
	# get a gene at random from the genes containing this intron
	@i2g = @{$intron2gene{$INTRON2}};
	$n = @i2g;
	$R2 = int(rand($n));
	$GENE = $i2g[$R2];
	@g2i = @{$gene2intron{$GENE}};
	$num_introns = @g2i;
	# get the intron number for this intron in this gene
	$intron_num = $gene2introncnt{$GENE}{$INTRON2};

	$SEQ = $seq{$GENE};
	$seqlength = length($SEQ);
	undef %geneWithIntron2indel;
	$SEQ = getpaddedintron($GENE, $intron_num);
	$seqlength = length($SEQ);

	@INDELS = @{$geneWithIntron2indel{$GENE}{$INTRON2}};
# the following fixes the starts/ends so they reflects the retained intron:
	$STARTS2 = $starts{$GENE};
	$STARTS2 =~ s/,\s*$//;
	$STARTS2 =~ s/^\s*,//;
	@S1 = split(/,/, $STARTS2);
	for($ii=0; $ii<@S1; $ii++) {
	    $S1[$ii]++;
	}
	$ENDS2 = $ends{$GENE};
	$ENDS2 =~ s/,\s*$//;
	$ENDS2 =~ s/^\s*,//;
	@E1 = split(/,/, $ENDS2);
	$STARTS = $S1[0];
	$ENDS = "";
	for($p=1; $p<@S1; $p++) {
	    if($strand{$GENE} eq "+") {
		if(!($intron_num == $p)) {
		    $STARTS = $STARTS . ",$S1[$p]";
		    $ENDS = $ENDS . ",$E1[$p-1]";
		}
	    }
	    if($strand{$GENE} eq "-") {
		if(!($intron_num == @S1 - $p)) {
		    $STARTS = $STARTS . ",$S1[$p]";
		    $ENDS = $ENDS . ",$E1[$p-1]";
		}
	    }
	}
	$ENDS = $ENDS . ",$E1[@S1-1]";
	$STARTS =~ s/,\s*$//;
	$STARTS =~ s/^\s*,//;
	$ENDS =~ s/,\s*$//;
	$ENDS =~ s/^\s*,//;
	@S1 = split(/,/, $STARTS);
	$STARTS2 = "";
	for($ii=0; $ii<@S1; $ii++) {
	    $S1[$ii]--;
	    $STARTS2 = $STARTS2 . "$S1[$ii],";
	}
	$STARTS2 =~ s/,\s*$//;
	$STARTS2 =~ s/^\s*,//;
	$STARTS = $STARTS2;

	@S1 = split(/,/, $starts{$GENE});
	@E1 = split(/,/, $ends{$GENE});

	@S1 = split(/,/, $STARTS);
	@E1 = split(/,/, $ENDS);

    }
    else {
	$R = rand(1);
	$c = 0;
	while($gene_distribution[$c] < $R && $c<(@gene_distribution-1)) {
	    $c++;
	}
	$GENE = $genes[$c];
	@INDELS = @{$gene2indel{$GENE}};
	$SEQ = $seq{$GENE};
	$STARTS = $starts{$GENE};
	$ENDS = $ends{$GENE};
    }
    $return_vector_ref = getreads($SEQ, \@INDELS, $STARTS, $ENDS, $CNT);

    @return_vector = @{$return_vector_ref};
    $fa = $return_vector[0];
    $bed = $return_vector[1];
    $mergedcoords = $return_vector[2];

    if($fa ne "none") {
	@F = split(/\n/,$fa);
	print SIMFAOUT "$F[0]\n";
	$stemp = &add_sequencing_error($F[1], $base_error);
	$stemp2 = &add_error_to_tail($stemp, $low_qual_tail_length, $percent_of_tails_that_are_low_qual, $quality_of_low_qual_tail);
	print SIMFAOUT "$stemp2\n";
	print SIMFAOUT "$F[2]\n";

	$stemp = &add_sequencing_error($F[3], $base_error);
	$stemp2 = &add_error_to_tail($stemp, $low_qual_tail_length, $percent_of_tails_that_are_low_qual, $quality_of_low_qual_tail);
	print SIMFAOUT "$stemp2\n";

	print SIMBEDOUT $bed;
	
	$CNT++;
	if($CNT > $num_reads) {
	    close(SIMBEDOUT);
	    close(SIMFAOUT);
	    exit(0);
	}
    }
}

sub getreads () {
    ($SEQ, $INDELS_ref, $STARTS, $ENDS, $CNT) = @_;

#    The following depends on:
#    1) %seq which maps gene ids to gene sequence
#    2) @INDELS a 2-d array, indel loc is $INDELS[$LOC][0], length is $INDELS[$LOC][1], sorted by $LOC
#    3) $starts and $ends, for the exon start/end positions, zero-based half open
#    4) The current read counter $CNT

    $seqlength = length($SEQ);
    if($seqlength < $READLENGTH) {
	$return_vector[0] = "none";
	$return_vector[1] = "none";
	$return_vector[2] = "none";
	return \@return_vector;
    }
    $fragmentlength = int(random_normal() * $READLENGTH/3 + $READLENGTH*3);
    if($fragmentlength < $READLENGTH) {
	$fragmentlength = $READLENGTH;
    }

    if($fragmentlength > $seqlength) {
	$fragmentlength = $seqlength;
    }
    if($seqlength > 2 * $READLENGTH) {
	$cutpoint = int(rand($seqlength-$fragmentlength+1) + $fragmentlength/2);
	$repeat = rand($seqlength / $cutpoint);
	$flip = int(rand(2));
	if($flip == 0) {
	    #print SIMFAOUT "flip = $flip\n";
	    $end = $cutpoint;
	    if($cutpoint >= $fragmentlength) {
		$start = $cutpoint - $fragmentlength + 1;
	    }
	    else {
		$start = 1;
	    }
	}
	if($flip == 1) {
	    $start = $cutpoint;
	    if($seqlength >= $cutpoint + $fragmentlength - 1) {
		$end = $cutpoint + $fragmentlength - 1;
	    }
	    else {
		$end = $seqlength;
	    }
	}
    }
    else {
	$start = 1;
	$end = $seqlength;
    }
    if($start == 1 && $end < $seqlength - $READLENGTH) {
	$flip2 = int(rand(4));
	if($flip2 == 1) {
	    $FF = int(rand($READLENGTH));
	    $start = $start + $FF;
	    $end = $end + $FF;
	}
    }
    if($start == $seqlength && $start > $READLENGTH) {
	$flip2 = int(rand(4));
	if($flip2 == 1) {
	    $FF = int(rand($READLENGTH));
	    $start = $start - $FF;
	    $end = $end - $FF;
	}
    }

    $start_forward = $start;
    $len = $end - $start + 1;
    if($len < $READLENGTH) {
	next;
    }

    $forward_read = substr($SEQ, $start_forward-1, $READLENGTH);
    $start_reverse = $end - $READLENGTH + 1;
    $reverse_read = substr($SEQ, $start_reverse-1, $READLENGTH);

    $rev_rev = reversecomplement($reverse_read);
    $fa = ">seq.$CNT" . "a\n$forward_read\n>seq.$CNT" . "b\n$rev_rev\n";

    # THE FOLLOWING GETS THE TRUE COORDINATES OF THE READ:

    #    ****   FORWARD READ   ****

    # fragment FORWARD read by its deletions and feed each piece separately to getcoords(),
    # then contcatenate coords together, also must correct the start and length of the FORWARD
    # read for the insertions (if any)

    $readlength = $READLENGTH;
    $checker_forward = 0;

    # this loop adjusts the start of the read for insertions/deletions upstream of the start
    for($ind=0; $ind<@INDELS; $ind++) {  # indels are sorted by location
	if($INDELS[$ind][0] < $start_forward) {
	    $start_forward = $start_forward - $INDELS[$ind][1]; 
	}
    }
    $end_forward = $start_forward + $readlength - 1;
    
    $offset = 0;
    $coords = "";
    # this loop does the fragmenting and feeding of each piece to getcoords()
    for($ind=0; $ind<@INDELS; $ind++) {  # indels are sorted by location
	if($start_forward <= $INDELS[$ind][0] && $INDELS[$ind][0] < ($start_forward+$readlength-1)) {
	    if($INDELS[$ind][1] > 0) {  # insertion w.r.t. reference
                # in case insertion goes beyond end of read, don't want
                # to overcorrect, the following 'if' takes care of that
		if( ($INDELS[$ind][0]+$INDELS[$ind][1]) > ($start_forward+$readlength-1) ) {
		    $adjustement_factor = (($INDELS[$ind][0]+$INDELS[$ind][1])-($start_forward+$readlength-1));
		    $readlength = $readlength + $adjustement_factor;
		    $end_forward = $end_forward + $adjustement_factor;
		    $checker_forward = 1;
		}
	    }
	    if($INDELS[$ind][1] < 0) {  # deletion w.r.t. reference
		$fragment_start = $start_forward + $offset;
		$fragment_length = $INDELS[$ind][0] - $fragment_start + 1;
		$coords = $coords . ", " . getcoords($fragment_start, $fragment_length, $STARTS, $ENDS);
		$offset = $offset + $fragment_length - $INDELS[$ind][1];
	    }
	    $readlength = $readlength - $INDELS[$ind][1];
	    $end_forward = $end_forward - $INDELS[$ind][1];
	}
    }
    
    $fragment_start = $start_forward + $offset;
    $fragment_length = $end_forward - $fragment_start + 1;
    $coords = $coords . ", " . getcoords($fragment_start, $fragment_length, $STARTS, $ENDS);
    $coords =~ s/^\s*,\s*//;
    $coords =~ s/\s*,\s*$//;
    $coords1 = $coords;
    $coords1 =~ s/.*\t//;

    #    ****   REVERSE READ   ****

    # fragment REVERSE read by its deletions and feed each piece separately to getcoords(),
    # then contcatenate coords together, also must correct the start and length of the REVERSE
    # read for the insertions (if any)

    $readlength = $READLENGTH;
    $checker_reverse = 0;

    # this loop adjusts the start of the read for insertions/deletions upstream of the start
    for($ind=0; $ind<@INDELS; $ind++) {  # indels are sorted by location
	if($INDELS[$ind][0] < $start_reverse) {
	    $start_reverse = $start_reverse - $INDELS[$ind][1]; 
	}
    }
    $end_reverse = $start_reverse + $readlength - 1;
    
    $offset = 0;
    $coords = "";
    # this loop does the fragmenting and feeding of each piece to getcoords()
    for($ind=0; $ind<@INDELS; $ind++) {  # indels are sorted by location
	if($start_reverse <= $INDELS[$ind][0] && $INDELS[$ind][0] < ($start_reverse+$readlength-1)) {
	    if($INDELS[$ind][1] > 0) {
		if( ($INDELS[$ind][0]+$INDELS[$ind][1]) > ($start_reverse+$readlength-1) ) {
		    $adjustement_factor = (($INDELS[$ind][0]+$INDELS[$ind][1])-($start_reverse+$readlength-1));
		    $readlength = $readlength + $adjustement_factor;
		    $end_reverse = $end_reverse + $adjustement_factor;
		    $checker_reverse = 1;
		}
	    }
	    if($INDELS[$ind][1] < 0) {
		$fragment_start = $start_reverse + $offset;
		$fragment_length = $INDELS[$ind][0] - $fragment_start + 1;
		$coords = $coords . ", " . getcoords($fragment_start, $fragment_length, $STARTS, $ENDS);
		$offset = $offset + $fragment_length - $INDELS[$ind][1];
	    }
	    $readlength = $readlength - $INDELS[$ind][1];
	    $end_reverse = $end_reverse - $INDELS[$ind][1];
	}
    }
    
    $fragment_start = $start_reverse + $offset;
    $fragment_length = $end_reverse - $fragment_start + 1;
    $coords = $coords . ", " . getcoords($fragment_start, $fragment_length, $STARTS, $ENDS);
    $coords =~ s/^\s*,\s*//;
    $coords =~ s/\s*,\s*$//;
    $coords2 = $coords;
    $coords2 =~ s/.*\t//;
    $mergedcoords = merge($coords1, $coords2);

    @D = split(/, /,$mergedcoords);
    $bed = "";
    for($dd=0; $dd<@D; $dd++) {
	$D[$dd] =~ s/-/\t/;
	$str = "$chr{$GENE}\t$D[$dd]\t+\n";
	if($str =~ /^\S+\t\d+\t\d+\t\+$/) {
	    if($seq_num_in_bedfile eq "true") {
		$bed = $bed . "seq.$CNT\t$str";
	    } else {
		$bed = $bed . "$str";
	    }
	}
	else {
	    print "-------\nstr = $str (something is wrong)\n";
	    print "readlength = $readlength\n";
	    print "coords1 = $coords1\n";
	    print "coords2 = $coords2\n";
	    print "mergedcoords = $mergedcoords\n";
	    print "3:GENE = $GENE\n";
	    print "SEQ = $SEQ\n";
	    print "$fa";
	    print "STARTS = $STARTS\n";
	    print "ENDS = $ENDS\n";
	    for($ind=0; $ind<@INDELS; $ind++) {
		print "INDELS[$ind][0] = $INDELS[$ind][0]\n";
		print "INDELS[$ind][1] = $INDELS[$ind][1]\n";
	    }
	    print "start = $start\n";
	    print "fragment end = $end\n";
	    print "start_forward = $start_forward\n";
	    print "end_forward = $end_forward\n";
	    print "start_reverse = $start_reverse\n";
	    print "end_reverse = $end_reverse\n";
	    print "fragment_start = $fragment_start\n";
	    print "fragment_length = $fragment_length\n";
	    print "fragmentlength = $fragmentlength\n";
	    print "cutpoint = $cutpoint\n";
	    print "seqlength = $seqlength\n";
	    print "flip = $flip\n";
	    print "checker_forward = $checker_forward\n";
	    print "checker_reverse = $checker_reverse\n";
	}
    }
    $return_vector[0] = $fa;
    $return_vector[1] = $bed;
    $return_vector[2] = $mergedspans;

    return \@return_vector;
}

sub merge () {
    ($aspans, $bspans) = @_;
    undef @astarts2;
    undef @aends2;
    undef @bstarts2;
    undef @bends2;

    $aspans =~ /(\d+)$/;
    $aspans_end = $1;
    $bspans =~ /^(\d+)/;
    $bspans_start = $1;
    if(!($aspans_end =~ /\S/) || !($bspans_start =~ /\S/)) {
	print "aspans = $aspans\n";
	print "bspans = $bspans\n";
    }

    @a = split(/, /, $aspans);
    for($i=0; $i<@a; $i++) {
	@b = split(/-/,$a[$i]);
	$astarts2[$i] = $b[0];
	$aends2[$i] = $b[1];
    }
    @a = split(/, /, $bspans);
    for($i=0; $i<@a; $i++) {
	@b = split(/-/,$a[$i]);
	$bstarts2[$i] = $b[0];
	$bends2[$i] = $b[1];
    }
    if($aends2[@aends2-1] + 1 < $bstarts2[0]) {
	$merged_spans = $aspans . ", " . $bspans;
    }
    if($aends2[@aends2-1] + 1 == $bstarts2[0]) {
	$aspans =~ s/-\d+$//;
	$bspans =~ s/^\d+-//;
	$merged_spans = $aspans . "-" . $bspans;
    }
    if($aends2[@aends2-1] + 1 > $bstarts2[0]) {
	$merged_spans = $aspans;
	for($i=0; $i<@bstarts2; $i++) {
	    if($aends2[@aends2-1] >= $bstarts2[$i] && ($aends2[@aends2-1] <= $bstarts2[$i+1] || $i == @bstarts2-1)) {
		$merged_spans =~ s/-\d+$//;
		$merged_spans = $merged_spans . "-" . $bends2[$i];
		for($j=$i+1; $j<@bstarts2; $j++) {
		    $merged_spans = $merged_spans . ", $bstarts2[$j]-$bends2[$j]";
		}
	    }
	}
    }
    return $merged_spans;
}

sub reversecomplement () {
    ($sq) = @_;
    @A = split(//,$sq);
    $rev = "";
    for($i=@A-1; $i>=0; $i--) {
	$flag = 0;
	if($A[$i] eq 'A') {
	    $rev = $rev . "T";
	    $flag = 1;
	}
	if($A[$i] eq 'T') {
	    $rev = $rev . "A";
	    $flag = 1;
	}
	if($A[$i] eq 'C') {
	    $rev = $rev . "G";
	    $flag = 1;
	}
	if($A[$i] eq 'G') {
	    $rev = $rev . "C";
	    $flag = 1;
	}
	if($flag == 0) {
	    $rev = $rev . $A[$i];
	}
    }

    return $rev;
}

sub getcoords () {
    ($readstart, $readlength2, $starts, $ends) = @_;
    $readend = $readstart + $readlength2 - 1;
    $COORDS = "";
    $starts =~ s/,\s*$//;
    $ends =~ s/,\s*$//;
    @S = split(/,/,$starts);
    @E = split(/,/,$ends);
    for($ii=0; $ii<@S; $ii++) {
	$S[$ii]++;
    }
    $cumlength = 0;
    $str = "";
    for($I=0; $I<@S; $I++) {
	$l = $E[$I] - $S[$I] + 1;
	$cumlength = $cumlength + $l;
	$str = $str . "$cumlength, ";
    }
    #print SIMFAOUT "INSIDE GET COORDS:\n";
    #print SIMFAOUT "cumlength: $str\n";
    #print SIMFAOUT "starts = $starts\n";
    #print SIMFAOUT "ends = $ends\n";
    $prefix = 0;
    $exonnum = 0;
    while($prefix < $readstart) {
	$prefix = $prefix + ($E[$exonnum] - $S[$exonnum] + 1);
	$exonnum++;
    }
    $offset1 = $prefix - ($E[$exonnum-1] - $S[$exonnum-1] + 1);
    $firstexon = $exonnum - 1; # the read starts in this exon
    $prefix = 0;
    $exonnum = 0;
    while($prefix < $readend) {
	$prefix = $prefix + ($E[$exonnum] - $S[$exonnum] + 1);
	$exonnum++;
    }
    $offset2 = $prefix - ($E[$exonnum-1] - $S[$exonnum-1] + 1);
    $lastexon = $exonnum - 1; # the read ends in this exon
    if($lastexon > @S) {
	$lastexon = @S - 1;
    }
    for($en = $firstexon; $en<=$lastexon; $en++) {
	if($en == $firstexon) {
	    $Z = $S[$en] + ($readstart - $offset1 - 1);
	    $COORDS = "$Z";
	}
	else {
	    $Z = $S[$en];
	    $COORDS = $COORDS . ", $Z";
	}
	if($en == $lastexon) {
	    $Z = $S[$en] + ($readend - $offset2 - 1);
	    $COORDS = $COORDS . "-$Z";
	}
	else {
	    $COORDS = $COORDS . "-$E[$en]";
	}
    }
    return $COORDS;
}

sub getrandombase () {
    $x = int(rand(4));
    if($x == 0) {
	return "A";
    }
    if($x == 1) {
	return "C";
    }
    if($x == 2) {
	return "G";
    }
    if($x == 3) {
	return "T";
    }
}

# Subroutine to pad each intron with the rest of the exons in the gene
sub getpaddedintron () {
    ($geneid, $intron) = @_;
    $intron--;
    if($intron >= 0) {
	$INTRON = $gene2intron{$geneid}[$intron];
	$INTRON =~ /^[^:]+:(\d+)-(\d+)/;
	$intronstart = $1;
	$intronend = $2;
    }
    else {
	$INTRON = "X";
    }
    @aa = @{$gene2exon{$geneid}};
    $numexons = @aa;
    $SEQ = "";
    $indelcounter=0; # counts the number of indels per *gene*
    $offset = 0;
    for($jj=0; $jj<$numexons; $jj++) {
	if($strand{$geneid} eq "+") {
	    $ii = $jj;
	} else {
	    $ii = $numexons - $jj - 1;
	}
	$exon = $aa[$jj];
	@aa2 = @{$indels{$aa[$ii]}};
	# add to hash that associates the *gene that has the intron $INTRON* to indels in gene w/intron specific coords (start of gene = 1)
	if($intron == $ii-1 && $strand{$geneid} eq "+") {
	    $offset = $offset + $intronend - $intronstart + 1;
	}
	for($k=0;$k<@aa2;$k++) {
	    $geneWithIntron2indel_temp{$geneid}{$INTRON}[$indelcounter][0] = $offset + $aa2[$k][0];
	    $geneWithIntron2indel_temp{$geneid}{$INTRON}[$indelcounter][1] = $aa2[$k][1];  # the length of the indel
	    $TTT = $geneWithIntron2indel_temp{$geneid}{$INTRON}[$indelcounter][0];
	    $indelcounter++;
	}
	$aa[$ii] =~ /^[^:]+:(\d+)-(\d+)/;
	$exonstart = $1;
	$exonend = $2;
	$offset = $offset + $exonend - $exonstart + 1;
	if($intron == $ii-1 && $strand{$geneid} eq "-") {
	    $offset = $offset + $intronend - $intronstart + 1;
	}
	if($strand{$geneid} eq "-") {
	    $SEQ = $exonseq{$aa[$jj]} . $SEQ;
	}
	else {
	    $SEQ = $SEQ . $exonseq{$aa[$jj]};
	}
	if($intron == $jj) {
	    if($strand{$geneid} eq "-") {
		$SEQ = $intronseq{$INTRON} . $SEQ;
	    } else {
		$SEQ = $SEQ . $intronseq{$INTRON};
	    }
	}
    }
    $FLAG = 0;
    @aa = @{$geneWithIntron2indel_temp{$geneid}{$INTRON}};
    while($FLAG == 0) {
	$FLAG = 1;
	for($ii=0; $ii<@aa-1; $ii++) {
	    if($aa[$ii][0] > $aa[$ii+1][0]) {
		$temp = $aa[$ii][0];
		$aa[$ii][0] = $aa[$ii+1][0];
		$aa[$ii+1][0] = $temp;
		$temp = $aa[$ii][1];
		$aa[$ii][1] = $aa[$ii+1][1];
		$aa[$ii+1][1] = $temp;
		$FLAG = 0;
	    }
	}
    }
    for($ii=0; $ii<@aa; $ii++) {
	$geneWithIntron2indel{$geneid}{$INTRON}[$ii][0] = $aa[$ii][0];
	$geneWithIntron2indel{$geneid}{$INTRON}[$ii][1] = $aa[$ii][1];
	$TTT = $geneWithIntron2indel{$geneid}{$INTRON}[$ii][0];
    }
    undef %geneWithIntron2indel_temp;

    return $SEQ;
}

sub add_sequencing_error () {
    ($read, $error_rate) = @_;

    $Rlength = length($read);
    $num_errors = random_binomial(1, $Rlength, $error_rate);
    for($ie=0; $ie<$num_errors; $ie++) {
	$loc = int(rand($Rlength));
	$origbase = substr($read,$loc,1);
	$B = getrandombase();
	while($B eq $origbase) {
	    $B = getrandombase();
	}
	if(length($read) > $loc) {
	    $Z = substr($read,$loc,1,$B);
	}
    }
    return $read;
}

sub add_error_to_tail () {
    ($read, $low_qual_tail_length, $percent_of_tails_that_are_low_qual, $quality_of_tail) = @_;

    $Rnum = rand(1);
    if($Rnum < $percent_of_tails_that_are_low_qual) {
	$Rlength = length($read);
	@S4 = split(//,$read);
	for($ie=$Rlength-$low_qual_tail_length; $ie<@S4; $ie++) {
	    $Rnum2 = rand(1);
	    if($Rnum2 > $quality_of_tail) { # the lower $quality_of_tail is, the more likely we are to substitute
		$origbase = substr($read,$ie,1);
		$B = getrandombase();
		while($B eq $origbase) {
		    $B = getrandombase();
		}
		if(length($read) > $loc) {
		    $Z = substr($read,$ie,1,$B);
		}
	    }
	}
    }
    return $read;
}
