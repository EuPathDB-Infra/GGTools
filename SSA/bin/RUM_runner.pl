$date = `date`;

if(@ARGV == 1 && @ARGV[0] eq "config") {
    die "
The following describes the configuration file:
Note: All entries can be absolute path, or relative path to where the RUM_runner.pl script is.

1) gene annotation file, can be relative path to where the RUM_runner.pl script is, or absolute path
   e.g.: indexes/mm9_ucsc_refseq_gene_info.txt
2) bowtie executable, can be relative path to where the RUM_runner.pl script is, or absolute path
   e.g.: bowtie/bowtie
3) blat executable, can be relative path to where the RUM_runner.pl script is, or absolute path
   e.g.: blat/blat
4) mdust executable, can be relative path to where the RUM_runner.pl script is, or absolute path
   e.g.: mdust/mdust
5) bowtie genome index, can be relative path to where the RUM_runner.pl script is, or absolute path
   e.g.: indexes/mm9
6) bowtie gene index, can be relative path to where the RUM_runner.pl script is, or absolute path
   e.g.: indexes/mm9_genes_ucsc_refseq
7) blat genome index, can be relative path to where the RUM_runner.pl script is, or absolute path
   e.g. indexes/mm9_genome_sequence_single-line-seqs.fa
8) perl scripts directory, can be relative path to where the RUM_runner.pl script is, or absolute path
   e.g.: scripts

";
}
if(@ARGV < 5) {
    die "

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 _   _   _   _   _   _   
               // \\// \\// \\// \\// \\// \\/                       
              //\\_//\\_//\\_//\\_//\\_//\\_// 
        o_O__O_ o      
       | ====== |       .-----------.
       `--------'       |||||||||||||
        || ~~ ||        |-----------|
        || ~~ ||        | .-------. |
        ||----||        ! | UPENN | !
       //      \\\\        \\`-------'/  
      // /!  !\\ \\\\        \\_     _/
     !!__________!!         \\   /  
     ||  ~~~~~~  ||          `-'
     || _        ||
     |||_|| ||\\/|||
     ||| \\|_||  |||
     ||          ||
     ||  ~~~~~~  ||   
     ||__________||    
.----|||        |||------------------.
     ||\\\\      //||                 /|
     |============|                //
     `------------'               //
---------------------------------'/
---------------------------------'

        ______________________________________
     -  The RNA-Seq Unified Mapper (RUM) Usage  -
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage: RUM_runner.pl <configfile> <reads file(s)> <output dir> <num chunks>
                     <name> [options]

<reads file(s)> :  For unpaired data, the single file of reads.
                   For paired data either: 
                      - The files of forward and reverse reads, separated
                        by three commas ',,,'.
                      - One file formatted using parse2fasta.pl.
                   Files can be either fasta or fastq.

<num chunks>    :  The number of pieces to break the job into.  Use one
                   chunk unless you are on a cluster, or have multiple
                   processors with lots of RAM.

<name>          :  A string that will identify this run - use only alphanumeric
                   and underscores, no whitespace or other characters.

Options: -single    : Data is single-end (default is paired-end).

         -fast      : Run with blat params that run about 3 times faster but
                      a tad less sensitive

         -limitNUhard x : Limits the number of ambiguous mappers to a max of x

         -limitNU   : Limits the number of ambiguous mappers to a max of 100
                      locations.  If you have short reads and a large genome
                      then this will probably be necessary (45 bases is short
                      for mouse, 70 bases is long, between it's hard to say).

         -dna   : Run in dna mode, meaning don't map across splice junctions.

         -quantify : use this if using the -dna flag but you still want quantified
                     features.  If this is set you must have the gene models file
                     specified in the rum config file.

         -junctions : use this if using the -dna flag but you still want junction
                      calls.

         -minidentity x : run blat with minIdentity=x (default x=93)

         -countmismatches : report in the last column the number of mismatches,
                            ignoring insertions

         -variable_read_lengths : set this if your reads are not all of the 
                                  same length

         -qsub      : Use qsub to fire the job off to multiple nodes.  This
                      means you're on a cluster that understands qsub.  If not
                      using -qsub, you can still break it into multiple chunks,
                      it will just fire each chunk off to a separate processor.
                      Don't use more chunks than you have processors though,
                      because that will just slow it down.

         -max_inertionss_per_read n : Allow at most n insertions in one read.  
                      The default is n=1.  Setting n>1 is only allowed for single
                      end reads.  Don't raise it unless you know what you are
                      doing, because it can greatly increase the false alignments.

         -kill      : To kill a job, run with all the same parameters but add
                      -kill.  Note: it is not sufficient to just terminate
                      RUM_runner.pl, that will leave other phantom processes.
                      Use -kill instead.

Running RUM_runner.pl with the one argument 'config' will explain how to make
the config file.

This program writes very large intermediate files.  If you have a large genome
such as mouse or human then it is recommended to run in chunks on a cluster, or
a machine with multiple processors.  Running with under five million reads per
chunk is usually best, and getting it under a million reads per chunk will speed
things considerably.

You can put an 's' after the number of chunks if they have already been broken
into chunks, so as to avoid repeating this time-consuming step.

Usage (again): RUM_runner.pl <configfile> <reads file(s)> <output dir> <num chunks>
                     <name> [options]

";
}

$configfile = $ARGV[0];
$readsfile = $ARGV[1];
$output_dir = $ARGV[2];
$output_dir =~ s!/$!!;
$numchunks = $ARGV[3];
$name = $ARGV[4];
if($name =~ /^-/) {
    die "\nError: The name '$name' is invalid, probably you forgot a required argument\n\n";
}

$name_o = $ARGV[4];
$name =~ s/\s+/_/g;
$name =~ s/[^a-zA-Z0-9_]//g;

if($name ne $name_o) {
    print STDERR "\nWarning: name changed from '$name_o' to '$name'.\n\n";
}
$paired_end = "true";
$fast = "false";
$dna = "false";
$limitNU = "false";
$limitNUhard = "false";
$qsub = "false";
$minidentity=93;
$kill = "false";
$variable_read_lengths = "false";
$countmismatches = "false";
$num_insertions_allowed = 1;
$junctions = "false";
$quatify = "false";
if(@ARGV > 5) {
    for($i=5; $i<@ARGV; $i++) {
	$optionrecognized = 0;
        if($ARGV[$i] eq "-max_ins_per_read") {
	    $i++;
	    $num_insertions_allowed = $ARGV[$i];
            if($ARGV[$i] =~ /^\d+$/) {
	        $optionrecognized = 1;
	    }
        }
	if($ARGV[$i] eq "-single") {
	    $paired_end = "false";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-junctions") {
	    $junctions = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-quantify") {
	    $quantify = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-countmismatches") {
	    $countmismatches = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-variable_read_lengths") {
	    $variable_read_lengths = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-fast") {
	    $fast = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-kill") {
	    $kill = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-dna") {
	    $dna = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-limitNU") {
	    $limitNU = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-qsub") {
	    $qsub = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-minidentity") {

	    $minidentity = $ARGV[$i+1];
	    $i++;
	    if(!($minidentity =~ /^\d+$/ && $minidentity <= 100)) {
		die "\nERROR: minidentity must be an integer between zero and 100.\nYou have given '$minidentity'.\n\n";
	    }
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-limitNUhard") {
	    $NU_limit = $ARGV[$i+1];
	    $i++;
	    $limitNUhard = "true";
	    if(!($NU_limit =~ /^\d+$/ && $NU_limit > 0)) {
		die "\nERROR: -limitNUhard must be an integer greater than zero.\nYou have given '$NU_limit'.\n\n";
	    }
	    $optionrecognized = 1;
	}
	if($optionrecognized == 0) {
	    print STDERR "\nERROR: option $ARGV[$i] not recognized.\n\n";
	    exit();
	}
    }
}
if($dna eq "false") {
    $junctions = "true";
    $quantify = "true";
}

if($kill eq "true") {
    $outdir = $output_dir;
    $str = `ps x | grep $outdir`;
    @candidates = split(/\n/,$str);
    for($i=0; $i<@candidates; $i++) {
	if($candidates[$i] =~ /^\s*(\d+)\s.*(\s|\/)$outdir\/pipeline.\d+.sh/) {
	    $pid = $1;
	    print STDERR "killing $pid\n";
	    `kill -9 $pid`;
	}
    }
    for($i=0; $i<@candidates; $i++) {
	if($candidates[$i] =~ /^\s*(\d+)\s.*(\s|\/)$outdir(\s|\/)/) {
	    if(!($candidates[$i] =~ /pipeline.\d+.sh/)) {
		$pid = $1;
		print STDERR "killing $pid\n";
		`kill -9 $pid`;
	    }
	}
    }
    exit();
}

print STDERR "

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 _   _   _   _   _   _   
               // \\// \\// \\// \\// \\// \\/                       
              //\\_//\\_//\\_//\\_//\\_//\\_// 
        o_O__O_ o      
       | ====== |       .-----------.
       `--------'       |||||||||||||
        || ~~ ||        |-----------|
        || ~~ ||        | .-------. |
        ||----||        ! | UPENN | !
       //      \\\\        \\`-------'/  
      // /!  !\\ \\\\        \\_     _/
     !!__________!!         \\   /  
     ||  ~~~~~~  ||          `-'
     || _        ||
     |||_|| ||\\/|||
     ||| \\|_||  |||
     ||          ||
     ||  ~~~~~~  ||   
     ||__________||    
.----|||        |||------------------.
     ||\\\\      //||                 /|
     |============|                //
     `------------'               //
---------------------------------'/
---------------------------------'
  ____________________________________________________________
- The RNA-Seq Unified Mapper (RUM) Pipeline has been initiated -
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
";

sleep(2);
print STDERR "Please wait while I check that everything is in order.\n\n";
sleep(2);
print STDERR "This could take a few minutes.\n\n";
sleep(2);

$check = `ps x | grep RUM_runner.pl`;
@a = split(/\n/,$check);
$CNT=0;
for($i=0; $i<@a; $i++) {
    chomp($a[$i]);
    $a[$i] =~ s/.*RUM_runner.pl *//;
    @b = split(/ +/,$a[$i]);
    if($b[2] eq $output_dir || $b[2] eq $ARGV[2]) {
	$CNT++;
	if($CNT > 1) {
	    die "\nERROR: You seem to already have an instance of RUM_runner.pl running on the\nsame working directory.  This will cause collisions of the temporary files.\n\nExiting...\n\n";
	}
    }
}

for($i=1; $i<=$numchunks; $i++) {
    $logfile = "$output_dir/rum_log.$i";
    if (-e $logfile) {
	unlink($logfile);
    }
}

open(INFILE, $configfile);
$gene_annot_file = <INFILE>;
chomp($gene_annot_file);
$bowtie_exe = <INFILE>;
chomp($bowtie_exe);
$blat_exe = <INFILE>;
chomp($blat_exe);
$mdust_exe = <INFILE>;
chomp($mdust_exe);
$genome_bowtie = <INFILE>;
chomp($genome_bowtie);
$transcriptome_bowtie = <INFILE>;
chomp($transcriptome_bowtie);
$genome_blat = <INFILE>;
chomp($genome_blat);
$scripts_dir = <INFILE>;
$scripts_dir =~ s!/$!!;
chomp($scripts_dir);
$genomefa = $genome_blat;
close(INFILE);

if($num_insertions_allowed > 1 && $paired_end eq "true") {
    die "\nError: for paired end data, you cannot set -max_ins_per_read to be greater than 1.\n\n";
}

if(($readsfile =~ /,,,/) && $paired_end eq "false") {
    die "\nError: You've given two files separated by three commas, this means we are expecting paired end\ndata but you set -single.\n\n";
}
if(!($readsfile =~ /,,,/) && !(-e $readsfile)) {
    die "\nError: The reads file '$readsfile' does not seem to exist\n\n";
}
$quals = "false";
if(($readsfile =~ /,,,/) && ($paired_end eq "true")) {
    @a = split(/,,,/, $readsfile);
    if(@a > 2) {
	die "\nError: You've given more than two files separated by three commas, should be at most two files.\n\n";
    }    
    if(!(-e $a[0])) {
	die "\nError: The reads file '$a[0]' does not seem to exist\n\n";
    }
    if(!(-e $a[1])) {
	die "\nError: The reads file '$a[1]' does not seem to exist\n\n";
    }
    if($a[0] eq $a[1]) {
	die "\nError: You specified the same file for the forward and reverse reads, must be an error...\n\n";
    }

    print STDERR "Reformatting reads file...\n";
    `perl $scripts_dir/parse2fasta.pl $a[0] $a[1] > $output_dir/reads.fa`;
    `perl $scripts_dir/fastq2qualities.pl $a[0] $a[1] > $output_dir/quals.fa`;
    $X = `head -2 $output_dir/quals.fa | tail -1`;
    if($X =~ /\S/ && !($X =~ /Sorry, can't figure these files out/)) {
	$quals = "true"
    }
    $readsfile = "$output_dir/reads.fa";
    $qualsfile = "$output_dir/quals.fa";
}

open(LOGFILE, ">$output_dir/rum.log");
print LOGFILE "config file: $configfile\n";
print LOGFILE "readsfile: $readsfile\n";
print LOGFILE "output_dir: $output_dir\n";
print LOGFILE "numchunks: $numchunks\n";
print LOGFILE "name: $name\n";
print LOGFILE "paired_end: $paired_end\n";
print LOGFILE "fast: $fast\n";
print LOGFILE "limitNU: $limitNU\n";
print LOGFILE "dna: $dna\n";
print LOGFILE "qsub: $qsub\n";
print LOGFILE "blat minidentity: $minidentity\n";
print LOGFILE "output junctions: $junctions\n";
print LOGFILE "output quantified values: $quantify\n";

if($numchunks =~ /(\d+)s/) {
    $numchunks = $1;
    $fasta_already_fragmented = "true";
} else {
    $fasta_already_fragmented = "false";
}

$head = `head -2 $readsfile | tail -1`;
chomp($head);
@a = split(//,$head);
if($variable_read_lengths eq "false") {
   $readlength = @a;
} else {
   $readlength = "v";
}
print LOGFILE "readlength: $readlength\n";
print LOGFILE "\nstart: $date\n";

$head = `head -4 $readsfile`;
$head =~ /seq.(\d+)(.).*seq.(\d+)(.)/s;
$num1 = $1;
$type1 = $2;
$num2 = $3;
$type2 = $4;
if($paired_end eq 'false') {
    if($type1 ne "a") {
	print STDERR "Reformatting reads file...\n";
	`perl scripts/parse2fasta.pl $readsfile > $output_dir/reads.fa`;
	`perl scripts/fastq2qualities.pl $readsfile > $output_dir/quals.fa`;
	$X = `head -2 $output_dir/quals.fa | tail -1`;
	if(!($X =~ /Sorry, can't figure/)) {
	    $quals = "true"
	}
	$readsfile = "$output_dir/reads.fa";
	$qualsfile = "$output_dir/quals.fa";
	$head = `head -4 $readsfile`;
	$head =~ /seq.(\d+)(.).*seq.(\d+)(.)/s;
	$num1 = $1;
	$type1 = $2;
	$num2 = $3;
	$type2 = $4;
    }
}

if($type1 ne "a") {
    print STDERR "ERROR: the fasta def lines are misformatted.  The first one should end in an 'a'.\n";
    print LOGFILE "Error: fasta file misformatted... The first line should end in an 'a'.\n";
    exit();
}
if($num2 ne "2" && $paired_end eq "false") {
    print STDERR "ERROR: the fasta def lines are misformatted.  The second one should be '2a' or '1b'\n";
    print STDERR "       depending on whether it is paired end or not.  ";
    print LOGFILE "Error: fasta file misformatted...  The second line should be '2a' or '1b' depending\n";
    print LOGFILE "on whether it is paired end or not..\n";
    if($paired_end eq "true") {
	print STDERR "Note: You are running in paired end mode.\n";
	print LOGFILE "Note: You ran in paired end mode.\n";
    }
    else {
	print STDERR "Note: You are not running in paired end mode.\n";
	print LOGFILE "Note: You ran in unpaired mode.\n";
    }
    exit();
}
if($type2 ne "b" && $paired_end eq "true") {
    print STDERR "ERROR: the fasta def lines are misformatted.  You are in paired end mode so the second\n";
    print STDERR "       one should end in a 'b'.\n";
    print LOGFILE "Error: fasta file misformatted...  You ran in paried end mode so the second\n";
    print LOGFILE "one should end in a 'b'.\n";
    exit();
}
if($type1 eq "a" && $type2 eq "a" && $paired_end eq "true") {
    print STDERR "ERROR: You are running in paired end mode, paired reads should have def\n";
    print STDERR "       lines '>seq.Na' and '>seq.Nb' for N = 1, 2, 3, ... ";
    print LOGFILE "Error: fasta file misformatted...\n";
    print LOGFILE "You are running in paired end mode, paired reads should have def\n";
    print LOGFILE "lines '>seq.Na' and '>seq.Nb' for N = 1, 2, 3, ... ";
    exit();
}
if($paired_end eq "false" && $type1 eq "a" && $type2 eq "a" && ($num1 ne "1" || $num2 ne "2")) {
    print STDERR "ERROR: You ran in unpaired mode, reads should have def\n";
    print STDERR "lines '>seq.Nb' for N = 1, 2, 3, ... ";
    print LOGFILE "Error: fasta file misformatted...\n";
    print LOGFILE "You ran in unpaired mode, reads should have def\n";
    print LOGFILE "lines '>seq.Nb' for N = 1, 2, 3, ... ";
    exit();
}
if($paired_end eq "true" && ($type1 ne "a" || $type2 ne "b")) {
    print STDERR "ERROR: You ran in paired mode, reads should have def\n";
    print STDERR "lines alternating '>seq.Na' and '>seq.Nb' for N = 1, 2, 3, ... ";
    print LOGFILE "Error: fasta file misformatted...\n";
    print LOGFILE "You ran in paired mode, reads should have def\n";
    print LOGFILE "lines alternating '>seq.Na' and '>seq.Nb' for N = 1, 2, 3, ... ";
    exit();
}

if($readlength ne "variable" && $readlength < 55 && $limitNU eq "false") {
    print STDERR "\n\nWARNING: you have pretty short reads ($readlength bases).  If you have a large\n";
    print STDERR "genome such as mouse or human then the files of ambiguous mappers could grow\n";
    print STDERR "very large, in this case it's recommended to run with the -limitNU option.  You\n";
    print STDERR "can watch the files that start with 'X' and 'Y' to see if they are growing\n";
    print STDERR "larger than 10 gigabytes per million reads at which point you might want to use.\n";
    print STDERR "-limitNU\n\n";
}

$pipeline_template = `cat pipeline_template.sh`;
if($dna eq "true") {
    $pipeline_template =~ s/# cp /cp /gs;
    $pipeline_template =~ s/xxx1.*xxx2//s;
}
if($fasta_already_fragmented eq "false") {
    print STDERR "Splitting files ...\n\n";
    $qualflag = 0;
    $x = breakup_file($readsfile, $numchunks);
    if($quals eq "true") {
        $qualflag = 1;
	$x = breakup_file($qualsfile, $numchunks);
    }
}

print STDERR "Reads fasta file already fragmented: $fasta_already_fragmented\n";
print STDERR "Number of Chunks: $numchunks\n";
print STDERR "Reads File: $readsfile\n";
print STDERR "Paired End: $paired_end\n";

print STDERR "\nEverything seems okay, I am going to fire off the job.\n\n";

for($i=1; $i<=$numchunks; $i++) {
    $pipeline_file = $pipeline_template;
    if($limitNUhard eq "true") {
	$pipeline_file =~ s!LIMITNUCUTOFF!$NU_limit!gs;
	$pipeline_file =~ s!sort_RUM.pl OUTDIR.RUM_NU_temp2.CHUNK!sort_RUM.pl OUTDIR/RUM_NU_temp3.CHUNK!gs;
    } else {
	$pipeline_file =~ s!perl SCRIPTSDIR/limit_NU.pl OUTDIR/RUM_NU_temp2.CHUNK LIMITNUCUTOFF > OUTDIR/RUM_NU_temp3.CHUNK\n!!gs;
    }
    if($num_insertions_allowed != 1) {
	$pipeline_file =~ s!MAXINSERTIONSALLOWED!-num_insertions_allowed $num_insertions_alllowed!gs;
    } else {
	$pipeline_file =~ s!MAXINSERTIONSALLOWED!!gs;
    }
    $pipeline_file =~ s!OUTDIR!$output_dir!gs;
    if($quals eq "false") {
	$pipeline_file =~ s!QUALSFILE.CHUNK!none!gs;
    } else {
	$pipeline_file =~ s!QUALSFILE!$qualsfile!gs;
    }
    $pipeline_file =~ s!CHUNK!$i!gs;
    $pipeline_file =~ s!MINIDENTITY!$minidentity!gs;
    $pipeline_file =~ s!BOWTIEEXE!$bowtie_exe!gs;
    $pipeline_file =~ s!GENOMEBOWTIE!$genome_bowtie!gs;
    $pipeline_file =~ s!BOWTIEEXE!$bowtie_exe!gs;
    $pipeline_file =~ s!READSFILE!$readsfile!gs;
    $pipeline_file =~ s!SCRIPTSDIR!$scripts_dir!gs;
    $pipeline_file =~ s!TRANSCRIPTOMEBOWTIE!$transcriptome_bowtie!gs;
    $pipeline_file =~ s!GENEANNOTFILE!$gene_annot_file!gs;
    $pipeline_file =~ s!BLATEXE!$blat_exe!gs;
    $pipeline_file =~ s!MDUSTEXE!$mdust_exe!gs;
    $pipeline_file =~ s!GENOMEBLAT!$genome_blat!gs;
    $pipeline_file =~ s!GENOMEFA!$genomefa!gs;
    $pipeline_file =~ s!READLENGTH!$readlength!gs;
    if($countmismatches eq "true") {
	$pipeline_file =~ s!COUNTMISMATCHES!-countmismatches!gs;
    } else {
	$pipeline_file =~ s!COUNTMISMATCHES!!gs;
    }
    if($limitNU eq "true") {
	$pipeline_file =~ s! -a ! -k 100 !gs;	
    }
    if($fast eq "false") {
	$pipeline_file =~ s!SPEED!-stepSize=5!gs;
    }
    else {
	$pipeline_file =~ s!SPEED!!gs;
    }
    if($paired_end eq "true") {
	$pipeline_file =~ s!PAIREDEND!paired!gs;
    } else {
	$pipeline_file =~ s!PAIREDEND!single!gs;
    }
    $outfile = "pipeline." . $i . ".sh";
    open(OUTFILE, ">$output_dir/$outfile") or die "\nError: cannot open '$output_dir/$outfile' for writing\n\n";
    print OUTFILE $pipeline_file;
    close(OUTFILE);

    if($qsub eq "true") {
	`qsub -l mem_free=6G -pe DJ 4 $output_dir/$outfile`;
    }
    else {
	system("/bin/bash $output_dir/$outfile &");
    }
    print STDERR "Chunk $i initiated\n";
    $status{$i} = 1;
}
if($numchunks > 1) {
    print STDERR "\nAll chunks initiated, now the long wait...\n";
    print STDERR "\nI'm going to watch for all chunks to finish, then I will merge everything...\n";
    sleep(2);
    if($qsub eq "false") {
	print STDERR "\nThe next thing to print here will be the status reports from bowtie.\n";
	print STDERR "     * don't be alarmed.\n\n";
    }
} else {
    print STDERR "\nThe job has been initiated, now the long wait...\n";
    sleep(2);
    if($qsub eq "false") {
	print STDERR "\nThe next thing to print here will be the status reports from bowtie.\n";
	print STDERR "     * don't be alarmed.\n\n";
    }
}

$currenttime = time();
$lastannouncetime = $currenttime;
$numannouncements = 0;
$doneflag = 0;

while($doneflag == 0) {
    $doneflag = 1;
    $numdone = 0;
    for($i=1; $i<=$numchunks; $i++) {
	$logfile = "$output_dir/rum_log.$i";
	if (-e $logfile) {
	    $x = `cat $logfile`;
	    if(!($x =~ /pipeline complete/s)) {
		$doneflag = 0;
	    } else {
		$numdone++;
		if($status{$i} == 1) {
		    $status{$i} = 2;
		    print STDERR "\n *** Chunk $i has finished.\n";
		}
	    }
	}
	else {
	    $doneflag = 0;
	}
    }
    if($doneflag == 0) {
	sleep(30);
	$currenttime = time();
	if($currenttime - $lastannouncetime > 3600) {
	    $lastannouncetime = $currenttime;
	    $numannouncements++;
	    if($numannouncements == 1) {
		if($numdone == 1) {
		    print STDERR "\nIt has been $numannouncements hour, $numdone chunk has finished.\n";
		} else {
		    print STDERR "\nIt has been $numannouncements hour, $numdone chunks have finished.\n";
		}
	    } else {
		if($numdone == 1) {
		    print STDERR "\nIt has been $numannouncements hours, $numdone chunk has finished.\n";
		} else {
		    print STDERR "\nIt has been $numannouncements hours, $numdone chunks have finished.\n";
		}
	    }
	}
    }
}

print STDERR "All chunks have finished, now to move on to the merging,\n";
print STDERR "creating the coverage plot and the quantified values, etc..\n\n";

$date = `date`;
print LOGFILE "finished creating RUM_Unique.*/RUM_NU.*: $date\n";
$x = `cp $output_dir/RUM_Unique.1 $output_dir/RUM_Unique`;
for($i=2; $i<=$numchunks; $i++) {
    $x = `cat $output_dir/RUM_Unique.$i >> $output_dir/RUM_Unique`;
}
$x = `cp $output_dir/RUM_NU.1 $output_dir/RUM_NU`;
for($i=2; $i<=$numchunks; $i++) {
    $x = `cat $output_dir/RUM_NU.$i >> $output_dir/RUM_NU`;
}
$x = `cp $output_dir/RUM.sam.1 $output_dir/RUM.sam`;
for($i=2; $i<=$numchunks; $i++) {
    $x = `cat $output_dir/RUM.sam.$i >> $output_dir/RUM.sam`;
}
print LOGFILE "finished creating RUM_Unique/RUM_NU/RUM.sam: $date\n";
print LOGFILE "starting M2C: $date\n";
$PPlog = "postprocessing_$name" . ".log";
$shellscript = "#!/bin/sh\n";
$shellscript = $shellscript . "perl $scripts_dir/count_reads_mapped.pl $output_dir/RUM_Unique $output_dir/RUM_NU > $output_dir/mapping_stats.txt\n";
$shellscript = $shellscript . "echo making bed > $output_dir/$PPlog\n";
$shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
$shellscript = $shellscript . "perl $scripts_dir/make_bed.pl $output_dir/RUM_Unique $output_dir/RUM_Unique.bed -zbho\n";
$shellscript = $shellscript . "echo starting M2C >> $output_dir/$PPlog\n";
$shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
$covfilename = $name . ".cov";
$logfilename = $name . ".log";
$shellscript = $shellscript . "java -Xmx2000m M2C $output_dir/RUM_Unique.bed $output_dir/RUM_$covfilename $output_dir/RUM_$logfilename -ucsc -name \"$name\" -start_coordinate_infile 0 -openinterval_infile -chunks 4\n";
if($quantify eq "true") {
   $shellscript = $shellscript . "echo starting to quantify features >> $output_dir/$PPlog\n";
   $shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
   $shellscript = $shellscript . "perl $scripts_dir/quantify_one_sample.pl $output_dir/RUM_$name";
   $shellscript = $shellscript . ".cov $gene_annot_file -zero -open > $output_dir/feature_quantifications_$name\n";
}
if($junctions eq "true") {
   $shellscript = $shellscript . "echo starting to compute junctions >> $output_dir/$PPlog\n";
   $shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
   $shellscript = $shellscript . "perl $scripts_dir/make_RUM_junctions_file.pl $output_dir/RUM_Unique $output_dir/RUM_NU $genomefa $gene_annot_file $output_dir/junctions_all.rum $output_dir/junctions_all.bed $output_dir/junctions_high-quality.bed -faok\n";
}
$shellscript = $shellscript . "echo finished >> $output_dir/$PPlog\n";
$shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
$str = "postprocessing_$name" . ".sh";
open(OUTFILE2, ">$output_dir/$str");
print OUTFILE2 $shellscript;
close(OUTFILE2);

if($qsub eq "true") {
    `qsub -l mem_free=6G -pe DJ 4 $output_dir/$str`;
}
else {
    system("/bin/bash $output_dir/$str");
}

print STDERR "\nWorking, now another wait...\n\n";
$doneflag = 0;
while($doneflag == 0) {
    $doneflag = 1;
    if (-e "$output_dir/$PPlog") {
	$x = `cat $output_dir/$PPlog`;
	if(!($x =~ /finished/s)) {
	    $doneflag = 0;
	}
    }
    else {
	$doneflag = 0;
    }
    if($doneflag == 0) {
	sleep(30);
    }
}

print STDERR "\nOkay, all finished.\n\n";

$date = `date`;
print LOGFILE "pipeline finished: $date\n";
close(LOGFILE);

sub breakup_file () {
    ($FILE, $numpieces) = @_;

    open(INFILE, $FILE);
    $filesize = `wc -l $FILE`;
    chomp($filesize);
    $filesize =~ s/^\s+//;
    $filesize =~ s/\s.*//;
    $numseqs = $filesize / 2;
    $piecesize = int($numseqs / $numpieces);
    $piecesize2 = int($numseqs / $numpieces / 2);
    if($numchunks > 1) {
	print LOGFILE "processing in $numchunks pieces of approx $piecesize2 reads each\n";
    } else {
	print LOGFILE "processing in one piece of approx $piecesize2 reads\n";
    }
    if($piecesize % 2 == 1) {
	$piecesize++;
    }
    $bflag = 0;

    for($i=1; $i<$numpieces; $i++) {
	$F2 = $FILE;
	$F2 =~ s!.*/!!;
	$outfilename = $output_dir . "/" . $F2 . "." . $i;

	open(OUTFILE, ">$outfilename");
	for($j=0; $j<$piecesize; $j++) {
	    $line = <INFILE>;
	    chomp($line);
	    $line =~ s/\^M$//s;
	    if($qualflag == 0) {
		$line =~ s/[^ACGTNab]$//s;
	    }
	    print OUTFILE "$line\n";
	    $line = <INFILE>;
	    chomp($line);
	    $line =~ s/\^M$//s;
	    if($qualflag == 0) {
		$line =~ s/[^ACGTNab]$//s;
	    }
	    print OUTFILE "$line\n";
	}
	close(OUTFILE);
    }
    $outfilename = $FILE . "." . $numpieces;
    open(OUTFILE, ">$outfilename");
    while($line = <INFILE>) {
	print OUTFILE $line;
    }
    close(OUTFILE);
    return 0;
}

sub checkstatus () {
    ($CHUNK) = @_;
    $log = `cat $output_dir/rum_log.$CHUNK`;
    @LOG = split(/\n/, $log);
    $LOG[1] =~ /(\d+)$/;
    $started_at = $1;
    if($log =~ /finished first bowtie run/s) {
	$LOG[3] =~ /(\d+)$/;
	$firstbowtie_at = $1;
	$firstbowtie_filesize = -s "$output_dir/X.$CHUNK";
	if($firstbowtie_filesize == 0) {
	    print STDERR "Warning: genome bowtie outfile for chunk $CHUNK is empty.\n";
	}
    } else {
	$firstbowtie_lastmodified =(stat ($output_dir/X.$CHUNK))[9];
	if($firstbowtie_lastmodified - $started_at > 600) {

	}
    }
    if($log =~ /finished parsing genome bowtie run/s) {
	$LOG[3] =~ /(\d+)$/;
	$firstbowtie_at = $1;
	$firstbowtie_filesize = -s "$output_dir/X.$CHUNK";
	if($firstbowtie_filesize == 0) {
	    print STDERR "Warning: genome bowtie outfile for chunk $CHUNK is empty.\n";
	}
    }

}
