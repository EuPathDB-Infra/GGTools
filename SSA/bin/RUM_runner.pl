
# Written by Gregory R Grant
# University of Pennsylvania, 2010

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
9) lib directory, this directory holds the config files and the pipeline_template.sh file.  It can
   be relative path to where the RUM_runner.pl script is, or absolute path
   e.g.: lib

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
      // /!  !\\ \\\\        \\_  O  _/
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
  .------------------------------------.
  | RNA-Seq Unified Mapper (RUM) Usage |
  ` ================================== '

Usage: RUM_runner.pl <config file> <reads file(s)> <output dir> <num chunks>
                     <name> [options]

<config file>   :  This file tells RUM where to find the various executables
                   and indexes.  This file is included in the 'lib' directory
                   when you download an organism, for example rum.config_mm9
                   for mouse build mm9, which will work if you leave everything
                   in its default location.  To modify or make your own config
                   file, run this program with the single argument 'config' for
                   more information on the config file.

<reads file(s)> :  What to put here depends on whether your data is paired or
                   unpaired:

                   1) For unpaired data, the single file of reads.
                      - Don't forget to set the -single option in this case
                   2) For paired data, either: 
                      - The files of forward and reverse reads, separated
                        by three commas ',,,' (with no spaces).
                      - One file formatted using parse2fasta.pl.

                   NOTE ON FILE FORMATS: Files can be either fasta or fastq,
                   the type is inferred.

<output dir>    :  Where to write the temp, interemediate, and results files.

<num chunks>    :  The number of pieces to break the job into.  Use one chunk
                   unless you are on a cluster, or have multiple cores
                   with lots of RAM.  Have at least one processing core per
                   chunk.  A genome like human will also need about 5 to 6 Gb
                   of RAM per chunk.  Even with a small genome, if you have
                   tens of millions of reads, you will still need a few Gb of
                   RAM to get through the post-processing.

<name>          :  A string to identify this run - use only alphanumeric,
                   underscores, and dashes.  No whitespace or other characters.

Options: There are many options, but RUM is typically run with the defaults.
         Just don't forget to set -single if you have single-end data.  The
         option -kill is also quite useful to stop a run, because killing just
         the main program will not always kill the spawned processes.

       -single    : Data is single-end (default is paired-end).

       -strandspecific : If the data are strand specific, then yo can use this
                         option to generate strand specific coverage plots and
                         quantified values.

       -dna       : Run in dna mode, meaning don't map across splice junctions.

       -fast      : Run with blat params that run about 3 times faster but
                    a tad less sensitive.

       -variable_read_lengths : Set this if your reads are not all of the same
                                length.

       -limitNU N : Limits the number of ambiguous mappers in the final output
                    by removing all reads that map to N locations or more.

       -limitBowtieNU : Limits the number of ambiguous mappers in the Bowtie
                        run to a max of 100.  If you have short reads and a
                        large genome, or a very repetitive genome, this might
                        be necessary to keep the bowtie files from getting out
                        of hand - 10Gb per lane is not abnormal but 100Gb might
                        be. (note: 45 bases is considered short for mouse, 70
                        bases considered long, between it's hard to say).

       -quantify : Use this *if* using the -dna flag and you still want quantified
                   features.  If this is set you must have the gene models file
                   specified in the rum config file.  Without the -dna flag
                   quantified features are generated by default so you don't
                   need to set this. 

       -junctions : Use this *if* using the -dna flag and you still want junction
                    calls.  If this is set you should have the gene models file
                    specified in the rum config file (if you have one).  Without
                    the -dna flag junctions generated by default so you don't
                    need to set this. 

       -ram n : Specify the max number of Gb of ram you want to dedicate to each
                chunk, if less than seven.  If you have at least 6 or 7 per chunk
                then don't worry about this.  The program will try to figure out
                the amount of available RAM automatically, so this option is
                mainly just a backup for when RAM is limited but the program is
                unable to figure that out - in this case a warning message will
                be generated.

       -minidentity x : run blat with minIdentity=x (default x=93).  You
                        shouldn't need to change this.

       -minlength x : don't report alignments less than this long.  The default
                      = 50 if the readlength >= 80, else = 35 if readlength >= 45
                      else = 0.8 * readlength.  Don't set this too low you will
                      start to pick up a lot of garbage.

       -countmismatches : Report in the last column the number of mismatches,
                          ignoring insertions

       -altgenes x : x is a file with gene models to use for calling junctions
                     novel.  If not specified will use the gene models file
                     specified in the config file. 

       -qsub      : Use qsub to fire the job off to multiple nodes on a cluster.
                    This means you're on a cluster that understands qsub.  

                      ** Note: without using -qsub, you can still specify more
                      than one chunk.  It should fire each chunk off to a
                      separate core.  But don't use more chunks than you have
                      cores, because that can slow things down considerable.

       -max_insertions_per_read n : Allow at most n insertions in one read.  
                    The default is n=1.  Setting n>1 is only allowed for single
                    end reads.  Don't raise it unless you know what you are
                    doing, because it can greatly increase the false alignments.

       -postprocess : Rerun just the post-processing steps, after the alignment,
                      if for some reason you need to do this.

       -noclean   : do not remove the intermediate and temp files after finishing.

       -kill      : To kill a job, run with all the same parameters but add
                    -kill.  Note: it is not sufficient to just terminate
                    RUM_runner.pl, that will leave other phantom processes.
                    Use -kill instead.

Default config files are supplied with each organism.  If you need to make or
modify one then running RUM_runner.pl with the one argument 'config' gives an
explaination of the the file.

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
if(!(-d $output_dir)) {
    die "\nError: The directory '$output_dir' does not seem to exists...\n\n";
}
$numchunks = $ARGV[3];
$name = $ARGV[4];
if($name =~ /^-/) {
    die "\nError: The name '$name' is invalid, probably you forgot a required argument\n\n";
}

$name_o = $ARGV[4];
$name =~ s/\s+/_/g;
$name =~ s/[^a-zA-Z0-9_-]//g;

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
$minlength=0;
$postprocess = "false";
$blatonly = "false";
$kill = "false";
$cleanup = "true";
$variable_read_lengths = "false";
$countmismatches = "false";
$num_insertions_allowed = 1;
$junctions = "false";
$quatify = "false";
$ram = 8;
$user_ram = "false";
$nocat = "false";
$quals_specified = "false";
$strandspecific = "false";
if(@ARGV > 5) {
    for($i=5; $i<@ARGV; $i++) {
	$optionrecognized = 0;

        if($ARGV[$i] eq "-max_insertions_per_read") {
	    $i++;
	    $num_insertions_allowed = $ARGV[$i];
            if($ARGV[$i] =~ /^\d+$/) {
	        $optionrecognized = 1;
	    }
        }
        if($ARGV[$i] eq "-ram") {
	    $i++;
	    $ram = $ARGV[$i];
            $user_ram = "true";
            if($ARGV[$i] =~ /^\d+$/) {
	        $optionrecognized = 1;
	    }
        }
	if($ARGV[$i] eq "-single") {
	    $paired_end = "false";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-strandspecific") {
	    $strandspecific = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-nocat") {
	    $nocat = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-noclean") {
	    $cleanup = "false";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-junctions") {
	    $junctions = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-postprocess") {
	    $postprocess = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-blatonly") {
	    $blatonly = "true";
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
	if($ARGV[$i] eq "-limitBowtieNU") {
	    $limitNU = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-qsub") {
	    $qsub = "true";
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-altgenes") {
	    $altgenes = "true";
            $i++;
            $altgene_file = $ARGV[$i];
            open(TESTIN, $altgene_file) or die "\nError: cannot open '$altgene_file' for reading.\n\n";
            close(TESTIN);
	    $optionrecognized = 1;
	}

	if($ARGV[$i] eq "-qualsfile" || $ARGV[$i] eq "-qualfile") {
	    $quals_specified = "true";
            $i++;
            $quals_file = $ARGV[$i];
            if($quals_file =~ /\//) {
               die "Error: do not specify -quals file with a full path, put it in the '$output_dir' directory.\n\n";
            }
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
	if($ARGV[$i] eq "-minlength") {
	    $minlength = $ARGV[$i+1];
	    $i++;
	    if(!($minlength =~ /^\d+$/ && $minlength >= 10)) {
		die "\nERROR: minlength must be an integer >= 10.\nYou have given '$minlength'.\n\n";
	    }
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-limitNU") {
	    $NU_limit = $ARGV[$i+1];
	    $i++;
	    $limitNUhard = "true";
	    if(!($NU_limit =~ /^\d+$/ && $NU_limit > 0)) {
		die "\nERROR: -limitNU must be an integer greater than zero.\nYou have given '$NU_limit'.\n\n";
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
      // /!  !\\ \\\\        \\_  O  _/
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

if($qsub eq "true") {
     print STDERR "\nWarning: You are using qsub - so if you have installed RUM somewhere other than your\nhome directory, then you will probably need to specify everything with full paths\nor this may not work.\n\n";
}

if($postprocess eq "false") {
     sleep(2);
     print STDERR "Please wait while I check that everything is in order.\n\n";
     sleep(2);
     print STDERR "This could take a few minutes.\n\n";
     sleep(2);

     if($qsub eq "true") {
         print STDERR "You have chosen to submit the jobs using 'qsub'.  I'm going to assume each node has\nsufficient RAM for this.  If you are running a mammalian genome then you should have\nat least 6 to 8 Gigs per node.\n\n";
     } else {
          print STDERR "I'm going to try to figure out how much RAM you have.\nIf you see some error messages here, don't worry, these are harmless.\n\n";
          sleep(2);
          # figure out how much RAM is available:
          if($user_ram eq "false") {
               $did_not_figure_out_ram = "false";
               $ramcheck = `free -g`;  # this should work on linux
               $ramcheck =~ /Mem:\s+(\d+)/s;
               $totalram = $1;
               if(!($totalram =~ /\d+/)) { # so above still didn't work, trying even harder
                   $x = `grep memory /var/run/dmesg.boot`; # this should work on freeBSD
                   $x =~ /avail memory = (\d+)/;
                   $totalram = int($1 / 1000000000);
                   if($totalram == 0) {
               	$totalram = "";
                   }
               }
               if(!($totalram =~ /\d+/)) { # so above didn't work, trying harder
                   $x = `top -l 1 | grep free`;  # this should work on a mac
                   $x =~ /(\d+)(.)\s+used, (\d+)(.) free/;
                   $used = $1;
                   $type1 = $2;
                   $free = $3;
                   $type2 = $4;
                   if($type1 eq "K" || $type1 eq "k") {
                    	$used = int($used / 1000000);
                   }
                   if($type2 eq "K" || $type2 eq "k") {
                    	$free = int($free / 1000000);
                        }
                   if($type1 eq "M" || $type1 eq "m") {
               	$used = int($used / 1000);
                   }
                   if($type2 eq "M" || $type2 eq "m") {
               	$free = int($free / 1000);
                   }
                   $totalram = $used + $free;
                   if($totalram == 0) {
               	$totalram = "";
                   }
               }
               if(!($totalram =~ /\d+/)) { # so above didn't work, warning user
                   $did_not_figure_out_ram = "true";
                   print "\nWarning: I could not determine how much RAM you have.  If you have less\nthan eight gigs per chunk and a large genome, or a lot of reads, then you shoud specify that\nusing the -ram option, or this might not work.\nFor a genome like human you'll need at least 5 or 6 Gb per chunk.\n\n";
                   $RAMperchunk = 6;
               } else {
                   $RAMperchunk = int($totalram / $numchunks);
                   if($RAMperchunk == 1) {
                       print "\nWarning: you have one gig of RAM per chunk.  If you\nhave a large genome or a lot of reads then this might not work.\nFor a genome like human you'll need at least 5 or 6 Gb per chunk.\n\n";
                   }elsif($RAMperchunk == 0) {
                       print "\nWarning: you have less than one gig of RAM per chunk.  If you\nhave a large genome or a lot of reads then this might not work.\nFor a genome like human you'll need at least 5 or 6 Gb per chunk.\n\n";
                      $RAMperchunk = 1;
                   }elsif($RAMperchunk < 7 && $RAMperchunk > 1) {
                       print "\nWarning: you have only $RAMperchunk gigs of RAM per chunk.  If you\nhave a large genome or a lot of reads then this might not work.\nFor a genome like human you'll need at least 5 or 6 Gb per chunk.\n\n";
                   }
              }
              $ram = $RAMperchunk;
          }
          if($did_not_figure_out_ram eq "false") {
              if($RAMperchunk >= 7) {
                  print STDERR "It seems like you have $totalram Gb of RAM on your machine.\n";
                  print STDERR "That's a generous amount, so unless you have too much other\nstuff running, RAM should not be a problem.\n";
              } else {
                  print STDERR "It seems like you have $totalram Gb of RAM on your machine.\n";
              }
              sleep(3);
              if($RAMperchunk >= 6) {
                   print STDERR "\nI'm going to try to use about $ram Gb of RAM per chunk.  Seems like that should work.\nIf that fails, try using the -ram option to lower it.  For a genome like human, you're\ngoing to need at least 5 or 6 Gb per chunk.\n\n";
              } else {
                   print STDERR "\nI'm going to try to use about $ram Gb of RAM per chunk.\nIf that fails, try using the -ram option to lower it.  For a genome like human, you're\ngoing to need at least 5 or 6 Gb per chunk.\n\n";
              }
          } else {
              print STDERR "\nI'm going to try to use about $ram Gb of RAM per chunk.  I couldn't figure out much you have so that's a (hopeful) guess.\nIf this fails, try using the -ram option to lower it.  For a genome like human, you're\ngoing to need at least 5 or 6 Gb per chunk.\n\n";
          }
     }
}

$check = `ps x | grep RUM_runner.pl`;
@a = split(/\n/,$check);
$CNT=0;
for($i=0; $i<@a; $i++) {
    chomp($a[$i]);
    next if $a[$i]=~ /screen.*RUM_runner/i;
    $a[$i] =~ s/.*RUM_runner.pl *//;
    @b = split(/ +/,$a[$i]);
    if($b[2] eq $output_dir || $b[2] eq $ARGV[2]) {
	$CNT++;
	if($CNT > 1) {
	    die "\nERROR: You seem to already have an instance of RUM_runner.pl running on the\nsame working directory.  This will cause collisions of the temporary files.\n\nExiting.\n\nTry killing by running the same command with -kill.\nIf that doesn't work use kill -9 on the process ID.\n\n";
	}
    }
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
if($kill eq "false") {
    sleep(2);
    print STDERR "Checking for phantom processes from prior runs that might need to be killed.\n\n";
    $outdir = $output_dir;
    $str = `ps x | grep $outdir | grep -v RUM_runner.pl`;
    @candidates = split(/\n/,$str);
    $cleanedflag = 0;
    for($i=0; $i<@candidates; $i++) {
	if($candidates[$i] =~ /^\s*(\d+)\s.*(\s|\/)$outdir\/pipeline.\d+.sh/) {
	    $pid = $1;
	    print STDERR "killing $pid\n";
	    `kill -9 $pid`;
            $cleanedflag = 1;
	}
    }
    for($i=0; $i<@candidates; $i++) {
	if($candidates[$i] =~ /^\s*(\d+)\s.*(\s|\/)$outdir(\s|\/)/) {
	    if(!($candidates[$i] =~ /pipeline.\d+.sh/)) {
		$pid = $1;
		print STDERR "killing $pid\n";
		`kill -9 $pid`;
                $cleanedflag = 1;
	    }
	}
    }
}
if($cleanedflag == 1) {
    sleep(2);
    print STDERR "OK there was some cleaning up to do, hopefully that worked.\n\n";
}
sleep(2);

for($i=1; $i<=$numchunks; $i++) {
    $logfile = "$output_dir/rum.log_chunk.$i";
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
$lib = <INFILE>;
$lib =~ s!/$!!;
chomp($lib);
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
if(($readsfile =~ /,,,/) && ($paired_end eq "true") && ($postprocess eq "false")) {
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

    `perl $scripts_dir/parse2fasta.pl $a[0] $a[1] | head -10000 > $output_dir/reads_temp.fa`;
    `perl $scripts_dir/fastq2qualities.pl $a[0] $a[1] | head -10000 > $output_dir/quals_temp.fa`;
    $X = `head -2 $output_dir/quals_temp.fa | tail -1`;
    if($X =~ /\S/ && !($X =~ /Sorry, can't figure these files out/)) {
        open(RFILE, "$output_dir/reads_temp.fa");
        open(QFILE, "$output_dir/quals_temp.fa");
        while($linea = <RFILE>) {
            $lineb = <QFILE>;
            $line1 = <RFILE>;
            $line2 = <QFILE>;
            chomp($line1);
            chomp($line2);
            if(length($line1) != length($line2)) {
               $readlength = length($line1);
               $quallength = length($line2);
               die "ERROR: It seems your read lengths differ from your quality string lengths.\nCheck line:\n$linea$line1\n$lineb$line2\n\n";
           }
        }
    }
    unlink("$output_dir/reads_temp.fa");
    unlink("$output_dir/quals_temp.fa");

    print STDERR "Reformatting reads file... please be patient.\n";
    `perl $scripts_dir/parse2fasta.pl $a[0] $a[1] > $output_dir/reads.fa`;
    `perl $scripts_dir/fastq2qualities.pl $a[0] $a[1] > $output_dir/quals.fa`;
    $X = `head -2 $output_dir/quals.fa | tail -1`;
    if($X =~ /\S/ && !($X =~ /Sorry, can't figure these files out/)) {
	 $quals = "true";
    }
    $qualsfile = "$output_dir/quals.fa";
    $readsfile = "$output_dir/reads.fa";
}

if($postprocess eq "true") {
    $readsfile = "$output_dir/reads.fa";
    if(!(-e $readsfile)) {
        $readsfile = $ARGV[1];
    }
    $qualsfile = "$output_dir/quals.fa";
    $quals = "false";
    if(-e $qualsfile) {
        $X = `head -2 $output_dir/quals.fa | tail -1`;
        if($X =~ /\S/ && !($X =~ /Sorry, can't figure these files out/)) {
	    $quals = "true"
        }
    }
}

$head = `head -2 $readsfile | tail -1`;
chomp($head);
$rl = length($head);
$tail = `tail -2 $readsfile | head -1`;
$tail =~ /seq.(\d+)/s;
$nr = $1;

if($minlength == 0) {
	if($rl < 80) {
	    if($match_length_cutoff == 0) {
		$match_length_cutoff = 35;
	    }
	} else {
	    if($match_length_cutoff == 0) {
		$match_length_cutoff = 50;
	    }
	}
	if($min_size_intersection_allowed >= .8 * $rl) {
	    if($match_length_cutoff == 0) {
		$match_length_cutoff = int(.6 * $rl);
	    }
	}
} else {
	$match_length_cutoff = $minlength;
}

open(LOGFILE, ">$output_dir/rum.log_master");
print LOGFILE "config file: $configfile\n";
print LOGFILE "readsfile: $ARGV[1]\n";
print LOGFILE "output_dir: $output_dir\n";
print LOGFILE "readlength = $rl\n";
$NR = &format_large_int($nr);
if($paired_end eq 'false') {
    print LOGFILE "number of reads: $NR\n";
} else {
    print LOGFILE "number of read pairs: $NR\n";
}
print LOGFILE "minimum length alignment to report = $match_length_cutoff\n";
print LOGFILE "numchunks: $numchunks\n";
print LOGFILE "name: $name\n";
print LOGFILE "paired_end: $paired_end\n";
print LOGFILE "fast: $fast\n";
print LOGFILE "limitBowtieNU: $limitNU\n";
print LOGFILE "limitNU: $limitNUhard\n";
print LOGFILE "dna: $dna\n";
print LOGFILE "qsub: $qsub\n";
print LOGFILE "blat minidentity: $minidentity\n";
print LOGFILE "output junctions: $junctions\n";
print LOGFILE "output quantified values: $quantify\n";
print LOGFILE "strand specific: $strandspecific\n";

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
   if($minlength > $readlength) {
       die "Error: you specified a minimum length alignment to report as '$minlength', however\nyour read length is only $readlength\n";
   }
} else {
   $readlength = "v";
}

print LOGFILE "\nstart: $date\n";

$head = `head -4 $readsfile`;
$head =~ /seq.(\d+)(.).*seq.(\d+)(.)/s;
$num1 = $1;
$type1 = $2;
$num2 = $3;
$type2 = $4;
if($postprocess eq "false") {
    if($paired_end eq 'false') {
        if($type1 ne "a") {
	   print STDERR "Reformatting reads file... please be patient.\n";
	   `perl $scripts_dir/parse2fasta.pl $readsfile > $output_dir/reads.fa`;
	   `perl $scripts_dir/fastq2qualities.pl $readsfile > $output_dir/quals.fa`;
	   $X = `head -2 $output_dir/quals.fa | tail -1`;
           if($X =~ /\S/ && !($X =~ /Sorry, can't figure these files out/)) {
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
}

if($quals_specified eq 'true') {
    open(TESTIN, "$output_dir/$quals_file") or die "\nError: cannot open '$quals_file' for reading, it should be in the '$output_dir' directory.\n\n";
    close(TESTIN);
    $qualsfile = "$output_dir/$quals_file";
    $quals = "true";
}

if($postprocess eq "true") {
    $head = `head -4 $readsfile`;
    $head =~ /seq.(\d+)(.).*seq.(\d+)(.)/s;
    $num1 = $1;
    $type1 = $2;
    $num2 = $3;
    $type2 = $4;
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
    print STDERR "very large, in this case it's recommended to run with the -limitBowtieNU option.  You\n";
    print STDERR "can watch the files that start with 'X' and 'Y' to see if they are growing\n";
    print STDERR "larger than 10 gigabytes per million reads at which point you might want to use.\n";
    print STDERR "-limitNU\n\n";
}

if($blatonly eq "true" && $dna eq "true") {
    $dna = "false";  # because blat only is dna only anyway, and setting them both breaks things below
}

if($postprocess eq "false") {
    $pipeline_template = `cat $lib/pipeline_template.sh`;
    if($cleanup eq 'false') {
        $pipeline_template =~ s/^.*unlink.*$//mg;
        $pipeline_template =~ s!if . -f OUTDIR.RUM_NU_temp3.CHUNK .\nthen\n\nfi\n!!gs;
    }
    if($dna eq "true") {
        $pipeline_template =~ s/# cp /cp /gs;
        $pipeline_template =~ s/xxx1.*xxx2//s;
    }
    if($blatonly eq "true") {
        $pipeline_template =~ s/xxx0.*xxx2//s;
        $pipeline_template =~ s!# cp OUTDIR/GU.CHUNK OUTDIR/BowtieUnique.CHUNK!echo `` >> OUTDIR/BowtieUnique.CHUNK!s;
        $pipeline_template =~ s!# cp OUTDIR/GNU.CHUNK OUTDIR/BowtieNU.CHUNK!echo `` >> OUTDIR/BowtieNU.CHUNK!s;
    }
    if($fasta_already_fragmented eq "false") {
        print STDERR "Splitting files ... please be patient.\n\n";
        $qualflag = 0;
        $x = breakup_file($readsfile, $numchunks);
        if($quals eq "true") {
            print STDERR "Half done splitting...\n\n";
            $qualflag = 1;
    	    $x = breakup_file($qualsfile, $numchunks);
        }
    } else {
        $readsfile =~ /([^\/]+)$/;
        $readsfile_nopath = $1;
        for($i=1; $i<=$numchunks; $i++) {
            $r = $readsfile_nopath . "." . $i;
            if(!(-e "$output_dir/$readsfile_nopath")) {
                die "\n-------------------------------------------------------------------\nError: You said the files were already broken up, but the file\n'$output_dir/$readsfile_nopath' does not seem to exist.\n\nNote: even if the files are already broken up, they still need to be\nin the <output dir> directory that you specified as '$output_dir'.\n-------------------------------------------------------------------\n\n";
            }
        }
	$quals = "true";
        for($i=1; $i<=$numchunks; $i++) {
            $qfile = "$output_dir/quals.$i";
            if(!(-e $qfile)) {
                $quals = "false";
                $i = $numchunks+1;
            }
            if($quals eq 'false') {
                print STDERR "Note: I did not find any quality files, I will assume you are\naware of that and just did not have them.\n\n";
            }
        }
    }

    print STDERR "Reads fasta file already fragmented: $fasta_already_fragmented\n";
    print STDERR "Number of Chunks: $numchunks\n";
    print STDERR "Reads File: $readsfile\n";
    print STDERR "Paired End: $paired_end\n";

    $readsfile =~ s!.*/!!;
    $readsfile = $output_dir . "/" . $readsfile;
    $t = `tail -2 $readsfile`;
    $t =~ /seq.(\d+)/;
    $NumSeqs = $1;
    $f = &format_large_int($NumSeqs);
    if($paired_end eq 'true') {
       print STDERR "Number of Read Pairs: $f\n";
    } else {
       print STDERR "Number of Reads: $f\n";
    }
    
    print STDERR "\nEverything seems okay, I am going to fire off the job.\n\n";
    
    for($i=1; $i<=$numchunks; $i++) {
        $pipeline_file = $pipeline_template;
        if($limitNUhard eq "true") {
    	   $pipeline_file =~ s!LIMITNUCUTOFF!$NU_limit!gs;
        } else {
    	   $pipeline_file =~ s!perl SCRIPTSDIR/limit_NU.pl OUTDIR/RUM_NU_temp3.CHUNK LIMITNUCUTOFF > OUTDIR/RUM_NU.CHUNK\n!mv OUTDIR/RUM_NU_temp3.CHUNK OUTDIR/RUM_NU.CHUNK\n!gs;
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
        if($strandspecific eq 'true') {
           $pipeline_file =~ s/STRAND1s/-strand p/gs;
           $pipeline_file =~ s/quant.S1s/quant.ps/gs;
           $pipeline_file =~ s/STRAND2s/-strand m/gs;
           $pipeline_file =~ s/quant.S2s/quant.ms/gs;
           $pipeline_file =~ s/STRAND1a/-strand p -anti/gs;
           $pipeline_file =~ s/quant.S1a/quant.pa/gs;
           $pipeline_file =~ s/STRAND2a/-strand m -anti/gs;
           $pipeline_file =~ s/quant.S2a/quant.ma/gs;
        } else {
           $pipeline_file =~ s/STRAND1s//sg;
           $pipeline_file =~ s/quant.S1s/quant/sg;
           $pipeline_file =~ s/[^\n]+quant.S2s[^\n]+\n//sg;
           $pipeline_file =~ s/[^\n]+quant.S1a[^\n]+\n//sg;
           $pipeline_file =~ s/[^\n]+quant.S2a[^\n]+\n//sg;
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
        if($minlength > 0) {
    	   $pipeline_file =~ s!MATCHLENGTHCUTOFF!-match_length_cutoff $minlength!gs;
        } else {
    	   $pipeline_file =~ s!MATCHLENGTHCUTOFF!!gs;
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
           $ofile = $output_dir . "/chunk.$i" . ".o";
           $efile = $output_dir . "/chunk.$i" . ".e";
    	   `qsub -l mem_free=7G -o $ofile -e $efile $output_dir/$outfile`;
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
        if($qsub eq "false" && $blatonly eq "false") {
    	    print STDERR "\nThe next thing to print here will (probably) be the status reports from bowtie for each chunk.\n";
    	    print STDERR "     * Don't be alarmed by the number of reads that 'failed to align'\n       that's just referring to one stage, the end result will be better.\n\n";
        }
    } else {
        print STDERR "\nThe job has been initiated, now the long wait...\n";
        sleep(2);
        if($qsub eq "false" && $blatonly eq "false") {
    	    print STDERR "\nThe next thing to print here will (probably) be the status reports from bowtie for each chunk.\n";
    	    print STDERR "     * Don't be alarmed by the number of reads that 'failed to align'\n       that's just referring to stage, the end result will be better.\n\n";
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
            $logfile = "$output_dir/rum.log_chunk.$i";
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
}
if($postprocess eq "false") {
     print STDERR "\nWhew, all chunks have finished.\n\nNext I will merge everything, create the coverage plots and\ncalculate the quantified values, etc.  This could take some time...\n\n";
} else {
     print STDERR "\nOK, will now merge everything, create the coverage plots and\ncalculate the quantified values, etc.  This could take some time...\n\n";
}

$t = `tail -2 $readsfile`;
$t =~ /seq.(\d+)/s;
$NumSeqs = $1;

if($nocat eq "false") {
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
    for($i=1; $i<=$numchunks; $i++) {
       open(SAMHEADER, "$output_dir/sam_header.$i");
       while($line = <SAMHEADER>) {
           chomp($line);
           $line =~ /SN:([^\s]+)\s/;
           $samheader{$1}=$line;
       }
       close(SAMHEADER);
    }
    open(SAMOUT, ">$output_dir/RUM.sam");
    foreach $key (sort {cmpChrs($a,$b)} keys %samheader) {
        $shout = $samheader{$key};
        print SAMOUT "$shout\n";
    }
    close(SAMOUT);
    for($i=1; $i<=$numchunks; $i++) {
        $x = `cat $output_dir/RUM.sam.$i >> $output_dir/RUM.sam`;
    }
}

print LOGFILE "finished creating RUM_Unique/RUM_NU/RUM.sam: $date\n";

if($cleanup eq 'true') {
   print STDERR "\nCleaning up some temp files...\n\n";
   `yes|rm $output_dir/RUM.sam.* $output_dir/RUM_Unique.* $output_dir/RUM_NU.* $output_dir/sam_header.* $output_dir/reads.fa.*`;
   if(-e "$output_dir/quals.fa.1") {
       `yes|rm $output_dir/quals.fa.*`;
   }
}

print LOGFILE "starting the post processing: $date\n";
$PPlog = "postprocessing_$name" . ".log";
$shellscript = "#!/bin/sh\n";
if($NumSeqs =~ /(\d+)/) {
    $shellscript = $shellscript . "perl $scripts_dir/count_reads_mapped.pl $output_dir/RUM_Unique $output_dir/RUM_NU -minseq 1 -maxseq $NumSeqs > $output_dir/mapping_stats.txt\n";
} else {
    $shellscript = $shellscript . "perl $scripts_dir/count_reads_mapped.pl $output_dir/RUM_Unique $output_dir/RUM_NU -minseq 1 > $output_dir/mapping_stats.txt\n";
}
if($quantify eq "true") {
    if($strandspecific eq 'true') {
        $shellscript = $shellscript . "perl $scripts_dir/merge_quants.pl $output_dir $numchunks $output_dir/feature_quantifications.ps -strand ps\n";
        $shellscript = $shellscript . "perl $scripts_dir/merge_quants.pl $output_dir $numchunks $output_dir/feature_quantifications.ms -strand ms\n";
        $shellscript = $shellscript . "perl $scripts_dir/merge_quants.pl $output_dir $numchunks $output_dir/feature_quantifications.pa -strand pa\n";
        $shellscript = $shellscript . "perl $scripts_dir/merge_quants.pl $output_dir $numchunks $output_dir/feature_quantifications.ma -strand ma\n";
        $shellscript = $shellscript . "perl $scripts_dir/merge_quants_strandspecific.pl $output_dir/feature_quantifications.ps $output_dir/feature_quantifications.ms $output_dir/feature_quantifications.pa $output_dir/feature_quantifications.ma $gene_annot_file $output_dir/feature_quantifications_$name\n";

    } else {
        $shellscript = $shellscript . "perl $scripts_dir/merge_quants.pl $output_dir $numchunks $output_dir/feature_quantifications_$name\n";
    }
}
$shellscript = $shellscript . "echo sorting RUM_Unique > $output_dir/$PPlog\n";
$shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
$shellscript = $shellscript . "perl $scripts_dir/sort_RUM_by_location.pl $output_dir/RUM_Unique $output_dir/RUM_Unique.sorted -ram $ram >> $output_dir/mapping_stats.txt\n";
$shellscript = $shellscript . "echo sorting RUM_NU >> $output_dir/$PPlog\n";
$shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
$shellscript = $shellscript . "perl $scripts_dir/sort_RUM_by_location.pl $output_dir/RUM_NU $output_dir/RUM_NU.sorted -ram $ram >> $output_dir/mapping_stats.txt\n";
$shellscript = $shellscript . "echo making coverage plots >> $output_dir/$PPlog\n";
$shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
$shellscript = $shellscript . "perl $scripts_dir/rum2cov.pl $output_dir/RUM_Unique.sorted $output_dir/RUM_Unique.cov -name \"$name Unique Mappers\"\n";
$shellscript = $shellscript . "perl $scripts_dir/rum2cov.pl $output_dir/RUM_NU.sorted $output_dir/RUM_NU.cov -name \"$name Non-Unique Mappers\"\n";
if($strandspecific eq 'true') {
      # breakup RUM_Unique and RUM_NU files into plus and minus
      $shellscript = $shellscript . "perl $scripts_dir/breakup_RUM_files_by_strand.pl $output_dir/RUM_Unique.sorted $output_dir/RUM_Unique.sorted.plus $output_dir/RUM_Unique.sorted.minus\n";
      $shellscript = $shellscript . "perl $scripts_dir/breakup_RUM_files_by_strand.pl $output_dir/RUM_NU.sorted $output_dir/RUM_NU.sorted.plus $output_dir/RUM_NU.sorted.minus\n";
      # run rum2cov on all four files
      $shellscript = $shellscript . "perl $scripts_dir/rum2cov.pl $output_dir/RUM_Unique.sorted.plus $output_dir/RUM_Unique.plus.cov -name \"$name Unique Mappers Plus Strand\"\n";
      $shellscript = $shellscript . "perl $scripts_dir/rum2cov.pl $output_dir/RUM_Unique.sorted.minus $output_dir/RUM_Unique.minus.cov -name \"$name Unique Mappers Minus Strand\"\n";
      $shellscript = $shellscript . "perl $scripts_dir/rum2cov.pl $output_dir/RUM_NU.sorted.plus $output_dir/RUM_NU.plus.cov -name \"$name Non-Unique Mappers Plus Strand\"\n";
      $shellscript = $shellscript . "perl $scripts_dir/rum2cov.pl $output_dir/RUM_NU.sorted.minus $output_dir/RUM_NU.minus.cov -name \"$name Non-Unique Mappers Minus Strand\"\n";
}
if($junctions eq "true") {
   $shellscript = $shellscript . "echo starting to compute junctions >> $output_dir/$PPlog\n";
   $shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
   if($altgenes eq "true") {
       $shellscript = $shellscript . "perl $scripts_dir/make_RUM_junctions_file.pl $output_dir/RUM_Unique $output_dir/RUM_NU $genomefa $altgene_file $output_dir/junctions_all.rum $output_dir/junctions_all.bed $output_dir/junctions_high-quality.bed -faok\n";
   } else {
       $shellscript = $shellscript . "perl $scripts_dir/make_RUM_junctions_file.pl $output_dir/RUM_Unique $output_dir/RUM_NU $genomefa $gene_annot_file $output_dir/junctions_all.rum $output_dir/junctions_all.bed $output_dir/junctions_high-quality.bed -faok\n";
   }
}
$shellscript = $shellscript . "echo finished >> $output_dir/$PPlog\n";
$shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
$str = "postprocessing_$name" . ".sh";
open(OUTFILE2, ">$output_dir/$str");
print OUTFILE2 $shellscript;
close(OUTFILE2);

if($qsub eq "true") {
    $ofile = $output_dir . "/postprocessing" . ".o";
    $efile = $output_dir . "/postprocessing" . ".e";
    `qsub -l mem_free=7G -o $ofile -e $efile $output_dir/$str`;
} else {
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

if($cleanup eq 'true') {
   `yes|rm $output_dir/quant.*`;
   `yes|rm $output_dir/pipeline.*`;
   if($strandspecific eq 'true') {
      `yes|rm $output_dir/feature_quantifications.ps`;
      `yes|rm $output_dir/feature_quantifications.ms`;
      `yes|rm $output_dir/feature_quantifications.pa`;
      `yes|rm $output_dir/feature_quantifications.ma`;
   }
}

print STDERR "\nOkay, all finished.\n\n";

$date = `date`;
print LOGFILE "pipeline finished: $date\n";
close(LOGFILE);

sub breakup_file () {
    ($FILE, $numpieces) = @_;

    open(INFILE, $FILE) or die "\nError: Cannot open '$FILE' for reading.\n\n";
    $tail = `tail -2 $FILE | head -1`;
    $tail =~ /seq.(\d+)/s;
    $numseqs = $1;
    $piecesize = int($numseqs / $numpieces);

    $t = `tail -2 $FILE`;
    $t =~ /seq.(\d+)/s;
    $NS = $1;
    $piecesize2 = &format_large_int($piecesize);
    if(!($FILE =~ /qual/)) {
	if($numchunks > 1) {
	    print LOGFILE "processing in $numchunks pieces of approx $piecesize2 reads each\n";
	} else {
	    $NS2 = &format_large_int($NS);
	    print LOGFILE "processing in one piece of $NS2 reads\n";
	}
    }
    if($piecesize % 2 == 1) {
	$piecesize++;
    }
    $bflag = 0;

    $F2 = $FILE;
    $F2 =~ s!.*/!!;

    if($paired_end eq 'true') {
	$PS = $piecesize * 2;
    } else {
	$PS = $piecesize;
    }

    for($i=1; $i<$numpieces; $i++) {
	$outfilename = $output_dir . "/" . $F2 . "." . $i;

	open(OUTFILE, ">$outfilename");
	for($j=0; $j<$PS; $j++) {
	    $line = <INFILE>;
	    chomp($line);
	    if($qualflag == 0) {
		$line =~ s/[^ACGTNab]$//s;
	    }
	    print OUTFILE "$line\n";
	    $line = <INFILE>;
	    chomp($line);
	    if($qualflag == 0) {
		$line =~ s/[^ACGTNab]$//s;
	    }
	    print OUTFILE "$line\n";
	}
	close(OUTFILE);
    }
    $outfilename = $output_dir . "/" . $F2 . "." . $numpieces;

    open(OUTFILE, ">$outfilename");
    while($line = <INFILE>) {
	print OUTFILE $line;
    }
    close(OUTFILE);
    return 0;
}

sub checkstatus () {
    ($CHUNK) = @_;
    $log = `cat $output_dir/rum.log_chunk.$CHUNK`;
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

sub merge() {
    $tempfilename1 = $CHR[$cnt] . "_temp.0";
    $tempfilename2 = $CHR[$cnt] . "_temp.1";
    $tempfilename3 = $CHR[$cnt] . "_temp.2";
    open(TEMPMERGEDOUT, ">$tempfilename3");
    open(TEMPIN1, $tempfilename1);
    open(TEMPIN2, $tempfilename2);
    $mergeFLAG = 0;
    getNext1();
    getNext2();
    while($mergeFLAG < 2) {
	chomp($out1);
	chomp($out2);
	if($start1 < $start2) {
	    if($out1 =~ /\S/) {
		print TEMPMERGEDOUT "$out1\n";
	    }
	    getNext1();
	} elsif($start1 == $start2) {
	    if($end1 <= $end2) {
		if($out1 =~ /\S/) {
		    print TEMPMERGEDOUT "$out1\n";
		}
		getNext1();
	    } else {
		if($out2 =~ /\S/) {
		    print TEMPMERGEDOUT "$out2\n";
		}
		getNext2();
	    }
	} else {
	    if($out2 =~ /\S/) {
		print TEMPMERGEDOUT "$out2\n";
	    }
	    getNext2();
	}
    }
    close(TEMPMERGEDOUT);
    `mv $tempfilename3 $tempfilename1`;
    unlink($tempfilename2);
}

sub getNext1 () {
    $line1 = <TEMPIN1>;
    chomp($line1);
    if($line1 eq '') {
	$mergeFLAG++;
	$start1 = 1000000000000;  # effectively infinity, no chromosome should be this large;
	return "";
    }
    @a = split(/\t/,$line1);
    $a[2] =~ /^(\d+)-/;
    $start1 = $1;
    if($a[0] =~ /a/ && $separate eq "false") {
	$a[0] =~ /(\d+)/;
	$seqnum1 = $1;
	$line2 = <TEMPIN1>;
	chomp($line2);
	@b = split(/\t/,$line2);
	$b[0] =~ /(\d+)/;
	$seqnum2 = $1;
	if($seqnum1 == $seqnum2 && $b[0] =~ /b/) {
	    if($a[3] eq "+") {
		$b[2] =~ /-(\d+)$/;
		$end1 = $1;
	    } else {
		$b[2] =~ /^(\d+)-/;
		$start1 = $1;
		$a[2] =~ /-(\d+)$/;
		$end1 = $1;
	    }
	    $out1 = $line1 . "\n" . $line2;
	} else {
	    $a[2] =~ /-(\d+)$/;
	    $end1 = $1;
	    # reset the file handle so the last line read will be read again
	    $len = -1 * (1 + length($line2));
	    seek(TEMPIN1, $len, 1);
	    $out1 = $line1;
	}
    } else {
	$a[2] =~ /-(\d+)$/;
	$end1 = $1;
	$out1 = $line1;
    }
}

sub getNext2 () {
    $line1 = <TEMPIN2>;
    chomp($line1);
    if($line1 eq '') {
	$mergeFLAG++;
	$start2 = 1000000000000;  # effectively infinity, no chromosome should be this large;
	return "";
    }
    @a = split(/\t/,$line1);
    $a[2] =~ /^(\d+)-/;
    $start2 = $1;
    if($a[0] =~ /a/ && $separate eq "false") {
	$a[0] =~ /(\d+)/;
	$seqnum1 = $1;
	$line2 = <TEMPIN2>;
	chomp($line2);
	@b = split(/\t/,$line2);
	$b[0] =~ /(\d+)/;
	$seqnum2 = $1;
	if($seqnum1 == $seqnum2 && $b[0] =~ /b/) {
	    if($a[3] eq "+") {
		$b[2] =~ /-(\d+)$/;
		$end2 = $1;
	    } else {
		$b[2] =~ /^(\d+)-/;
		$start2 = $1;
		$a[2] =~ /-(\d+)$/;
		$end2 = $1;
	    }
	    $out2 = $line1 . "\n" . $line2;
	} else {
	    $a[2] =~ /-(\d+)$/;
	    $end2 = $1;
	    # reset the file handle so the last line read will be read again
	    $len = -1 * (1 + length($line2));
	    seek(TEMPIN2, $len, 1);
	    $out2 = $line1;
	}
    } else {
	$a[2] =~ /-(\d+)$/;
	$end2 = $1;
	$out2 = $line1;
    }
}

sub format_large_int () {
    ($int) = @_;
    @a = split(//,"$int");
    $j=0;
    $newint = "";
    $n = @a;
    for($i=$n-1;$i>=0;$i--) {
	$j++;
	$newint = $a[$i] . $newint;
	if($j % 3 == 0) {
	    $newint = "," . $newint;
	}
    }
    $newint =~ s/^,//;
    return $newint;
}

sub cmpChrs () {
    $a2_c = lc($b);
    $b2_c = lc($a);
    if($a2_c =~ /^\d+$/ && !($b2_c =~ /^\d+$/)) {
        return 1;
    }
    if($b2_c =~ /^\d+$/ && !($a2_c =~ /^\d+$/)) {
        return -1;
    }
    if($a2_c =~ /^[ivxym]+$/ && !($b2_c =~ /^[ivxym]+$/)) {
        return 1;
    }
    if($b2_c =~ /^[ivxym]+$/ && !($a2_c =~ /^[ivxym]+$/)) {
        return -1;
    }
    if($a2_c eq 'm' && ($b2_c eq 'y' || $b2_c eq 'x')) {
        return -1;
    }
    if($b2_c eq 'm' && ($a2_c eq 'y' || $a2_c eq 'x')) {
        return 1;
    }
    if($a2_c =~ /^[ivx]+$/ && $b2_c =~ /^[ivx]+$/) {
        $a2_c = "chr" . $a2_c;
        $b2_c = "chr" . $b2_c;
    }
    if($a2_c =~ /$b2_c/) {
	return -1;
    }
    if($b2_c =~ /$a2_c/) {
	return 1;
    }
    # dealing with roman numerals starts here
    if($a2_c =~ /chr([ivx]+)/ && $b2_c =~ /chr([ivx]+)/) {
	$a2_c =~ /chr([ivx]+)/;
	$a2_roman = $1;
	$b2_c =~ /chr([ivx]+)/;
	$b2_roman = $1;
	$a2_arabic = arabic($a2_roman);
    	$b2_arabic = arabic($b2_roman);
	if($a2_arabic > $b2_arabic) {
	    return -1;
	} 
	if($a2_arabic < $b2_arabic) {
	    return 1;
	}
	if($a2_arabic == $b2_arabic) {
	    $tempa = $a2_c;
	    $tempb = $b2_c;
	    $tempa =~ s/chr([ivx]+)//;
	    $tempb =~ s/chr([ivx]+)//;
	    undef %temphash;
	    $temphash{$tempa}=1;
	    $temphash{$tempb}=1;
	    foreach $tempkey (sort {cmpChrs($a,$b)} keys %temphash) {
		if($tempkey eq $tempa) {
		    return 1;
		} else {
		    return -1;
		}
	    }
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
    # roman numerals ends here
    if($a2_c =~ /chr(\d+)$/ && $b2_c =~ /chr.*_/) {
        return 1;
    }
    if($b2_c =~ /chr(\d+)$/ && $a2_c =~ /chr.*_/) {
        return -1;
    }
    if($a2_c =~ /chr([a-z])$/ && $b2_c =~ /chr.*_/) {
        return 1;
    }
    if($b2_c =~ /chr([a-z])$/ && $a2_c =~ /chr.*_/) {
        return -1;
    }
    if($a2_c =~ /chr(\d+)/) {
        $numa = $1;
        if($b2_c =~ /chr(\d+)/) {
            $numb = $1;
            if($numa < $numb) {return 1;}
	    if($numa > $numb) {return -1;}
	    if($numa == $numb) {
		$tempa = $a2_c;
		$tempb = $b2_c;
		$tempa =~ s/chr\d+//;
		$tempb =~ s/chr\d+//;
		undef %temphash;
		$temphash{$tempa}=1;
		$temphash{$tempb}=1;
		foreach $tempkey (sort {cmpChrs($a,$b)} keys %temphash) {
		    if($tempkey eq $tempa) {
			return 1;
		    } else {
			return -1;
		    }
		}
	    }
        } else {
            return 1;
        }
    }
    if($a2_c =~ /chrx(.*)/ && ($b2_c =~ /chr(y|m)$1/)) {
	return 1;
    }
    if($b2_c =~ /chrx(.*)/ && ($a2_c =~ /chr(y|m)$1/)) {
	return -1;
    }
    if($a2_c =~ /chry(.*)/ && ($b2_c =~ /chrm$1/)) {
	return 1;
    }
    if($b2_c =~ /chry(.*)/ && ($a2_c =~ /chrm$1/)) {
	return -1;
    }
    if($a2_c =~ /chr\d/ && !($b2_c =~ /chr[^\d]/)) {
	return 1;
    }
    if($b2_c =~ /chr\d/ && !($a2_c =~ /chr[^\d]/)) {
	return -1;
    }
    if($a2_c =~ /chr[^xy\d]/ && (($b2_c =~ /chrx/) || ($b2_c =~ /chry/))) {
        return -1;
    }
    if($b2_c =~ /chr[^xy\d]/ && (($a2_c =~ /chrx/) || ($a2_c =~ /chry/))) {
        return 1;
    }
    if($a2_c =~ /chr(\d+)/ && !($b2_c =~ /chr(\d+)/)) {
        return 1;
    }
    if($b2_c =~ /chr(\d+)/ && !($a2_c =~ /chr(\d+)/)) {
        return -1;
    }
    if($a2_c =~ /chr([a-z])/ && !($b2_c =~ /chr(\d+)/) && !($b2_c =~ /chr[a-z]+/)) {
        return 1;
    }
    if($b2_c =~ /chr([a-z])/ && !($a2_c =~ /chr(\d+)/) && !($a2_c =~ /chr[a-z]+/)) {
        return -1;
    }
    if($a2_c =~ /chr([a-z]+)/) {
        $letter_a = $1;
        if($b2_c =~ /chr([a-z]+)/) {
            $letter_b = $1;
            if($letter_a lt $letter_b) {return 1;}
	    if($letter_a gt $letter_b) {return -1;}
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
    if($a2_c le $b2_c) {
	return 1;
    }
    if($b2_c le $a2_c) {
	return -1;
    }


    return 1;
}

sub isroman($) {
    $arg = shift;
    $arg ne '' and
      $arg =~ /^(?: M{0,3})
                (?: D?C{0,3} | C[DM])
                (?: L?X{0,3} | X[LC])
                (?: V?I{0,3} | I[VX])$/ix;
}

sub arabic($) {
    $arg = shift;
    %roman2arabic = qw(I 1 V 5 X 10 L 50 C 100 D 500 M 1000);
    %roman_digit = qw(1 IV 10 XL 100 CD 1000 MMMMMM);
    @figure = reverse sort keys %roman_digit;
    $roman_digit{$_} = [split(//, $roman_digit{$_}, 2)] foreach @figure;
    isroman $arg or return undef;
    ($last_digit) = 1000;
    $arabic = 0;
    ($arabic);
    foreach (split(//, uc $arg)) {
        ($digit) = $roman2arabic{$_};
        $arabic -= 2 * $last_digit if $last_digit < $digit;
        $arabic += ($last_digit = $digit);
    }
    $arabic;
}

sub Roman($) {
    $arg = shift;
    %roman2arabic = qw(I 1 V 5 X 10 L 50 C 100 D 500 M 1000);
    %roman_digit = qw(1 IV 10 XL 100 CD 1000 MMMMMM);
    @figure = reverse sort keys %roman_digit;
    $roman_digit{$_} = [split(//, $roman_digit{$_}, 2)] foreach @figure;
    0 < $arg and $arg < 4000 or return undef;
    $roman = "";
    ($x, $roman);
    foreach (@figure) {
        ($digit, $i, $v) = (int($arg / $_), @{$roman_digit{$_}});
        if (1 <= $digit and $digit <= 3) {
            $roman .= $i x $digit;
        } elsif ($digit == 4) {
            $roman .= "$i$v";
        } elsif ($digit == 5) {
            $roman .= $v;
        } elsif (6 <= $digit and $digit <= 8) {
            $roman .= $v . $i x ($digit - 5);
        } elsif ($digit == 9) {
            $roman .= "$i$x";
        }
        $arg -= $digit * $_;
        $x = $i;
    }
    $roman;
}

sub roman($) {
    lc Roman shift;
}
