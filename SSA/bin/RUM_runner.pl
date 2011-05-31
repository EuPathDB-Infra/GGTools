
# Written by Gregory R Grant
# University of Pennsylvania, 2010

$version = "1.09.  Released May 28th, 2011";

$| = 1;

if($ARGV[0] eq '-version' || $ARGV[0] eq '-v' || $ARGV[0] eq '--version' || $ARGV[0] eq '--v') {
    die "RUM version: $version\n";
}

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
RUM version: $version
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

<reads file(s)> :  1) For unpaired data, the single file of reads.
                   2) For paired data the files of forward and reverse reads,
                      separated by three commas ',,,' (with no spaces).

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

Options: There are many options, but RUM is typically run with the defaults. The
         option -kill is also quite useful to stop a run, because killing just
         the main program will not always kill the spawned processes.

       -strandspecific : If the data are strand specific, then you can use this
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

       -noclean   : do not remove the intermediate and temp files after finishing.

       -kill      : To kill a job, run with all the same parameters but add
                    -kill.  Note: it is not sufficient to just terminate
                    RUM_runner.pl, that will leave other phantom processes.
                    Use -kill instead.

       -ram n : On some systems it might not be able to determine the amount of
                RAM you have.  In that case, with this option you can specify
                the number of Gb of ram you want to dedicate to each chunk.
                This is rarely necessary and never necessary if you have at
                least 6 Gb per chunk.

       -version : Returns the current version.

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

$JID = int(rand(10000000)) + 10;

$configfile = $ARGV[0];
$readsfile = $ARGV[1];
$output_dir = $ARGV[2];
$output_dir =~ s!/$!!;
if(!(-d $output_dir)) {
    die "\nERROR: The directory '$output_dir' does not seem to exists...\n\n";
}
open(ERRORLOG, ">$output_dir/rum.error-log");
print ERRORLOG "\n--------------------\n";
print ERRORLOG "Job ID: $JID\n";
print ERRORLOG "--------------------\n";

$numchunks = $ARGV[3];
$NUMCHUNKS = $ARGV[3];
if($numchunks =~ /(\d+)s/) {
    $numchunks = $1;
}
$name = $ARGV[4];
if($name =~ /^-/) {
    print ERRORLOG "\nERROR: The name '$name' is invalid, probably you forgot a required argument\n\n";
    die "\nERROR: The name '$name' is invalid, probably you forgot a required argument\n\n";
}

$name_o = $ARGV[4];
$name =~ s/\s+/_/g;
$name =~ s/^[^a-zA-Z0-9_.-]//;
$name =~ s/[^a-zA-Z0-9_.-]$//g;
$name =~ s/[^a-zA-Z0-9_.-]/_/g;

if($name ne $name_o) {
    print "\nWarning: name changed from '$name_o' to '$name'.\n\n";
}
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
$ram = 6;
$user_ram = "false";
$nocat = "false";
$quals_specified = "false";
$strandspecific = "false";
$quantify = "false";
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
	if($ARGV[$i] eq "-strandspecific") {
	    $strandspecific = "true";
	    $optionrecognized = 1;
	}
        if($ARGV[$i] eq "-ram") {
	    $i++;
	    $ram = $ARGV[$i];
            $user_ram = "true";
            if($ARGV[$i] =~ /^\d+$/) {
	        $optionrecognized = 1;
	    }
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
	if($ARGV[$i] eq "-variable_read_lengths" || $ARGV[$i] eq "-variable_length_reads") {
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
            if(!(open(TESTIN, $altgene_file))) {
                print ERRORLOG "\nERROR: cannot open '$altgene_file' for reading.\n\n";
                die "\nERROR: cannot open '$altgene_file' for reading.\n\n";
            }
            close(TESTIN);
	    $optionrecognized = 1;
	}

	if($ARGV[$i] eq "-qualsfile" || $ARGV[$i] eq "-qualfile") {
	    $quals_specified = "true";
            $i++;
            $quals_file = $ARGV[$i];
            $quals = "true";
            if($quals_file =~ /\//) {
               print ERRORLOG "ERROR: do not specify -quals file with a full path, put it in the '$output_dir' directory.\n\n";
               die "ERROR: do not specify -quals file with a full path, put it in the '$output_dir' directory.\n\n";
            }
	    $optionrecognized = 1;
	}

	if($ARGV[$i] eq "-minidentity") {
	    $minidentity = $ARGV[$i+1];
	    $i++;
	    if(!($minidentity =~ /^\d+$/ && $minidentity <= 100)) {
                print ERRORLOG "\nERROR: minidentity must be an integer between zero and 100.\nYou have given '$minidentity'.\n\n";
		die "\nERROR: minidentity must be an integer between zero and 100.\nYou have given '$minidentity'.\n\n";
	    }
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-minlength") {
	    $minlength = $ARGV[$i+1];
	    $i++;
	    if(!($minlength =~ /^\d+$/ && $minlength >= 10)) {
                print ERRORLOG "\nERROR: minlength must be an integer >= 10.\nYou have given '$minlength'.\n\n";
		die "\nERROR: minlength must be an integer >= 10.\nYou have given '$minlength'.\n\n";
	    }
	    $optionrecognized = 1;
	}
	if($ARGV[$i] eq "-limitNU") {
	    $NU_limit = $ARGV[$i+1];
	    $i++;
	    $limitNUhard = "true";
	    if(!($NU_limit =~ /^\d+$/ && $NU_limit > 0)) {
                print ERRORLOG "\nERROR: -limitNU must be an integer greater than zero.\nYou have given '$NU_limit'.\n\n";
		die "\nERROR: -limitNU must be an integer greater than zero.\nYou have given '$NU_limit'.\n\n";
	    }
	    $optionrecognized = 1;
	}
	if($optionrecognized == 0) {
            print ERRORLOG "\nERROR: option $ARGV[$i] not recognized.\n\n";
	    die "\nERROR: option $ARGV[$i] not recognized.\n\n";
	}
    }
}

$H = `hostname`;
if($H =~ /beta.genomics.upenn.edu/ && $qsub eq "false") {
    print ERRORLOG "ERROR: you cannot run RUM on the PGFI cluster without using the -qsub option.\n\n";
    die "ERROR: you cannot run RUM on the PGFI cluster without using the -qsub option.\n\n";
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
	    print "killing $pid\n";
	    `kill -9 $pid`;
	}
    }
    for($i=0; $i<@candidates; $i++) {
	if($candidates[$i] =~ /^\s*(\d+)\s.*(\s|\/)$outdir(\s|\/)/) {
	    if(!($candidates[$i] =~ /pipeline.\d+.sh/)) {
		$pid = $1;
		print "killing $pid\n";
		`kill -9 $pid`;
	    }
	}
    }
    exit();
}


print "

RUM version: $version

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
     print "\nWarning: You are using qsub - so if you have installed RUM somewhere other than your\nhome directory, then you will probably need to specify everything with full paths,\including in the <rum config> file, nor this may not work.\n\n";
}

if($postprocess eq "false") {
     sleep(1);
     print "Please wait while I check that everything is in order.\n\n";
     sleep(1);
     print "This could take a few minutes.\n\n";
     sleep(1);
}

open(INFILE, $configfile);
$gene_annot_file = <INFILE>;
chomp($gene_annot_file);
if($dna eq "false") {
    if(!(-e $gene_annot_file)) {
       print ERRORLOG "\nERROR: the file '$gene_annot_file' does not seem to exist.\n\n";
       die "\nERROR: the file '$gene_annot_file' does not seem to exist.\n\n";
    }
}
$bowtie_exe = <INFILE>;
chomp($bowtie_exe);
if(!(-e $bowtie_exe)) {
    print ERRORLOG "\nERROR: the executable '$bowtie_exe' does not seem to exist.\n\n";
}
$blat_exe = <INFILE>;
chomp($blat_exe);
if(!(-e $blat_exe)) {
    print ERRORLOG "\nERROR: the executable '$blat_exe' does not seem to exist.\n\n";
}
$mdust_exe = <INFILE>;
chomp($mdust_exe);
if(!(-e $mdust_exe)) {
    print ERRORLOG "\nERROR: the executable '$mdust_exe' does not seem to exist.\n\n";
}
$genome_bowtie = <INFILE>;
chomp($genome_bowtie);
$transcriptome_bowtie = <INFILE>;
chomp($transcriptome_bowtie);
$genome_blat = <INFILE>;
chomp($genome_blat);
if(!(-e $genome_blat)) {
    print ERRORLOG "\nERROR: the file '$genome_blat' does not seem to exist.\n\n";
}
$scripts_dir = <INFILE>;
$scripts_dir =~ s!/$!!;
chomp($scripts_dir);
$lib = <INFILE>;
$lib =~ s!/$!!;
chomp($lib);
$genomefa = $genome_blat;
close(INFILE);
$genome_size = -s $genomefa / 1000000000;
$min_ram = int($genome_size * 1.67)+1;

if($postprocess eq "false") {
     if($qsub eq "true") {
         print "You have chosen to submit the jobs using 'qsub'.  I'm going to assume each node has\nsufficient RAM for this.  If you are running a mammalian genome then you should have\nat least 6 Gigs per node.\n\n";
     } else {
          print "I'm going to try to figure out how much RAM you have.\nIf you see some error messages here, don't worry, these are harmless.\n\n";
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
                   print "\nWarning: I could not determine how much RAM you have.  If you have less\nthan $min_ram gigs per chunk this might not work.  I'm going to proceed with fingers crossed.\n\n";
                   $ram = $min_ram;
               } else {
                   $RAMperchunk = int($totalram / $numchunks);
               }
          }
          if($did_not_figure_out_ram eq "false") {
              if($RAMperchunk >= $min_ram) {
                  print "It seems like you have $totalram Gb of RAM on your machine.\n";
                  print "\nUnless you have too much other stuff running, RAM should not be a problem.\n";
              } else {
                  print "\nWarning: you have only $RAMperchunk Gb of RAM per chunk.  Based on the\nsize of your genome you will probably need more like $min_ram Gb per chunk.\nAnyway I will try and see what happens.\n\n";
              }
              $ram = $min_ram;
              if($ram < 6 && $ram < $RAMperchunk) {
                   $ram = $RAMperchunk;
                   if($ram > 6) {
                       $ram = 6;
                   }
              }
              sleep(1);
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
            print ERRORLOG "\nERROR: You seem to already have an instance of RUM_runner.pl running on the\nsame working directory.  This will cause collisions of the temporary files.\n\nExiting.\n\nTry killing by running the same command with -kill.\nIf that doesn't work use kill -9 on the process ID.\n\n";
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
	    print "killing $pid\n";
	    `kill -9 $pid`;
	}
    }
    for($i=0; $i<@candidates; $i++) {
	if($candidates[$i] =~ /^\s*(\d+)\s.*(\s|\/)$outdir(\s|\/)/) {
	    if(!($candidates[$i] =~ /pipeline.\d+.sh/)) {
		$pid = $1;
		print "killing $pid\n";
		`kill -9 $pid`;
	    }
	}
    }
    exit();
}
if($kill eq "false") {
    sleep(1);
    print "\nChecking for phantom processes from prior runs that might need to be killed.\n\n";
    $outdir = $output_dir;
    $str = `ps x | grep $outdir | grep -v RUM_runner.pl`;
    @candidates = split(/\n/,$str);
    $cleanedflag = 0;
    for($i=0; $i<@candidates; $i++) {
	if($candidates[$i] =~ /^\s*(\d+)\s.*(\s|\/)$outdir\/pipeline.\d+.sh/) {
	    $pid = $1;
	    print "killing $pid\n";
	    `kill -9 $pid`;
            $cleanedflag = 1;
	}
    }
    for($i=0; $i<@candidates; $i++) {
	if($candidates[$i] =~ /^\s*(\d+)\s.*(\s|\/)$outdir(\s|\/)/) {
	    if(!($candidates[$i] =~ /pipeline.\d+.sh/)) {
		$pid = $1;
		print "killing $pid\n";
		`kill -9 $pid`;
                $cleanedflag = 1;
	    }
	}
    }
}
if($cleanedflag == 1) {
    sleep(2);
    print "OK there was some cleaning up to do, hopefully that worked.\n\n";
}
sleep(2);

for($i=1; $i<=$numchunks; $i++) {
    $logfile = "$output_dir/rum.log_chunk.$i";
    if (-e $logfile) {
	unlink($logfile);
    }
}


$paired_end = "";
if(($readsfile =~ /,,,/)) {
   $paired_end = "true";
}
$file_needs_splitting = "false";
$preformatted = "false";
if(!($readsfile =~ /,,,/)) {
    if(!(-e $readsfile)) {
        print ERRORLOG "\nERROR: The reads file '$readsfile' does not seem to exist\n\n";
        die "\nERROR: The reads file '$readsfile' does not seem to exist\n\n";
    }
    $head = `head -4 $readsfile`;
    $head =~ /seq.(\d+)(.).*seq.(\d+)(.)/s;
    $num1 = $1;
    $type1 = $2;
    $num2 = $3;
    $type2 = $4;

    $paired_end = "false";
    if($num1 == 1 && $num2 == 1 & $type1 eq 'a' && $type2 eq 'b') {
         $paired_end = "true";
         $file_needs_splitting = "true";
         $preformatted = "true";
    } 
    if($num1 == 1 && $num2 == 2 & $type1 eq 'a' && $type2 eq 'a') {
         $paired_end = "false";
         $file_needs_splitting = "true";
         $preformatted = "true";
    } 
}

if($num_insertions_allowed > 1 && $paired_end eq "true") {
    print ERRORLOG "\nERROR: for paired end data, you cannot set -max_insertions_per_read to be greater than 1.\n\n";
    die "\nERROR: for paired end data, you cannot set -max_insertions_per_read to be greater than 1.\n\n";
}

if($paired_end eq "true") {
     print "Processing as paired-end data\n";
} else {
     print "Processing as single-end data\n";
}

$quals = "false";
if($readsfile =~ /,,,/ && $postprocess eq "false") {
    @a = split(/,,,/, $readsfile);
    if(@a > 2) {
        print ERRORLOG "\nERROR: You've given more than two files separated by three commas, should be at most two files.\n\n";
	die "\nERROR: You've given more than two files separated by three commas, should be at most two files.\n\n";
    }
    if(!(-e $a[0])) {
        print ERRORLOG "\nERROR: The reads file '$a[0]' does not seem to exist\n\n";
	die "\nERROR: The reads file '$a[0]' does not seem to exist\n\n";
    }
    if(!(-e $a[1])) {
        print ERRORLOG "\nERROR: The reads file '$a[1]' does not seem to exist\n\n";
	die "\nERROR: The reads file '$a[1]' does not seem to exist\n\n";
    }
    if($a[0] eq $a[1]) {
        print ERRORLOG "\nERROR: You specified the same file for the forward and reverse reads, must be an error...\n\n";
	die "\nERROR: You specified the same file for the forward and reverse reads, must be an error...\n\n";
    }

    # Going to figure out here if these are standard fastq files

    $head40 = `head -40 $a[0]`;
    $head40 =~ s/^\s*//s;
    $head40 =~ s/\s*$//s;
    @b = split(/\n/, $head40);
    $fastq = "true";
    for($i=0; $i<10; $i++) {
        if(!($b[$i*4] =~ /^@/)) {
            $fastq = "false";
        }
        if(!($b[$i*4+1] =~ /^[acgtnACGTN.]+$/)) {
            $fastq = "false";
        }
        if(!($b[$i*4+2] =~ /^\+/)) {
            $fastq = "false";
        }
    }

   # Check to see if it's fasta

    $fasta = "true";
    for($i=0; $i<10; $i++) {
        if(!($b[$i*2] =~ /^>/)) {
            $fasta = "false";
        }
        if(!($b[$i*2+1] =~ /^[acgtnACGTN.]+$/)) {
            $fasta = "false";
        }
    }

    # Check here that the quality scores are the same length as the reads.

    $FL = `head -50000 $a[0] | wc -l`;
    chomp($FL);
    $FL =~ s/[^\d]//gs;

    `perl $scripts_dir/parse2fasta.pl $a[0] $a[1] | head -$FL > $output_dir/reads_temp.fa 2>> $output_dir/rum.error-log`;
    `perl $scripts_dir/fastq2qualities.pl $a[0] $a[1] | head -$FL > $output_dir/quals_temp.fa 2>> $output_dir/rum.error-log`;
    $X = `head -20 $output_dir/quals_temp.fa`;
    if($X =~ /\S/s && !($X =~ /Sorry, can't figure these files out/s)) {
        open(RFILE, "$output_dir/reads_temp.fa");
        open(QFILE, "$output_dir/quals_temp.fa");
        while($linea = <RFILE>) {
            $lineb = <QFILE>;
            $line1 = <RFILE>;
            $line2 = <QFILE>;
            chomp($line1);
            chomp($line2);
            if(length($line1) != length($line2)) {
               print ERRORLOG "ERROR: It seems your read lengths differ from your quality string lengths.\nCheck line:\n$linea$line1\n$lineb$line2\n\n";
               die "ERROR: It seems your read lengths differ from your quality string lengths.\nCheck line:\n$linea$line1\n$lineb$line2\n\n";
            }
        }
    }

    # Check that reads are not variable length

    if($X =~ /\S/s) {
        open(RFILE, "$output_dir/reads_temp.fa");
        $length_flag = 0;
        while($linea = <RFILE>) {
            $line1 = <RFILE>;
            chomp($line1);
            if($length_flag == 0) {
                 $length_hold = length($line1);
                 $length_flag = 1;
            }
            if(length($line1) != $length_hold && $variable_read_lengths eq 'false') {
               print ERRORLOG "\nWARNING: It seems your read lengths vary, but you didn't set -variable_length_reads.\nI'm going to set it for you, but it's generally safer to set it on the command-line since\nI only spot check the file.\n\n";
               print "\nWARNING: It seems your read lengths vary, but you didn't set -variable_length_reads.\nI'm going to set it for you, but it's generally safer to set it on the command-line since\nI only spot check the file.\n\n";
               $variable_read_lengths = "true";
            }
            $length_hold = length($line1);
        }
    }

    # Clean up:

    unlink("$output_dir/reads_temp.fa");
    unlink("$output_dir/quals_temp.fa");

    # Done checking

    print "\nReformatting reads file... please be patient.\n";

    if($fastq eq "true" && $variable_read_lengths eq "false" && $FL == 50000) {
        `perl $scripts_dir/parsefastq.pl $a[0],,,$a[1] $numchunks $output_dir/reads.fa $output_dir/quals.fa 2>> $output_dir/rum.error-log`;
   	$quals = "true";
        $file_needs_splitting = "false";
    } elsif($fasta eq "true" && $variable_read_lengths eq "false" && $FL == 50000 && $preformatted eq "false") {
        `perl $scripts_dir/parsefasta.pl $a[0],,,$a[1] $numchunks $output_dir/reads.fa 2>> $output_dir/rum.error-log`;
   	$quals = "false";
        $file_needs_splitting = "false";
    } elsif($preformatted eq "false") {
        `perl $scripts_dir/parse2fasta.pl $a[0] $a[1] > $output_dir/reads.fa 2>> $output_dir/rum.error-log`;
        `perl $scripts_dir/fastq2qualities.pl $a[0] $a[1] > $output_dir/quals.fa 2>> $output_dir/rum.error-log`;
        $file_needs_splitting = "true";
        $X = `head -20 $output_dir/quals.fa`;
        if($X =~ /\S/s && !($X =~ /Sorry, can't figure these files out/s)) {
    	     $quals = "true";
        }
    }
    if($preformatted eq "false") {
       $qualsfile = "$output_dir/quals.fa";
       $readsfile = "$output_dir/reads.fa";
    } else {
       $file_needs_splitting = "true";
    }
}

if($postprocess eq "true") {
    $readsfile = "$output_dir/reads.fa";
    if(!(-e $readsfile)) {
        $readsfile = $ARGV[1];
    }
    $qualsfile = "$output_dir/quals.fa";
    $quals = "false";
    if(-e $qualsfile) {
        $X = `head -20 $output_dir/quals.fa`;
        if($X =~ /\S/s && !($X =~ /Sorry, can't figure these files out/s)) {
	    $quals = "true";
        }
    }
}


$head = `head -4 $readsfile`;
$head =~ /seq.(\d+)(.).*seq.(\d+)(.)/s;
$num1 = $1;
$type1 = $2;
$num2 = $3;
$type2 = $4;

if($paired_end eq 'false' && $postprocess eq "false") {
    if($type1 ne "a" || $type2 ne "a") {
        print "\nReformatting reads file... please be patient.\n";

        $head40 = `head -40 $readsfile`;
        $head40 =~ s/^\s*//s;
        $head40 =~ s/\s*$//s;
        @b = split(/\n/, $head40);
        $fastq = "true";
        for($i=0; $i<10; $i++) {
            if(!($b[$i*4] =~ /^@/)) {
                $fastq = "false";
            }
            if(!($b[$i*4+1] =~ /^[acgtnACGTN.]+$/)) {
                $fastq = "false";
            }
            if(!($b[$i*4+2] =~ /^\+/)) {
                $fastq = "false";
            }
        }

        # Check to see if it's fasta

         $fasta = "true";
         for($i=0; $i<10; $i++) {
             if(!($b[$i*2] =~ /^>/)) {
                 $fasta = "false";
             }
             if(!($b[$i*2+1] =~ /^[acgtnACGTN.]+$/)) {
                 $fasta = "false";
             }
         }

        # Check here that the quality scores are the same length as the reads.

        $FL = `head -50000 $readsfile | wc -l`;
        chomp($FL);
        $FL =~ s/[^\d]//gs;

        `perl $scripts_dir/parse2fasta.pl $readsfile | head -$FL > $output_dir/reads_temp.fa 2>> $output_dir/rum.error-log`;
        `perl $scripts_dir/fastq2qualities.pl $readsfile | head -$FL > $output_dir/quals_temp.fa 2>> $output_dir/rum.error-log`;
        $X = `head -20 $output_dir/quals_temp.fa`;
        if($X =~ /\S/s && !($X =~ /Sorry, can't figure these files out/s)) {
            open(RFILE, "$output_dir/reads_temp.fa");
            open(QFILE, "$output_dir/quals_temp.fa");
            while($linea = <RFILE>) {
                $lineb = <QFILE>;
                $line1 = <RFILE>;
                $line2 = <QFILE>;
                chomp($line1);
                chomp($line2);
                if(length($line1) != length($line2)) {
                   print ERRORLOG "ERROR: It seems your read lengths differ from your quality string lengths.\nCheck line:\n$linea$line1\n$lineb$line2\n\n";
                   die "ERROR: It seems your read lengths differ from your quality string lengths.\nCheck line:\n$linea$line1\n$lineb$line2\n\n";
                }
            }
        }
    
        # Check that reads are not variable length
    
        if($X =~ /\S/s) {
            open(RFILE, "$output_dir/reads_temp.fa");
            $length_flag = 0;
            while($linea = <RFILE>) {
                $line1 = <RFILE>;
                chomp($line1);
                if($length_flag == 0) {
                     $length_hold = length($line1);
                     $length_flag = 1;
                }
                if(length($line1) != $length_hold && $variable_read_lengths eq 'false') {
                   print ERRORLOG "\nWARNING: It seems your read lengths vary, but you didn't set -variable_length_reads.\nI'm going to set it for you, but it's generally safer to set it on the command-line since\nI only spot check the file.\n\n";
                   print "\nWARNING: It seems your read lengths vary, but you didn't set -variable_length_reads.\nI'm going to set it for you, but it's generally safer to set it on the command-line since\nI only spot check the file.\n\n";
                   $variable_read_lengths = "true";
                }
                $length_hold = length($line1);
            }
        }

        # Clean up:
        unlink("$output_dir/reads_temp.fa");
        unlink("$output_dir/quals_temp.fa");

        # Done checking

        if($fastq eq "true" && $variable_read_lengths eq "false" && $FL == 50000) {
           `perl $scripts_dir/parsefastq.pl $readsfile $numchunks $output_dir/reads.fa $output_dir/quals.fa 2>> $output_dir/rum.error-log`;
           $quals = "true";
           $file_needs_splitting = "false";
        } elsif($fasta eq "true" && $variable_read_lengths eq "false" && $FL == 50000 && $preformatted eq "false") {
            `perl $scripts_dir/parsefasta.pl $readsfile $numchunks $output_dir/reads.fa 2>> $output_dir/rum.error-log`;
       	    $quals = "false";
            $file_needs_splitting = "false";
        } elsif($preformatted eq "false") {
           `perl $scripts_dir/parse2fasta.pl $readsfile > $output_dir/reads.fa 2>> $output_dir/rum.error-log`;
           `perl $scripts_dir/fastq2qualities.pl $readsfile > $output_dir/quals.fa 2>> $output_dir/rum.error-log`;
           $file_needs_splitting = "true";
           $X = `head -20 $output_dir/quals.fa`;
           if($X =~ /\S/s && !($X =~ /Sorry, can't figure these files out/s)) {
	       $quals = "true";
	   }
        }
        if($preformatted eq "true") {
            $file_needs_splitting = "true";
        } else {
            $readsfile = "$output_dir/reads.fa";
            $qualsfile = "$output_dir/quals.fa";
        }
    }
}

$head = `head -2 $readsfile | tail -1`;
chomp($head);
@a = split(//,$head);
if($variable_read_lengths eq "false") {
   $readlength = @a;
   if($minlength > $readlength) {
       print ERRORLOG "ERROR: you specified a minimum length alignment to report as '$minlength', however\nyour read length is only $readlength\n";
       die "ERROR: you specified a minimum length alignment to report as '$minlength', however\nyour read length is only $readlength\n";
   }
} else {
   $readlength = "v";
}

if($NUMCHUNKS =~ /(\d+)s/) {
    $file_needs_splitting = "false";
}

if($file_needs_splitting eq "true") {
    $x = &breakup_file($readsfile, $numchunks);
    print "Splitting files ... please be patient.\n\n";
    $qualflag = 0;
    if($quals eq "true" || $quals_specified eq "true") {
        print "Half done splitting...\n\n";
        $qualflag = 1;
        if($quals_specified eq 'true') {
       	    $x = &breakup_file("$output_dir/$quals_file", $numchunks);
        } else {
       	    $x = &breakup_file($qualsfile, $numchunks);
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
   if($match_length_cutoff >= .8 * $rl) {
      if($match_length_cutoff == 0) {
         $match_length_cutoff = int(.6 * $rl);
       }
   }
} else {
	$match_length_cutoff = $minlength;
}

if($quals_specified eq 'true') {
    if(!(open(TESTIN, "$output_dir/$quals_file"))) {
       print ERRORLOG "\nERROR: cannot open '$quals_file' for reading, it should be in the '$output_dir' directory.\n\n";
       die "\nERROR: cannot open '$quals_file' for reading, it should be in the '$output_dir' directory.\n\n";
    } 
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

open(LOGFILE, ">$output_dir/rum.log_master");
print LOGFILE "RUM version: $version\n";
print LOGFILE "\nJob ID: $JID\n";
print LOGFILE "\nstart: $date\n";
print LOGFILE "config file: $configfile\n";
print LOGFILE "readsfile: $ARGV[1]\n";
print LOGFILE "output_dir: $output_dir\n";
if($variable_read_lengths eq "false") {
    print LOGFILE "readlength: $rl\n";
} else {
    print LOGFILE "readlength: variable\n";
}
$NR = &format_large_int($nr);
if($paired_end eq 'false') {
    print LOGFILE "number of reads: $NR\n";
} else {
    print LOGFILE "number of read pairs: $NR\n";
}

if($variable_read_lengths eq "false" || $minlength > 0) {
    if($minlength == 0) {
        print LOGFILE "minimum length alignment to report: $match_length_cutoff\n  *** NOTE: If you want shorter alignments reported, use the -minlength option.\n";
    } else {
        print LOGFILE "minimum length alignment to report: $match_length_cutoff.\n";
    }
    print "\n *** NOTE: I am going to report alginments of length $match_length_cutoff.\n";
    print "If you want shorter alignments to be reported, use the -minlength option.\n\n";
} else {
    print LOGFILE "minimum length alignment to report: NA since read length is variable\n";
}
$nc = $numchunks;
$nc =~ s/s//;
print LOGFILE "numchunks: $nc\n";
print LOGFILE "name: $name\n";
print LOGFILE "paired_end: $paired_end\n";
print LOGFILE "fast: $fast\n";
print LOGFILE "ram per chunk: $ram\n";
print LOGFILE "limitBowtieNU: $limitNU\n";
print LOGFILE "limitNU: $limitNUhard\n";
print LOGFILE "dna: $dna\n";
print LOGFILE "qsub: $qsub\n";
print LOGFILE "blat minidentity: $minidentity\n";
print LOGFILE "output junctions: $junctions\n";
print LOGFILE "output quantified values: $quantify\n";
print LOGFILE "strand specific: $strandspecific\n";
print LOGFILE "number insertions allowed per read: $num_insertions_allowed\n";
print LOGFILE "count mismatches: $countmismatches\n";

print ERRORLOG "\nNOTE: I am going to report alginments of length $match_length_cutoff.\n";
print ERRORLOG "  *** If you want shorter alignments to be reported, use the -minlength option.\n\n";

if($readlength ne "v" && $readlength < 55 && $limitNU eq "false") {
    print ERRORLOG "\nWARNING: you have pretty short reads ($readlength bases).  If you have a large\n";
    print ERRORLOG "genome such as mouse or human then the files of ambiguous mappers could grow\n";
    print ERRORLOG "very large, in this case it's recommended to run with the -limitBowtieNU option.  You\n";
    print ERRORLOG "can watch the files that start with 'X' and 'Y' to see if they are growing\n";
    print ERRORLOG "larger than 10 gigabytes per million reads at which point you might want to use.\n";
    print ERRORLOG "-limitNU\n\n";

    print "\n\nWARNING: you have pretty short reads ($readlength bases).  If you have a large\n";
    print "genome such as mouse or human then the files of ambiguous mappers could grow\n";
    print "very large, in this case it's recommended to run with the -limitBowtieNU option.  You\n";
    print "can watch the files that start with 'X' and 'Y' to see if they are growing\n";
    print "larger than 10 gigabytes per million reads at which point you might want to use.\n";
    print "-limitNU\n\n";
}

if($blatonly eq "true" && $dna eq "true") {
    $dna = "false";  # because blat only is dna only anyway, and setting them both breaks things below
}

open(OUT, ">$output_dir/restart.ids");
print OUT "";
close(OUT);

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

    print "Number of Chunks: $numchunks\n";
    print "Reads File: $ARGV[1]\n";
    print "Paired End: $paired_end\n";

    $readsfile =~ s!.*/!!;
    $readsfile = $output_dir . "/" . $readsfile;
    $t = `tail -2 $readsfile`;
    $t =~ /seq.(\d+)/;
    $NumSeqs = $1;
    $f = &format_large_int($NumSeqs);
    if($paired_end eq 'true') {
       print "Number of Read Pairs: $f\n";
    } else {
       print "Number of Reads: $f\n";
    }
    
    print "\nEverything seems okay, I am going to fire off the job.\n\n";
    
    for($i=1; $i<=$numchunks; $i++) {
        if(!(open(EOUT, ">$output_dir/errorlog.$i"))) {
            print ERRORLOG "\nERROR: cannot open '$output_dir/errorlog.$i' for writing\n\n";
            die "\nERROR: cannot open '$output_dir/errorlog.$i' for writing\n\n";
        }
        close(EOUT);
        $pipeline_file = $pipeline_template;
        $pipeline_file =~ s!ERRORFILE!$output_dir/errorlog!gs;
        if($limitNUhard eq "true") {
    	   $pipeline_file =~ s!LIMITNUCUTOFF!$NU_limit!gs;
        } else {
    	   $pipeline_file =~ s!perl SCRIPTSDIR/limit_NU.pl OUTDIR/RUM_NU_temp3.CHUNK LIMITNUCUTOFF > OUTDIR/RUM_NU.CHUNK[^\n]*\n!mv OUTDIR/RUM_NU_temp3.CHUNK OUTDIR/RUM_NU.CHUNK\n!gs;
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
           $pipeline_file =~ s/[^\n]+quant.S2s[^\n]*\n//sg;
           $pipeline_file =~ s/[^\n]+quant.S1a[^\n]*\n//sg;
           $pipeline_file =~ s/[^\n]+quant.S2a[^\n]*\n//sg;
        }
        if($ram != 6) {
           $pipeline_file =~ s!RAM!-ram $ram!gs;
        } else {
           $pipeline_file =~ s! RAM!!gs;
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
        if($dna eq 'true') {
    	   $pipeline_file =~ s!DNA!-dna!gs;
        } else {
    	   $pipeline_file =~ s!DNA!!gs;
        }
        if($minlength > 0) {
    	   $pipeline_file =~ s!MATCHLENGTHCUTOFF!-match_length_cutoff $minlength!gs;
           $pipeline_file =~ s!MINOVERLAP!$minlength!gs;
        } else {
    	   $pipeline_file =~ s!MATCHLENGTHCUTOFF!!gs;
           $pipeline_file =~ s!-minoverlap MINOVERLAP!!gs;
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
        if(!(open(OUTFILE, ">$output_dir/$outfile"))) {
            print ERRORLOG "\nERROR: cannot open '$output_dir/$outfile' for writing\n\n";
            die "\nERROR: cannot open '$output_dir/$outfile' for writing\n\n";
        }
        if($qsub eq "true") {
    	   $pipeline_file =~ s!2>>\s*[^\s]*!!gs;
    	   $pipeline_file =~ s!2>\s*[^\s]*!!gs;
        }
        print OUTFILE $pipeline_file;

        # Add postprocessing steps to the last chunk only:

        if($i == $numchunks) {
            $t = `tail -2 $readsfile`;
            $t =~ /seq.(\d+)/s;
            $NumSeqs = $1;
            $PPlog = "postprocessing_$name" . ".log";
            $shellscript = "\n\n# Postprocessing stuff starts here...\n\n";
            $shellscript = $shellscript . "perl $scripts_dir/wait.pl $output_dir $JID 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
            if($NumSeqs =~ /(\d+)/) {
                $shellscript = $shellscript . "echo 'computing mapping statistics' > $output_dir/$PPlog\n";
                $shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
                $shellscript = $shellscript . "perl $scripts_dir/count_reads_mapped.pl $output_dir/RUM_Unique $output_dir/RUM_NU -minseq 1 -maxseq $NumSeqs > $output_dir/mapping_stats.txt 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
            } else {
                $shellscript = $shellscript . "perl $scripts_dir/count_reads_mapped.pl $output_dir/RUM_Unique $output_dir/RUM_NU -minseq 1 > $output_dir/mapping_stats.txt 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
            }
            if($quantify eq "true") {
                $shellscript = $shellscript . "echo 'merging feature quantifications' >> $output_dir/$PPlog\n";
                $shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
                if($strandspecific eq 'true') {
                    $shellscript = $shellscript . "perl $scripts_dir/merge_quants.pl $output_dir $numchunks $output_dir/feature_quantifications.ps -strand ps -chunk_ids_file $output_dir/restart.ids 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
                    $shellscript = $shellscript . "perl $scripts_dir/merge_quants.pl $output_dir $numchunks $output_dir/feature_quantifications.ms -strand ms -chunk_ids_file $output_dir/restart.ids 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
                    $shellscript = $shellscript . "perl $scripts_dir/merge_quants.pl $output_dir $numchunks $output_dir/feature_quantifications.pa -strand pa -chunk_ids_file $output_dir/restart.ids 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
                    $shellscript = $shellscript . "perl $scripts_dir/merge_quants.pl $output_dir $numchunks $output_dir/feature_quantifications.ma -strand ma -chunk_ids_file $output_dir/restart.ids 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
                    $shellscript = $shellscript . "perl $scripts_dir/merge_quants_strandspecific.pl $output_dir/feature_quantifications.ps $output_dir/feature_quantifications.ms $output_dir/feature_quantifications.pa $output_dir/feature_quantifications.ma $gene_annot_file $output_dir/feature_quantifications_$name 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
                } else {
                    $shellscript = $shellscript . "perl $scripts_dir/merge_quants.pl $output_dir $numchunks $output_dir/feature_quantifications_$name -chunk_ids_file $output_dir/restart.ids 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
                }
            }

            $string = "$output_dir/RUM_Unique.sorted";
            for($j=1; $j<$numchunks+1; $j++) {
                $string = $string . " $output_dir/RUM_Unique.sorted.$j";
            }
            $shellscript = $shellscript . "echo 'merging RUM_Unique.sorted files' > $output_dir/$PPlog\n";
            $shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
            $shellscript = $shellscript . "perl $scripts_dir/merge_sorted_RUM_files.pl $string -chunk_ids_file $output_dir/restart.ids 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
            $string = "$output_dir/RUM_NU.sorted";
            for($j=1; $j<$numchunks+1; $j++) {
                $string = $string . " $output_dir/RUM_NU.sorted.$j";
            }
            $shellscript = $shellscript . "echo 'merging RUM_NU.sorted files' >> $output_dir/$PPlog\n";
            $shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
            $shellscript = $shellscript . "perl $scripts_dir/merge_sorted_RUM_files.pl $string -chunk_ids_file $output_dir/restart.ids 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
            
            $string = "$output_dir/mapping_stats.txt";
            for($j=1; $j<$numchunks+1; $j++) {
                $string = $string . " $output_dir/chr_counts_u.$j";
            }
            $shellscript = $shellscript . "echo '' >> $output_dir/mapping_stats.txt\n";
            $shellscript = $shellscript . "echo 'RUM_Unique reads per chromosome' >> $output_dir/mapping_stats.txt\n";
            $shellscript = $shellscript . "perl $scripts_dir/merge_chr_counts.pl $string -chunk_ids_file $output_dir/restart.ids 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
            
            $string = "$output_dir/mapping_stats.txt";
            for($j=1; $j<$numchunks+1; $j++) {
                $string = $string . " $output_dir/chr_counts_nu.$j";
            }
            $shellscript = $shellscript . "echo '' >> $output_dir/mapping_stats.txt\n";
            $shellscript = $shellscript . "echo 'RUM_NU reads per chromosome' >> $output_dir/mapping_stats.txt\n";
            $shellscript = $shellscript . "perl $scripts_dir/merge_chr_counts.pl $string -chunk_ids_file $output_dir/restart.ids 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
            
            if($junctions eq "true") {
               $shellscript = $shellscript . "echo 'computing junctions' >> $output_dir/$PPlog\n";
               $shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
               if($altgenes eq "true") {
                   $shellscript = $shellscript . "perl $scripts_dir/make_RUM_junctions_file.pl $output_dir/RUM_Unique $output_dir/RUM_NU $genomefa $altgene_file $output_dir/junctions_all.rum $output_dir/junctions_all.bed $output_dir/junctions_high-quality.bed -faok 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
               } else {
                   $shellscript = $shellscript . "perl $scripts_dir/make_RUM_junctions_file.pl $output_dir/RUM_Unique $output_dir/RUM_NU $genomefa $gene_annot_file $output_dir/junctions_all.rum $output_dir/junctions_all.bed $output_dir/junctions_high-quality.bed -faok 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
               }
            }
            
            $shellscript = $shellscript . "echo 'making coverage plots' >> $output_dir/$PPlog\n";
            $shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
            $shellscript = $shellscript . "perl $scripts_dir/rum2cov.pl $output_dir/RUM_Unique.sorted $output_dir/RUM_Unique.cov -name \"$name Unique Mappers\" 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
            $shellscript = $shellscript . "perl $scripts_dir/rum2cov.pl $output_dir/RUM_NU.sorted $output_dir/RUM_NU.cov -name \"$name Non-Unique Mappers\" 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
            if($strandspecific eq 'true') {
                  # breakup RUM_Unique and RUM_NU files into plus and minus
                  $shellscript = $shellscript . "perl $scripts_dir/breakup_RUM_files_by_strand.pl $output_dir/RUM_Unique.sorted $output_dir/RUM_Unique.sorted.plus $output_dir/RUM_Unique.sorted.minus 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
                  $shellscript = $shellscript . "perl $scripts_dir/breakup_RUM_files_by_strand.pl $output_dir/RUM_NU.sorted $output_dir/RUM_NU.sorted.plus $output_dir/RUM_NU.sorted.minus 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
                  # run rum2cov on all four files
                  $shellscript = $shellscript . "perl $scripts_dir/rum2cov.pl $output_dir/RUM_Unique.sorted.plus $output_dir/RUM_Unique.plus.cov -name \"$name Unique Mappers Plus Strand\" 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
                  $shellscript = $shellscript . "perl $scripts_dir/rum2cov.pl $output_dir/RUM_Unique.sorted.minus $output_dir/RUM_Unique.minus.cov -name \"$name Unique Mappers Minus Strand\" 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
                  $shellscript = $shellscript . "perl $scripts_dir/rum2cov.pl $output_dir/RUM_NU.sorted.plus $output_dir/RUM_NU.plus.cov -name \"$name Non-Unique Mappers Plus Strand\" 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
                  $shellscript = $shellscript . "perl $scripts_dir/rum2cov.pl $output_dir/RUM_NU.sorted.minus $output_dir/RUM_NU.minus.cov -name \"$name Non-Unique Mappers Minus Strand\" 2>> $output_dir/PostProcessing-errorlog || exit 1\n";
            }
            $shellscript = $shellscript . "echo finished >> $output_dir/$PPlog\n";
            $shellscript = $shellscript . "echo `date` >> $output_dir/$PPlog\n";
            if($qsub eq "true") {
        	   $shellscript =~ s!2>>\s*[^\s]*!!gs;
        	   $shellscript =~ s!2>\s*[^\s]*!!gs;
            }
            print OUTFILE $shellscript;
        }
        close(OUTFILE);
    
        if($qsub eq "true") {
           $ofile = $output_dir . "/chunk.$i" . ".o";
           $efile = $output_dir . "/errorlog.$i";
           $MEM = $ram . "G";
    	   $Q = `qsub -l mem_free=$MEM -o $ofile -e $efile $output_dir/$outfile`;
           $Q =~ /Your job (\d+)/;
           $jobid{$i} = $1;
        } else {
    	   system("/bin/bash $output_dir/$outfile &");
        }
        print "Chunk $i initiated\n";
        $status{$i} = 1;
    }
    if($numchunks > 1) {
        print "\nAll chunks initiated, now the long wait...\n";
        print "\nI'm going to watch for all chunks to finish, then I will merge everything...\n\n";
        sleep(2);
    } else {
        print "\nThe job has been initiated, now the long wait...\n";
        sleep(2);
    }
    
    $currenttime = time();
    $lastannouncetime = $currenttime;
    $numannouncements = 0;
    $doneflag = 0;
    
    while($doneflag == 0) {
        sleep(30);
        $doneflag = 1;
        $numdone = 0;
        for($i=1; $i<=$numchunks; $i++) {
            if($restarted{$i} =~ /\S/) {
                $logfile = "$output_dir/rum.log_chunk.$i.$restarted{$i}";
            } else {
                $logfile = "$output_dir/rum.log_chunk.$i";
            }
            if (-e $logfile) {
    	        $x = `cat $logfile`;
    	        if(!($x =~ /pipeline complete/s)) {
    		    $doneflag = 0;
    	        } else {
    		    $numdone++;
    		    if($status{$i} == 1) {
    		        $status{$i} = 2;
		        print "\n *** Chunk $i has finished.\n";
                        if($i != $numchunks) {
                            delete $jobid{$i};
                        }
		    }
	        }
	    }
    	    else {
	        $doneflag = 0;
	    }
        }
        if($doneflag == 0) {
	    $currenttime = time();
    	    if($currenttime - $lastannouncetime > 3600) {
	        $lastannouncetime = $currenttime;
	        $numannouncements++;
	        if($numannouncements == 1) {
		    if($numdone == 1) {
		        print "\nIt has been $numannouncements hour, $numdone chunk has finished.\n";
		    } else {
		        print "\nIt has been $numannouncements hour, $numdone chunks have finished.\n";
		    }
	        } else {
		    if($numdone == 1) {
		        print "\nIt has been $numannouncements hours, $numdone chunk has finished.\n";
		    } else {
		        print "\nIt has been $numannouncements hours, $numdone chunks have finished.\n";
		    }
	        }
    	    }
        }
        for($i=1; $i<=$numchunks; $i++) {
          # Check here to make sure node still running
              if($qsub eq 'true') {
                  if($restarted{$i} =~ /\S/) {
                      $logfile = "$output_dir/rum.log_chunk.$i.$restarted{$i}";
                  } else {
                      $logfile = "$output_dir/rum.log_chunk.$i";
                  }
                  $x = "";
                  if (-e $logfile) {
        	        $x = `cat $logfile`;
                  }
                  $Jobid = $jobid{$i};
                  $X = `qstat -j $Jobid | grep job_number 2> $output_dir/temp.12321`;
                  unlink("$output_dir/temp.12321");
                  if(!($X =~ /job_number:\s+$Jobid/s) && (!($x =~ /pipeline complete/s) || ($x =~ /pipeline complete/s && $i == $numchunks && $status{$i} == 2))) {
                       $DATE = `date`;
                       $DATE =~ s/^\s+//;
                       $DATE =~ s/\s+$//;
                       print ERRORLOG "\n *** Chunk $i seems to have failed sometime around $DATE!  Trying to restart it...\n";
                       print "\n *** Chunk $i seems to have failed sometime around $DATE!\nDon't panic, I'm going to try to restart it.\n";
                       $ofile = $output_dir . "/chunk.restart.$i" . ".o";
                       $efile = $output_dir . "/chunk.restart.$i" . ".e";
                       $outfile = "pipeline." . $i . ".sh";
                       $FILE = `cat $output_dir/$outfile`;
                       $restarted{$i}++;
                       open(OUT, ">$output_dir/restart.ids");
                       foreach $key (keys %restarted) {
                           print OUT "$key\t$restarted{$key}\n";
                       }
                       close(OUT);
                       # changing the names of the files of this chunk to avoid possible collision with
                       # phantom processes that didn't die properly..
                       if($restarted{$i} == 1) {
                           $FILE =~ s/\.$i/.$i.1/g;
                           $J1 = $i;
                           $J3 = $i;
                           `mv $output_dir/reads.fa.$i $output_dir/reads.fa.$i.1`;
                           if(-e "$output_dir/quals.fa.$i") {
                               `mv $output_dir/quals.fa.$i $output_dir/quals.fa.$i.1`;
                           }
                       } else {
                           $J1 = $restarted{$i} - 1;
                           $J2 = $restarted{$i};
                           $FILE =~ s/\.$i\.$J1/.$i.$J2/g;
                           $J3 = "$i.$J1";
                           `mv $output_dir/reads.fa.$i.$J1 $output_dir/reads.fa.$i.$J2`;
                           if(-e "$output_dir/quals.fa.$i.$J1") {
                              `mv $output_dir/quals.fa.$i.$J1 $output_dir/quals.fa.$i.$J2`;
                           }
                       }
                       open(OUTX, ">$output_dir/$outfile");
                       print OUTX $FILE;
                       close(OUTX);

# Note, can't modify the postprocessing chunk to reflect the new file names, since it
# has already been submmitted.  Instead the postprocessing scripts that need file names
# will recover the correct ones from the restart.ids file

                       # remove the old files...  
                       open(OUT, ">>$output_dir/restart_deleted_logs");
                       print OUT "------ chunk $i restarted, here is its error log before it was deleted --------\n";
                       close(OUT);
                       `cat $output_dir/errorlog.$i >> $output_dir/restart_deleted_logs`;

                       &deletefiles($output_dir, $J3);

                       if($i == $numchunks) {
                           # this is the post-processing node.  Check if it finished up to the
                           # post-processing, if so then remove that part so as not to repeat it.
                           if($status{$i} == 2) {  # it has finished
                               $FILE =~ s/# xxx0.*Postprocessing stuff starts here.../\n/s;
                               open(OUTFILE, ">$output_dir/$outfile");
                               print OUTFILE $FILE;
                               close(OUTFILE);
                           }
                       }
                       $MEM = $ram . "G";
                       $Q = `qsub -l mem_free=$MEM -o $ofile -e $efile $output_dir/$outfile`;
                       $Q =~ /Your job (\d+)/;
                       $jobid{$i} = $1;
                       if($jobid{$i} =~ /^\d+$/) {
                             $DATE = `date`;
                             $DATE =~ s/^\s+//;
                             $DATE =~ s/\s+$//;
                             sleep(2);
                             print ERRORLOG " *** Chunk $i seems to have restarted successfully at $DATE.\n\n";
                             print " *** OK chunk $i seems to have restarted.\n\n";
                       } else {
                             print ERRORLOG " *** Hmph, that didn't seem to work.  I'm going to try again in 30 seconds.\nIf this keeps happening then something bigger might be wrong.  If you\ncan't figure it out, write ggrant@pcbi.upenn.edu and let him know.\n\n";
                             print " *** Hmph, that didn't seem to work.  I'm going to try again in 30 seconds.\nIf this keeps happening then something bigger might be wrong.  If you\ncan't figure it out, write ggrant@pcbi.upenn.edu and let him know.\n\n";
                       }
                  }
             } else {
                  # still need to implement this for the non qsub case
             }
         }
    }
}

if($postprocess eq "false") {
     print "\nWhew, all chunks have finished.\n\nNext I will merge everything, create the coverage plots and\ncalculate the quantified values, etc.  This could take some time...\n\n";
} else {
     print "\nOK, will now merge everything, create the coverage plots and\ncalculate the quantified values, etc.  This could take some time...\n\n";
}

if($qsub eq "true") {
    $efile = $output_dir . "/errorlog.$numchunks";
    open(EFILE, ">>$efile");
    print EFILE "\nPost-Processing Log Starts Here\n";
    close(EFILE);
}

if($nocat eq "false") {
    $date = `date`;
    if(defined $restarted{1}) {
        $R = $restarted{1};
        $x = `cp $output_dir/RUM_Unique.1.$R $output_dir/RUM_Unique`;
        $x = `cp $output_dir/RUM_NU.1.$R $output_dir/RUM_NU`;
    } else {
        $x = `cp $output_dir/RUM_Unique.1 $output_dir/RUM_Unique`;
        $x = `cp $output_dir/RUM_NU.1 $output_dir/RUM_NU`;
    }
    for($i=2; $i<=$numchunks; $i++) {
        if(defined $restarted{$i}) {
            $R = $restarted{$i};
            $x = `cat $output_dir/RUM_Unique.$i.$R >> $output_dir/RUM_Unique`;
            $x = `cat $output_dir/RUM_NU.$i.$R >> $output_dir/RUM_NU`;
        } else {
            $x = `cat $output_dir/RUM_Unique.$i >> $output_dir/RUM_Unique`;
            $x = `cat $output_dir/RUM_NU.$i >> $output_dir/RUM_NU`;
        }
    }
    for($i=1; $i<=$numchunks; $i++) {
       if(defined $restarted{$i}) {
           $R = $restarted{$i};
           if(!(open(SAMHEADER, "$output_dir/sam_header.$i.$R"))) {
              print ERRORLOG "\nERROR: Cannot open '$output_dir/sam_header.$i.$R' for reading.\n\n";
              die "\nERROR: Cannot open '$output_dir/sam_header.$i.$R' for reading.\n\n";
           }
       } else {
           if(!(open(SAMHEADER, "$output_dir/sam_header.$i"))) {
              print ERRORLOG "\nERROR: Cannot open '$output_dir/sam_header.$i' for reading.\n\n";
              die "\nERROR: Cannot open '$output_dir/sam_header.$i' for reading.\n\n";
           }
       }
       while($line = <SAMHEADER>) {
           chomp($line);
           $line =~ /SN:([^\s]+)\s/;
           $samheader{$1}=$line;
       }
       close(SAMHEADER);
    }
    if(!(open(SAMOUT, ">$output_dir/RUM.sam"))) {
        print ERRORLOG "\nERROR: Cannot open '$output_dir/RUM.sam' for writing.\n\n";
        die "\nERROR: Cannot open '$output_dir/RUM.sam' for writing.\n\n";
    }
    foreach $key (sort {cmpChrs($a,$b)} keys %samheader) {
        $shout = $samheader{$key};
        print SAMOUT "$shout\n";
    }
    close(SAMOUT);
    for($i=1; $i<=$numchunks; $i++) {
       if(defined $restarted{$i}) {
           $R = $restarted{$i};
           $x = `cat $output_dir/RUM.sam.$i.$R >> $output_dir/RUM.sam`;
       } else {
           $x = `cat $output_dir/RUM.sam.$i >> $output_dir/RUM.sam`;
       }
    }
}

print "Finished creating RUM_Unique/RUM_NU/RUM.sam: $date\n";

if($cleanup eq 'true') {
   print "\nCleaning up some temp files...\n\n";
   `yes|rm $output_dir/RUM.sam.* $output_dir/sam_header.*`;
}

# XXX Need to make a separate shell script now to handle the option -postprocess
#   - not yet implemented..

print "\nStarted postprocessing at $date\n";
print LOGFILE "\nStarted postprocessing at $date\n";

# Write file that wait.pl is watching for, in the shell script for the last chunk.
# Once that is written, wait.pl finishes and the postprocessing can start.

if(!(open(OUTFILE, ">$output_dir/$JID"))) {
    print ERRORLOG "\nERROR: Cannot open '$output_dir/$JID' for writing.\n\n";
    die "\nERROR: Cannot open '$output_dir/$JID' for writing.\n\n";
}
print OUTFILE "$JID\n";
close(OUTFILE);

print "\nWorking, now another wait...\n\n";

$doneflag = 0;

while($doneflag == 0) {
    $doneflag = 1;
    $x = "";
    if (-e "$output_dir/$PPlog") {
	$x = `cat $output_dir/$PPlog`;
	if(!($x =~ /finished/s)) {
	    $doneflag = 0;
	}
    } else {
	$doneflag = 0;
    }

# check here for node failure, and restart if necessary


    if($qsub eq 'true' && $doneflag == 0) {
        $Jobid = $jobid{$numchunks};
        $X = `qstat -j $Jobid | grep job_number`;
        if(!($X =~ /job_number:\s+$Jobid/s)) {
            $DATE = `date`;
            $DATE =~ s/^\s+//;
            $DATE =~ s/\s+$//;
            print ERRORLOG "\n *** Chunk $numchunks seems to have failed during post-processing, sometime around $DATE!\nI'm going to try to restart it.\n";
            print "\n *** Chunk $numchunks seems to have failed during post-processing, sometime around $DATE!\nDon't panic, I'm going to try to restart it.\n";
            $ofile = $output_dir . "/chunk.restart.$numchunks" . ".o";
            $efile = $output_dir . "/errorlog.restart.$numchunks";
            $outfile = "pipeline." . $numchunks . ".sh";

            # first remove the pre-post-processing stuff from the shell script
            $FILE = `cat $output_dir/$outfile`;
            $FILE =~ s/# xxx0.*Postprocessing stuff starts here.../\n/s;
            $FILE =~ s/perl scripts.wait.pl [^\s]+ \d+[^\n]*//s;
            open(OUTFILE, ">$output_dir/$outfile");
            print OUTFILE $FILE;
            close(OUTFILE);

            $MEM = $ram . "G";
            $Q = `qsub -l mem_free=$MEM -o $ofile -e $efile $output_dir/$outfile`;
            $Q =~ /Your job (\d+)/;
            $jobid{$numchunks} = $1;
            if($jobid{$numchunks} =~ /^\d+$/) {
                  sleep(2);
                  print ERRORLOG " *** OK, post-processing seems to have restarted.\n\n";
                  print " *** OK, post-processing seems to have restarted.\n\n";
            } else {
                  print ERRORLOG " *** Hmph, that didn't seem to work.  I'm going to try again in 30 seconds.\nIf this keeps happening then something bigger might be wrong.  If you\ncan't figure it out, write ggrant@pcbi.upenn.edu and let him know.\n\n";
                  print " *** Hmph, that didn't seem to work.  I'm going to try again in 30 seconds.\nIf this keeps happening then something bigger might be wrong.  If you\ncan't figure it out, write ggrant@pcbi.upenn.edu and let him know.\n\n";
            }
       }
   }
   if($doneflag == 0) {
      sleep(30);
   }
}

# Check RUM_Unique and RUM_Unique.sorted are the same size
$filesize1 = -s "$output_dir/RUM_Unique";
$filesize2 = -s "$output_dir/RUM_Unique.sorted";
if($filesize1 != $filesize2) {
    print ERRORLOG "ERROR: RUM_Unique and RUM_Unique.sorted are not the same size.  This probably indicates a problem.\n";
    print "ERROR: RUM_Unique and RUM_Unique.sorted are not the same size.  This probably indicates a problem.\n";
}

# Check RUM_NU and RUM_NU.sorted are the same size
$filesize1 = -s "$output_dir/RUM_NU";
$filesize2 = -s "$output_dir/RUM_NU.sorted";
if($filesize1 != $filesize2) {
    print ERRORLOG "ERROR: RUM_NU and RUM_NU.sorted are not the same size.  This could indicates a problem.\n";
    print "ERROR: RUM_NU and RUM_NU.sorted are not the same size.  This could indicates a problem.\n";
}

# XXX   More error checks to implement:
#
# Find last chr in RUM_Unique and RUM_NU.
# Make sure thr right one of those last chrs is the last chr in RUM.sam.
# Make sure that last chr is the last chr in RUM_Unique.cov and RUM_NU.cov.
# If any of these fail, report them to ERRORLOG.

$check_if_any_errors_already_reported = `grep -i error $output_dir/rum.error-log`;
if(!($check_if_any_errors_already_reported =~ /\S/)) {
   $noerrors = "true";
} else {
   $noerrors = "false";
}

if($qsub eq "false") {
    $E = `cat $output_dir/PostProcessing-errorlog`;
    $E =~ s/^\s*//s;
    $E =~ s/\s*$//s;
} else {
    if(-e "$output_dir/errorlog.restart.$numchunks") {
        $E1 = `cat $output_dir/errorlog.$numchunks`;
        $E1 =~ s/^.*Post-Processing Log Starts Here//s;
        $E = `cat $output_dir/errorlog.restart.$numchunks`;
        $E =~ s/^.*Post-Processing Log Starts Here//s;
        $E = $E1 . "\n" . $E;
    } else {
        $E = `cat $output_dir/errorlog.$numchunks`;
        $E =~ s/^.*Post-Processing Log Starts Here//s;
    }
}
if($E =~ /\S/) {
    print ERRORLOG "\n------- Post Processing Errors -------\n";
    print ERRORLOG "$E\n";
    $noerrors = "false";
}
for($i=1; $i<=$numchunks; $i++) {
    $E = `cat $output_dir/errorlog.$i`;
    $E =~ s/# reads[^\n]+\n//sg;
    $E =~ s/Reported \d+ [^\n]+\n//sg;
    $E =~ s/^\s*//s;
    $E =~ s/\s*$//s;
    if($qsub eq "true") {
        $E =~ s/Post-Processing Log Starts Here.*$//s;
        if(-e "$output_dir/errorlog.restart.$i") {
            $E1 = `cat $output_dir/errorlog.restart.$i`;
            $E1 =~ s/Post-Processing Log Starts Here.*$//s;
            $E1 =~ s/# reads[^\n]+\n//sg;
            $E1 =~ s/Reported \d+ [^\n]+\n//sg;
            $E = $E1 . "\n" . $E;
        }
    }
    if($E =~ /\S/) {
        print ERRORLOG "\n------- errors from chunk $i -------\n";
        print ERRORLOG "$E\n";
        `yes|rm $output_dir/errorlog.$i`;
        $noerrors = "false";
    }
    if(defined $restarted{$i}) {
        $R = $restarted{$i};
        $E = `grep \"$output_dir\" $output_dir/rum.log_chunk.$i.$R | grep -v finished`;
    } else {
        $E = `grep \"$output_dir\" $output_dir/rum.log_chunk.$i | grep -v finished`;
    }
    $E =~ s/^\s*//s;
    $E =~ s/\s*$//s;
    @a = split(/\n/,$E);
    $flag = 0;
    for($j=0; $j<@a; $j++) {
        @b = split(/\s+/,$a[$j]);
        if($b[4] == 0) {
            $file = $b[@b-1];
            if($flag == 0) {
                print ERRORLOG "\n";
                $flag = 1;
            }
            print ERRORLOG "WARNING: temp file '$file' had size zero.\n  *  Could be no mappers in that step, but this often indicates an error.\n";
            $noerrors = "false";
        }
    }
}
if($noerrors eq "true") {
    print ERRORLOG "\nNo Errors. Very good!\n\n";
}
print ERRORLOG "--------------------------------------\n";

if($cleanup eq 'true') {
   print "\nCleaning up some more temp files...\n\n";
   for($i=1; $i<=$numchunks; $i++) {
      if(defined $restarted{$i}) {
         $ext = ".$restarted{$i}";
      } else {
         $ext = "";
      }
      `yes|rm $output_dir/RUM_Unique.$i$ext $output_dir/RUM_NU.$i$ext`;
      `yes|rm $output_dir/RUM_Unique.sorted.$i$ext $output_dir/RUM_NU.sorted.$i$ext`;
      `yes|rm $output_dir/reads.fa.$i$ext`;
      `yes|rm $output_dir/quals.fa.$i$ext`;
   }
   `yes|rm $output_dir/chr_counts*`;
   `yes|rm $output_dir/quant.*`;
   `yes|rm $output_dir/pipeline.*`;
   if($strandspecific eq 'true') {
      `yes|rm $output_dir/feature_quantifications.ps`;
      `yes|rm $output_dir/feature_quantifications.ms`;
      `yes|rm $output_dir/feature_quantifications.pa`;
      `yes|rm $output_dir/feature_quantifications.ma`;
   }
   `yes|rm $output_dir/$JID`;
}

print "\nOkay, all finished.\n\n";

$date = `date`;
print LOGFILE "pipeline finished: $date\n";
close(LOGFILE);
close(ERRORLOG);

sub breakup_file () {
    ($FILE, $numpieces) = @_;

    if(!(open(INFILE, $FILE))) {
       print ERRORLOG "\nERROR: Cannot open '$FILE' for reading.\n\n";
       die "\nERROR: Cannot open '$FILE' for reading.\n\n";
    }
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

sub deletefiles () {
    ($dir, $suffix) = @_;

    $dir =~ s!/$!!;
    $suffix =~ s/^\.//;

    $file[0] = "BlatNU.XXX";
    $file[1] = "BlatUnique.XXX";
    $file[2] = "BowtieNU.XXX";
    $file[3] = "BowtieUnique.XXX";
    $file[4] = "chr_counts_nu.XXX";
    $file[5] = "chr_counts_u.XXX";
    $file[6] = "CNU.XXX";
    $file[7] = "GNU.XXX";
    $file[8] = "GU.XXX";
    $file[9] = "quant.XXX";
    $file[10] = "R.XXX";
    $file[11] = "R.mdust.XXX";
    $file[12] = "rum.log_chunk.XXX";
    $file[13] = "RUM_NU.XXX";
    $file[14] = "RUM_NU_idsorted.XXX";
    $file[15] = "RUM_NU.sorted.XXX";
    $file[16] = "RUM_NU_temp.XXX";
    $file[17] = "RUM_NU_temp2.1";
    $file[18] = "RUM_NU_temp2.XXX";
    $file[19] = "RUM_NU_temp2.3";
    $file[20] = "RUM_NU_temp2.4";
    $file[21] = "RUM_NU_temp3.XXX";
    $file[22] = "RUM.sam.XXX";
    $file[23] = "R.XXX.blat";
    $file[24] = "RUM_Unique.XXX";
    $file[25] = "RUM_Unique.sorted.XXX";
    $file[26] = "RUM_Unique_temp.XXX";
    $file[27] = "RUM_Unique_temp2.XXX";
    $file[28] = "sam_header.XXX";
    $file[29] = "TNU.XXX";
    $file[30] = "TU.XXX";
    $file[31] = "X.XXX";
    $file[32] = "Y.XXX";

    for($i_d=0; $i_d<@file; $i_d++) {
	$F = $dir . "/" . $file[$i_d];
	$F =~ s/XXX/$suffix/;
	if(-e $F) {
	    `yes|rm $F`;
	}
    }
    return "";
}
