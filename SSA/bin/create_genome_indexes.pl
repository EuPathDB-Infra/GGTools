#!/usr/bin/perl

# Written by Gregory R. Grant
# Universiry of Pennsylvania, 2010

if(@ARGV < 1) {
    die "
Usage: create_genome_indexes.pl <genome file>

This script is part of the pipeline of scripts used to create RUM indexes.
For more information see the library file: 'how2setup_genome-indexes_forPipeline.txt'.

Genome fasta file must be formatted as described in:
'how2setup_genome-indexes_forPipeline.txt'.

";
}

$infile = $ARGV[0];
$F1 = $infile;
$F1 =~ s/.txt$/.fa/;
$F2 = $infile;
$F2 =~ s/.txt$/_one-line-seqs_temp.fa/;
$F3 = $infile;
$F3 =~ s/.txt$/_one-line-seqs.fa/;

`perl modify_fasta_header_for_genome_seq_database.pl $infile > $F1`;
`perl modify_fa_to_have_seq_on_one_line.pl $F1 > $F2`;
`perl sort_genome_fa_by_chr.pl $F2 >  $F3`;

unlink($F1);
unlink($F2);
