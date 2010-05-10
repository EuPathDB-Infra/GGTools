#!/usr/bin/perl

# Written by Gregory R. Grant
# Universiry of Pennsylvania, 2010

if(@ARGV < 1) {
    die "
Usage: create_gene_indexes.pl <name>

This script is part of the pipeline of scripts used to create RUM indexes.
For more information see the library file: 'how2setup_genome-indexes_forPipeline.txt'.

";
}


$NAME = $ARGV[0];

$N1 = $NAME . "_gene_info_orig.txt";
$N2 = $NAME . "_genome.fa";
$N3 = $NAME . "_genes.fa";
$N4 = $NAME . "_gene_info.txt";

`perl make_master_file_of_genes.pl gene_info_files > gene_info_merged_unsorted.txt`;
`perl fix_geneinfofile_for_neg_introns.pl gene_info_merged_unsorted.txt > gene_info_merged_unsorted_fixed.txt`;
`perl sort_geneinfofile.pl gene_info_merged_unsorted_fixed.txt > gene_info_merged_sorted_fixed.txt`;
`perl make_ids_unique4geneinfofile.pl gene_info_merged_sorted_fixed.txt $N1;`;
`perl get_master_list_of_exons_from_geneinfofile.pl $N1`;
`perl modify_fa_to_have_seq_on_one_line.pl $N2 > temp.fa`;
`perl make_fasta_files_for_master_list_of_genes.pl temp.fa master_list_of_exons.txt $N1 $N4 > $N3`;
