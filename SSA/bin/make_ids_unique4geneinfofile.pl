#!/usr/bin/perl

# Written by Gregory R. Grant
# University of Pennsylvania, 2010

if(@ARGV < 2) {
    die "
Usage: make_ids_unique4geneinfofile.pl <gene info file> <output file>

This script is part of the pipeline of scripts used to create RUM indexes.
You should probably not be running it alone.  See the library file:
'how2setup_genome-indexes_forPipeline.txt'.

";
}

open(INFILE, $ARGV[0]);
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    @b = split(/::::/,$a[7]);
    for($i=0; $i<@b; $i++) {
        $b[$i] =~ /(.*)\(([^\)]+)\)$/;
        $id = $1;
        $type = $2;
        $id =~ s/.*://;
        $id =~ s/\(.*//;
        $idcount{$type}{$id}++;
    }
}
close(INFILE);

open(INFILE, $ARGV[0]);
open(OUTFILE, ">$ARGV[1]");
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    @b = split(/::::/,$a[7]);
    for($i=0; $i<@b; $i++) {
        $b[$i] =~ /(.*)\(([^\)]+)\)$/;
        $id = $1;
        $type = $2;
        $id =~ s/.*://;
        $id =~ s/\(.*//;
        $idcount2{$type}{$id}++;
        if($idcount{$type}{$id} > 1) {
            $j = $idcount2{$type}{$id};
            $id_with_number = $id . "[[$j]]";
            $line =~ s/$id\($type\)/$id_with_number\($type\)/;      
        }
    }
    print OUTFILE "$line\n";
}
close(INFILE);
close(OUTFILE);
