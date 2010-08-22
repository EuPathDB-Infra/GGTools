#!/usr/bin/perl

# Written by Gregory R. Grant
# Universiity of Pennsylvania, 2010

if(@ARGV < 2) {
    die "
Usage: sort_RUM.pl <RUM file> <sorted outfile>

This script sorts a RUM output file by sequence number.  It keeps
consistent pairs together.
";
}

$|=1;
open(INFILE, $ARGV[0]);
my @tmp;
while(<INFILE>){
  if(/seq.(\d+)/){
    push(@tmp,[$_,$1]);
  }
}
close INFILE;

print STDERR "Read in ".scalar(@tmp)." alignments\n";

open(OUT,">$ARGV[1]");
foreach my $s (sort{$a->[1] <=> $b->[1]}@tmp){
  print OUT $s->[0];
}

close OUT;

