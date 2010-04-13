#!/usr/bin/perl

use strict;
use File::Basename;

&usage unless scalar(@ARGV >= 2);

my $type = shift @ARGV;

&usage unless ($type eq 'single' || $type eq 'paired');

# before we get started, make sure all files exist
foreach my $file (@ARGV) {
  die "input file '$file' does not exist" unless -e $file;
  my $tmpfile = $file;
  $tmpfile =~ s/\.gz$//;
  my ($basename, $path, $suffix) = fileparse($tmpfile, qr/\.[^.]*/);
  die "input files may not have .fa suffixes" if ($suffix eq '.fa');
}

foreach my $file (@ARGV) {
  print STDERR "processing file '$file'\n";
  if ($file =~ /\.gz$/) {
    system("gunzip $file");
    $file =~ s/\.gz$//;
  }
  my ($basename, $path, $suffix) = fileparse($file, qr/\.[^.]*/);
  my $cmd = "steve-parse2fasta_${type}-end.pl $file > $path$basename.fa"; 
  print STDERR "$cmd\n";
  system($cmd) == 0 or die "Failed running $cmd";
  system("gzip $file");
}

sub usage {
  print "
parse2fasta_multi.pl <single|paired> file1 file2 ...

Parse one or more short read sequence files into GGTools compatible
fasta format.

Each input file of the form xxxxxx.yyyy will produce a file of the form
xxxxxx.fa.

Depending on the specified type, calles parse2fasta_single.pl or parse2fasta_paired.pl on each file
";

  exit(1);
}
