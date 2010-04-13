#!/usr/bin/perl
use strict;
use Data::Dumper;

# before we get started, make sure all files exist
foreach my $file (@ARGV) {
  die "" unless -e $file;
}

# each block of N lines should have one line of sequence
# in the same place within each block.
# print out that line in fasta format, with a generated ID.
# also validate that each block has the same number of lines
# and that the line w/ sequence is in the same place in the block
foreach my $file (@ARGV) {
  open(F, $file) || die "can't open file '$file'\n";
  my $blockSizeHash = {};
  my $seqCount = 0;
  my $prevSeqLineNum;
  my $initialLinesCount = 1;
  my $foundFirstSeq = 0;
  while (<F>) {
    s/\s//;			# lose all white space
    s/[^ACGTN]$//;		# lose trailing bogus character?
    if (!$_ || /[^ACGTN]/) {
      $initialLinesCount++ unless $foundFirstSeq;
      next;
    }
    validate($file, $blockSizeHash, $. - $prevSeqLineNum) if $foundFirstSeq;
    $foundFirstSeq = 1;
    $seqCount++;
    print ">seq$seqCount\n$_\n";
    $prevSeqLineNum = $.;
  }
  validate($file, $blockSizeHash, $initialLinesCount + $. - $prevSeqLineNum);

  die "Didn't find any sequences\n" if (!$seqCount);

  close(F);
}

sub validate {
  my ($file, $blockSizeHash, $key) = @_;
  $blockSizeHash->{$key}++;

  # a valid file will have tracked all blocks into one hash key of the
  # standard size block
  if (scalar(keys(%$blockSizeHash)) != 1) {
    print STDERR "File $file seems to be invalid at line $..\n block count hash:\n";
    print STDERR Dumper($blockSizeHash) . "\n\n";
    exit(1);
  }
}

