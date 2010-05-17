#!/usr/bin/perl

# Written by Gregory R Grant
# University of Pennsylvania, 2010

if(@ARGV < 1 || $ARGV[0] =~ /help/ || $ARGV[0] eq "?") {
    die "
Usage: rumazon_setup.pl -go

This script sets up the rum pipeline on an amazon
cloud instance.  Run with -go in order to do the
full install.  Takes no other parameters.  You
will be queried for the right organism to install.

";
}
if(!($ARGV[1] eq "-go")) {
    exit(0);
}

$STR = "alias rum=\"cd /mnt/vol2/\"\n";
$STR = $STR . "alias lsr=\"ls -ltr\"\n";
$STR = $STR . "alias lsd=\"ls -lrt | grep ^d\"\n";
$STR = $STR . "alias s3cmd=\"~/s3cmd-0.9.9.91/s3cmd\"\n";

$x = `cat .bash_profile`;
$x = $x . "\n$STR";

open(OUTFILE, ">.bash_profile");
print OUTFILE $x;
close(OUTFILE);

`mkdir /mnt/vol2/scripts`;
`mkdir /mnt/vol2/bin`;
`mkdir /mnt/vol2/indexes`;
`mkdir /mnt/vol2/data`;
print STDERR "installing java, please wait...\n";
`yes|yum install java`;
print STDERR "installing emacs, please wait...\n";
`yes|yum install emacs`;

`wget http://sourceforge.net/projects/s3tools/files/s3cmd/0.9.9.91/s3cmd-0.9.9.91.tar.gz/download`;
`gunzip s3cmd-0.9.9.91.tar.gz`;
`tar -xvf s3cmd-0.9.9.91.tar`;
`wget http://itmat.rum.s3.amazonaws.com/.emacs`;
`wget http://itmat.rum.s3.amazonaws.com/rum_pipeline.tar`;
`mv rum_pipeline.tar /mnt/vol2/`;
`tar -C /mnt/vol2 -xvf /mnt/vol2/rum_pipeline.tar`;
`rm /mnt/vol2/rum_pipeline.tar`;
`wget http://itmat.rum.s3.amazonaws.com/organisms.txt`;

$str = `grep "start \-\-" organisms.txt`;
@organisms = split(/\n/,$str);
$num_organisms = @organisms;
print "\n";
print "--------------------------------------\n";
print "The following organisms are available:\n\n";
for($i=0; $i<@organisms; $i++) {
    $organisms[$i] =~ s/^-- //;
    $organisms[$i] =~ s/ start --$//;
    $j = $i+1;
    print "($j) $organisms[$i]\n";
}
print "--------------------------------------\n";
print "\n";
print "Enter the number of the organism you want to install: ";
$orgnumber = <STDIN>;
chomp($orgnumber);
print "\n";
while(!($orgnumber =~ /^\d+$/) || ($orgnumber <= 0) || ($orgnumber > $num_organisms)) {
    print "Please enter a number between 1 and $num_organisms: ";
    $orgnumber = <STDIN>;
}
$orgnumber--;
print "You have chosen organism $organisms[$orgnumber]\n\n";

open(INFILE, "organisms.txt");
$line = <INFILE>;
chomp($line);
$org = "$organisms[$orgnumber]";
until($line =~ /-- $org start --/) {
    $line = <INFILE>;
    chomp($line);
}
$line = <INFILE>;
chomp($line);
until($line =~ /-- $org end --/) {
    print "$line\n";
    $file = $line;
    $file =~ s!.*/!!;
    if($line =~ /rum.config/) {
	`wget $line`;
	`mv $file /mnt/vol2/$file`;
    } else {
	`wget $line`;
	`mv $file /mnt/vol2/indexes/$file`;
    }
    $line = <INFILE>;
    chomp($line);
}
print "\n";
$x = `ls /mnt/vol2/indexes/ | grep gz`;
if($x =~ /\S/) {
    `gunzip /mnt/vol2/indexes/*gz`;
}
