#!/usr/bin/perl

# Written by Gregory R Grant
# Universtiy of Pennsylvania, 2010

use strict;

if(@ARGV < 2) {
    die "
Usage: smoother.pl <infile> <outfile> [option]

Where:
<infile> is a spreadsheet of numerical data with the a first column of id's.

Option:
   -header  : infile has a header line

";
}
my $header = 'false';
if($ARGV[2] eq '-header') {
    $header = 'true';
}

open (OUTFILE, ">$ARGV[1]");
my @test_vector;
my $infile = $ARGV[0];
my $a=`cat $infile`;
my @originaldata = split(/\n/,$a);
my $numberdatarows=@originaldata;

my $start;
if($header eq 'true') {
    print "$originaldata[1]\n";
    $start = 2;
} else {
    $start = 1;
}

for(my $i=$start; $i<$numberdatarows; $i++) {
    my @X_vect;
    chomp($originaldata[$i]);
    my @rowdata=split(/\t/,$originaldata[$i]);
    my $rowlength=@rowdata;
    my $k=0;
    my @inputvector;
    my $line;
    for(my $y=1; $y<$rowlength; $y++) {
	if($rowdata[$y] ne "NA") {
	    $inputvector[$k]=$rowdata[$y];
	    $X_vect[$k]=$k+1;
	    $k++;
	}
    }
    my $smoothed_vector_ref= Smooth(\@inputvector);    
    my @smoothed_vector = @{$smoothed_vector_ref};
    
    $line="$rowdata[0]";
    $k=0;
    for(my $y=1; $y<$rowlength; $y++) {
	if($rowdata[$y] ne "NA") {
	    my $val = int(10000*$smoothed_vector[$k])/10000;
	    $line = $line . "\t$val";
	    $k++;
	}
	else {
	    $line = $line . "\tNA";
	}
    }
    $line = $line . "\n";

    print OUTFILE "$line";
#    print "$line";

}


sub Smooth {
    my ($input_vector_ref) = @_;
    my @input_vector = @{$input_vector_ref};
    
    my $vector_length=@input_vector;

    if($vector_length<6) {
	die "vector is too short to smooth, must be at least length six\n";
    }
    my @temp1;
    my @temp2;
    my $a;
    my $b;
    my @smoothed_vector;

#  the case 0

    $temp1[0]=$input_vector[0];
    $temp1[1]=$input_vector[1];
    $temp1[2]=$input_vector[2];
    $temp2[0]=1;
    $temp2[1]=2;
    $temp2[2]=3;

    ($a, $b) = LS(\@temp2, \@temp1);
    $smoothed_vector[0]=$a+$b;
    
#  the case n (n=vector length)

    $temp1[0]=$input_vector[$vector_length-3];
    $temp1[1]=$input_vector[$vector_length-2];
    $temp1[2]=$input_vector[$vector_length-1];
    ($a, $b) = LS(\@temp2, \@temp1);
    $smoothed_vector[$vector_length-1]=$a*3+$b;

#  the case 1

    $temp1[0]=$input_vector[0];
    $temp1[1]=$input_vector[1];
    $temp1[2]=$input_vector[2];
    $temp1[3]=$input_vector[3];
    $temp2[0]=1;
    $temp2[1]=2;
    $temp2[2]=3;
    $temp2[3]=4;

    ($a, $b) = LS(\@temp2, \@temp1);

    $smoothed_vector[1]=$a*2+$b;

#  the case n-1 (n=vector length)

    $temp1[0]=$input_vector[$vector_length-4];
    $temp1[1]=$input_vector[$vector_length-3];
    $temp1[2]=$input_vector[$vector_length-2];
    $temp1[3]=$input_vector[$vector_length-1];
    ($a, $b) = LS(\@temp2, \@temp1);
    $smoothed_vector[$vector_length-2]=$a*3+$b;

# the rest of the cases

    for(my $i=2; $i<$vector_length-2; $i++) {

	$temp1[0]=$input_vector[$i-2];
	$temp1[1]=$input_vector[$i-1];
	$temp1[2]=$input_vector[$i];
	$temp1[3]=$input_vector[$i+1];
	$temp1[4]=$input_vector[$i+2];
	$temp2[4]=5;

	($a, $b) = LS(\@temp2, \@temp1);
	$smoothed_vector[$i]=$a*3+$b;
    }

    return (\@smoothed_vector);
}    
    
sub LS {
    my ($X_vect_ref, $Y_vect_ref)=@_;
    my @X_vect=@{$X_vect_ref};
    my @Y_vect=@{$Y_vect_ref};
    
    my $X_vect_length=@X_vect;
    my $Y_vect_length=@Y_vect;

    if($X_vect_length != $Y_vect_length) {
	"ERROR: tried to run the LS subroutine with vectors of differing lengths\n";
    }

    my $y=0;
    my $x=0;
    my $yy=0;
    my $xx=0;
    my $xy=0;

    for(my $i=0; $i<$X_vect_length;$i++) {
	$x=$x+$X_vect[$i];
	$y=$y+$Y_vect[$i];
	$xx=$xx+($X_vect[$i])*($X_vect[$i]);
	$yy=$yy+($Y_vect[$i])*($Y_vect[$i]);
	$xy=$xy+($X_vect[$i]*$Y_vect[$i]);
    }    

    my $a=($y*$xx-$x*$xy)/($X_vect_length*$xx-$x*$x);
    my $b=($X_vect_length*$xy-$x*$y)/($X_vect_length*$xx-$x*$x);

    return($b,$a);
}
