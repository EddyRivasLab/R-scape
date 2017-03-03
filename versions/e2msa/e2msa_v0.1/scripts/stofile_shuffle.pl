#!/usr/bin/perl -w
#stofile_shuffle.pl
 
use strict;
use Class::Struct;
use lib '/groups/eddy/home/rivase/projects/evohmm/scripts';
use lib '/Users/rivase/projects/evohmm/scripts';
use FUNCS;
use List::Util qw(shuffle);

use vars qw ();  # required if strict used
use Getopt::Std;
getopts ('');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  stofile_shuffle.pl [options] <stofile>  \n\n";
        print "options:\n";
	exit;
}

my $stofile = shift;

my $nmsa = 0;
my @msa = ();
my @order = ();

open(STO,   "$stofile")   || die;
while (<STO>) {
    if (/\# STOCKHOLM 1.0/) {
	$msa[$nmsa] = $_;
	$order[$nmsa] = $nmsa;
    }
    elsif (/^\/\//) { 
	$msa[$nmsa] .= $_; 
	$nmsa ++;
    }
    else {
	$msa[$nmsa] .= $_;
    }
}
close(STO);

@order =  shuffle @order;

my $newn = $#order+1;
if ($nmsa != $newn) { print " $newn shoudl be $nmsa\n"; die; }

for (my $n = 0; $n < $nmsa; $n ++)
{
     print $msa[$order[$n]];
}
