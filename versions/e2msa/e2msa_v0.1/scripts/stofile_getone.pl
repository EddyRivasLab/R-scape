#!/usr/bin/perl -w
#stofile_getone.pl
 
use strict;
use Class::Struct;
use lib '/groups/eddy/home/rivase/projects/evohmm/scripts';
use lib '/Users/rivase/projects/evohmm/scripts';
use FUNCS;

use vars qw ();  # required if strict used
use Getopt::Std;
getopts ('');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  stofile_getone.pl [options] which <stofile>  \n\n";
        print "options:\n";
	exit;
}

my $which = shift;
my $stofile = shift;

my $nmsa = 0;
my $thisone = 0;

open(STO,   "$stofile")   || die;
while (<STO>) {
    if (/\# STOCKHOLM 1.0/) {
	$nmsa ++;

	if ($nmsa == $which) {
	    $thisone = 1;
	    print $_;
	}

    }
    elsif (/^\/\//) { 
	if ($thisone) { print $_; }
	$thisone = 0;
    }
    elsif ($thisone) {
	print $_;
    }
}

close(STO);
