#!/usr/bin/perl -w
#pdb_parse.pl

use strict;
use Class::Struct;

use vars qw ($opt_C $opt_M $opt_R $opt_W $opt_D $opt_L $opt_v);  # required if strict used
use Getopt::Std;
getopts ('C:M:RW:D:L:v');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  pdb_parse.pl [options] <pdbfile> <stofile> <rscapebin> <gnuplotdir> \n\n";
        print "options:\n";
 	exit;
}

my $pdbfile   = shift;
my $stofile   = shift;
my $rscapebin = shift;
my $gnuplot   = shift;

#my $lib;
#BEGIN { $lib = "$rscapebin/../scripts" };
#use lib '$lib';
use lib '/Users/rivase/src/src/mysource/scripts';
use PDBFUNCS;
use FUNCS;
use constant GNUPLOT => '/usr/local/bin/gnuplot';

my $coorfile = "";
if ($opt_C) { $coorfile = "$opt_C";}
my $mapallfile = "";
if ($opt_M) { $mapallfile = "$opt_M"; }

my $maxD = 8;
if ($opt_D) { $maxD = $opt_D; }
$maxD = int($maxD*100)/100;
my $minL = 1;
if ($opt_L) { $minL = $opt_L; }

my $dornaview = 0;
if ($opt_R) { $dornaview = 1; }

#options: CA CB C MIN AVG NOH / C1' (for RNA suggested by Westhof)
my $which = "MIN";
if ($opt_W) { $which = "$opt_W"; }

my $seeplots = 1;

my $ncnt_t = 0; ## total contacts from all chains
my @cnt_t;
my $msalen;

PDBFUNCS::contacts_from_pdbfile ($gnuplot, $rscapebin, $pdbfile, $stofile, \$msalen, \$ncnt_t, \@cnt_t, $maxD, $minL, 
				 $which, $dornaview, $coorfile, $mapallfile, $seeplots);
