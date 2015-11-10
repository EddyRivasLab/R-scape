#!/usr/bin/perl -w
#simR-scape.pl

use strict;
use Class::Struct;
use lib '/Users/rivase/projects/R-scape/scripts';
use FUNCS;
#use constant GNUPLOT => '/usr/bin/gnuplot';
use constant GNUPLOT => '/opt/local/bin/gnuplot';

use vars qw ($opt_c $opt_C $opt_E $opt_f $opt_t $opt_T $opt_v );  # required if strict used
use Getopt::Std;
getopts ('cCE:ft:T:v');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  simR-scape.pl [options] <msafile> \n\n";
        print "options:\n";
 	print "-c <s>    : C2  type\n";
 	print "-C <s>    : C16 type\n";
 	print "-E <x>    : E-value is <x>\n";
	print "-t <s>    : covariation type is <s>\n";
 	print "-T <s>    : tree type is <s> (options are: 'sim' 'star' or 'given'\n";
 	print "-v    :  be verbose\n";
 	exit;
}
my $omsafile  = shift;
my $simsafile = "$omsafile.sim";

#programs
my $rscape     = "~/src/src/mysource/src/R-scape";
my $rscape_sim = "~/src/src/mysource/src/R-scape-sim";

# options for R-scape-sim
my $N;
my $atbl;
my $treetype = "sim";
if ($opt_T) { $treetype = "$opt_T"; }
if ($treetype =~ /^sim$/ || $treetype =~ /^star$/ $treetype =~ /^given$/) { next; }
else { print "bad treetype $treetype\n"; die; }
my $noss = 0;
my $noindels = 0;

# options for R-scape
my $Eval = 0.05;
if ($opt_E) { $Eval = $opt_E; }
if ($Eval < 0) { print "bad Evalue $Eval\n"; die; }
my $covtype = "GTp";
if ($opt_t) { $covtype = "$opt_t"; }
my $nofigures = 1;
if ($opt_f) { $nofigures = 0; }
my $isC16 = 0; if ($opt_C) { $isC16 = 1; }
my $isC2  = 0; if ($opt_c) { $isC2  = 1; }

# options for both
my $verbose = 0;
if ($opt_v) { $verbose = 1; }

run_rscape_sim($N, $atbl, $treetype, $noss, $noindels, $omsafile, $simsafile, $verbose);
run_rscape_sim($Eval, $covtype, $nofigures, $isC2, $isC16, $simsafile, $verbose);



sub
run_rscapesim {
    my ($N, $atbl, $treetype, $noss, $noindels, $nofigures, $omsafile, $simsafile, $verbose) = @_;
}



sub
run_rscape {
    my ($Eval, $covtype, $nofigures, $msafile, $isC2, $isC16, $verbose) = @_;
}
