#!/usr/bin/perl -w
#simR-scape.pl

use strict;
use Class::Struct;
use lib '/Users/rivase/projects/R-scape/scripts';
use FUNCS;
#use constant GNUPLOT => '/usr/bin/gnuplot';
use constant GNUPLOT => '/opt/local/bin/gnuplot';

use vars qw ($opt_a $opt_A $opt_c $opt_C $opt_E $opt_f $opt_g $opt_i $opt_K $opt_N $opt_s $opt_t $opt_T $opt_u $opt_v );  # required if strict used
use Getopt::Std;
getopts ('a:A:cCE:fgiK:N:st:T:u:v');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  simR-scape.pl [options] <msafile> \n\n";
        print "options:\n";
 	print "-a <x>    : set abl to <x>\n";
 	print "-A <x>    : set atbl to <x> \n";
 	print "-c        : C2  type\n";
 	print "-C        : C16 type\n";
 	print "-E <x>    : E-value is <x>\n";
	print "-i        : noindels flag is on\n";
 	print "-f        : draw figures\n";
 	print "-g        : debugging flag\n";
 	print "-K <n>    : # of simulations is <n>\n";
	print "-N <n>    : # seqs in alignment is <n>\n";
	print "-t <s>    : covariation type is <s>\n";
 	print "-T <s>    : tree type is <s> (options are: 'sim' 'star' or 'given'\n";
	print "-u <n>    : use sq number <n> as root\n";
 	print "-s        : do not use ss\n";
 	print "-v        : be verbose\n";
 	exit;
}
my $omsafile  = shift;
my $msaname = "$omsafile"; if ($msaname =~ /([^\/]+).sto$/) { $msaname = $1; }

#programs
my $rscape     = "~/src/src/mysource/src/R-scape";
my $rscape_sim = "~/src/src/mysource/src/R-scape-sim";

# options for R-scape-sim
my $N         = 40;    if ($opt_N) { $N = $opt_N; }
my $atbl      = -1;    if ($opt_A) { $atbl = $opt_A; }
my $abl       = -1;    if ($opt_a) { $abl  = $opt_a; }
my $treetype  = "sim"; if ($opt_T) { $treetype = "$opt_T"; } if ($treetype =~ /^sim$/ || $treetype =~ /^star$/ || $treetype =~ /^rand$/ || $treetype =~ /^given$/) {  } else { print "bad treetype $treetype\n"; die; }
my $noss      = 0;     if ($opt_s) { $noss = 1; }
my $noindels  = 0;     if ($opt_i) { $noindels = 1; }
my $gdb       = 0;     if ($opt_g) { $gdb = 1; }

# options for R-scape
my $Eval      = 0.05;  if ($opt_E) { $Eval = $opt_E; } if ($Eval < 0) { print "bad Evalue $Eval\n"; die; }
my $covtype   = "GTp"; if ($opt_t) { $covtype = "$opt_t"; }
my $nofigures = 1;     if ($opt_f) { $nofigures = 0; }
my $isC16     = 0;     if ($opt_C) { $isC16 = 1; }
my $isC2      = 0;     if ($opt_c) { $isC2  = 1; }
my $usesq     = -1;    if ($opt_u) { $usesq = $opt_u; }

# options for both
my $verbose  = 0;      if ($opt_v) { $verbose = 1; }
 
my $simsafile = "$msaname.synthetic_N$N\_$treetype.sto"; 
my $outfile   = "$msaname.synthetic_N$N\_$treetype.out"; 

my $K = 20;
my $tau;
my $Fm;
my $mean_tau = 0.0;
my $stdv_tau = 0.0;
my $mean_F   = 0.0;
my $stdv_F   = 0.0;

if ($opt_K) { $K = $opt_K; }
for (my $k = 0; $k < $K; $k ++) {
    run_rscapesim($N, $abl, $atbl, $treetype, $noss, $noindels, $usesq, $gdb, $omsafile, $simsafile, $verbose);
    run_rscape ($Eval, $covtype, $nofigures, $isC2, $isC16, $outfile, $simsafile, $verbose);
    if ($verbose) { system("more $outfile\n"); }

    parse_outfile($outfile, \$tau, \$Fm);
    if ($tau < 0) { print "bad tau $tau\n"; }
    if ($Fm  < 0) { print "bad Fm  $Fm\n"; die; }
    $mean_F   += $Fm;
    $stdv_F   += $Fm*$Fm;
    $mean_tau += $tau;
    $stdv_tau += $tau*$tau;
}

FUNCS::calculate_averages(\$mean_F,   \$stdv_F, $K);
FUNCS::calculate_averages(\$mean_tau, \$stdv_tau, $K);

printf("tau %f +/- %f\n", $mean_tau, $stdv_tau);
printf("F   %f +/- %f\n", $mean_F, $stdv_F);

sub
run_rscapesim {
    my ($N, $abl, $atbl, $treetype, $noss, $noindels, $usesq, $gdb, $omsafile, $simsafile, $verbose) = @_;
    
    my $cmd = ($gdb)? "gdb --args $rscape_sim ": "$rscape_sim ";
    if ($N     > 0) { $cmd .= "-N $N "; }
    if ($abl   > 0) { $cmd .= "--abl $abl "; }
    if ($atbl  > 0) { $cmd .= "--atbl $atbl "; }
    if ($usesq > 0) { $cmd .= "--usesq $usesq "; }
    if ($noss)      { $cmd .= "--noss "; }
    if ($noindels)  { $cmd .= "--noindels "; }

    $cmd .= "--$treetype ";
    $cmd .= "-o $simsafile $omsafile";

    printf("$cmd\n");
    system("$cmd\n");
}



sub
run_rscape {
    my ($Eval, $covtype, $nofigures, $isC2, $isC16, $outfile, $msafile, $verbose) = @_;
    my $cmd = "$rscape ";
    $cmd .= "-E $Eval "; 
    if ($isC2)     { $cmd .= "--C2 "; }
    if ($isC16)    { $cmd .= "--C16 "; }
    if ($nofigures) { $cmd .= "--nofigures "; }
    $cmd .= "-o $outfile $msafile";

    system("$cmd\n");
}

sub 
parse_outfile {
    my ($outfile, $ret_tau, $ret_F) = @_;

    my $tau = -1;
    my $F   = -1;

    open (OUT, "$outfile") || die; 
    while (<OUT>) {
	if (/\# GammaFIT:.+tau (\S+)\s*/)  { $tau = $1; }
	if (/\#.+thresh.+\s+(\S+)]\s*$/)   { $F   = $1; }
    }

    $$ret_F   = $F;
    $$ret_tau = $tau;
}
