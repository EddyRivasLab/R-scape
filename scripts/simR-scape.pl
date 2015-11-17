#!/usr/bin/perl -w
#simR-scape.pl

use strict;
use Class::Struct;
use lib '/Users/rivase/projects/R-scape/scripts';
use FUNCS;
#use constant GNUPLOT => '/usr/bin/gnuplot';
use constant GNUPLOT => '/opt/local/bin/gnuplot';

use vars qw ($opt_a $opt_A $opt_c $opt_C $opt_E $opt_f $opt_g $opt_i $opt_K $opt_N $opt_s $opt_t $opt_T $opt_u $opt_v $opt_V);  # required if strict used
use Getopt::Std;
getopts ('a:A:cCE:fgiK:N:st:T:u:vV');

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
 	print "-V        : viewplots\n";
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
my $treetype  = "sim"; if ($opt_T) { $treetype = "$opt_T"; } if ($treetype =~ /^sim$/   || 
								 $treetype =~ /^star$/  || 
								 $treetype =~ /^rand$/  || 
								 $treetype =~ /^given$/ || 
								 $treetype =~ /^all$/     ) {  } else { print "bad treetype $treetype\n"; die; }
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
my $K         = 20;    if ($opt_K) { $K = $opt_K; }
my $verbose   = 0;     if ($opt_v) { $verbose = 1; }
my $viewplots = 0;     if ($opt_V) { $viewplots = 1; }

if ($treetype =~ /^all$/) {
    if ($usesq < 0) { $usesq = 1 + int(rand($N)); } # so all treetype methods start from the same ancestral
    
    run_for_treetype("rand", $K, $omsafile, $verbose,
		     $N, $abl, $atbl, $noss, $noindels, $usesq, $gdb, 
		     $Eval, $covtype, $nofigures, $isC2, $isC16);
    run_for_treetype("star", $K, $omsafile, $verbose,
		     $N, $abl, $atbl, $noss, $noindels, $usesq, $gdb, 
		     $Eval, $covtype, $nofigures, $isC2, $isC16);
    run_for_treetype("sim", $K, $omsafile, $verbose,
		     $N, $abl, $atbl, $noss, $noindels, $usesq, $gdb, 
		     $Eval, $covtype, $nofigures, $isC2, $isC16);
    run_for_treetype("given", $K, $omsafile, $verbose,
		     $N, $abl, $atbl, $noss, $noindels, $usesq, $gdb, 
		     $Eval, $covtype, $nofigures, $isC2, $isC16);
}
else {
    run_for_treetype($treetype, $K, $omsafile, $verbose,
		     $N, $abl, $atbl, $noss, $noindels, $usesq, $gdb, 
		     $Eval, $covtype, $nofigures, $isC2, $isC16);
}


sub run_for_treetype {
    my ($treetype, $K, $omsafile, $verbose,
	$N, $abl, $atbl, $noss, $noindels, $usesq, $gdb, 
	$Eval, $covtype, $nofigures, $isC2, $isC16) = @_;
    
    my $simsafile = ($noss)? "$msaname\_synthetic_N$N\_$treetype.noss.sto"   : "$msaname\_synthetic_N$N\_$treetype.sto"; 
    my $outfile   = ($noss)? "$msaname\_synthetic_N$N\_$treetype.noss.out"   : "$msaname\_synthetic_N$N\_$treetype.out"; 
    my $hisfile   = ($noss)? "$msaname\_synthetic_N$N\_$treetype.cyk.his.ps" : "$msaname\_synthetic_N$N\_$treetype.his.ps"; 

    my $taufile   = ($noss)? "$msaname\_synthetic_N$N\_$treetype.noss.tauhis"    : "$msaname\_synthetic_N$N\_$treetype.tauhis";
    my $taups     = ($noss)? "$msaname\_synthetic_N$N\_$treetype.noss.tauhis.ps" : "$msaname\_synthetic_N$N\_$treetype.tauhis.ps";

    my $maxscfile = ($noss)? "$msaname\_synthetic_N$N\_$treetype.noss.maxschis"    : "$msaname\_synthetic_N$N\_$treetype.maxschis";
    my $maxscps   = ($noss)? "$msaname\_synthetic_N$N\_$treetype.noss.maxschis.ps" : "$msaname\_synthetic_N$N\_$treetype.maxschis.ps";

    my $tau;
    my $Fm;
    my $maxsc;
    my $avgid;
    my $mean_tau   = 0.0;
    my $stdv_tau   = 0.0;
    my $min_tau    = 123456789;
    my $max_tau    = -123456789;
    my $mean_F     = 0.0;
    my $stdv_F     = 0.0;
    my $mean_maxsc = 0.0;
    my $stdv_maxsc = 0.0;
    my $min_maxsc  = 123456789;
    my $max_maxsc  = -123456789;
    my $mean_avgid = 0.0;
    my $stdv_avgid = 0.0;

    my @his_tau = ();
    my $Nt = 40;
    my $kt = 5;
    FUNCS::init_histo_array($Nt, $kt, \@his_tau);

    my @his_maxsc = ();
    my $Nsc = 500;
    my $ksc = 1;
    FUNCS::init_histo_array($Nsc, $ksc, \@his_maxsc);
    
    for (my $x = 0; $x < $K; $x ++) {
	run_rscapesim($N, $abl, $atbl, $treetype, $noss, $noindels, $usesq, $gdb, $omsafile, $simsafile, $verbose);
	run_rscape ($Eval, $covtype, $nofigures, $isC2, $isC16, $outfile, $simsafile, $verbose);
	if ($verbose) { system("more $outfile\n"); }
	
	parse_outfile($outfile, $Nsc, $ksc, \@his_maxsc, \$maxsc, $Nt, $kt, \@his_tau, \$tau, \$Fm, \$avgid);
	if ($tau < 0) { print "bad tau $tau\n"; }
	if ($Fm  < 0) { print "bad Fm  $Fm\n"; die; }
	$mean_F     += $Fm;
	$stdv_F     += $Fm*$Fm;
	$mean_tau   += $tau;
	$stdv_tau   += $tau*$tau;
 	$mean_maxsc += $maxsc;
	$stdv_maxsc += $maxsc*$maxsc;
	$mean_avgid += $avgid;
	$stdv_avgid += $avgid*$avgid;
	if ($tau   < $min_tau)   { $min_tau   = $tau;   }
	if ($tau   > $max_tau)   { $max_tau   = $tau;   }
	if ($maxsc < $min_maxsc) { $min_maxsc = $maxsc; }
	if ($maxsc > $max_maxsc) { $max_maxsc = $maxsc; }
    }
    FUNCS::calculate_averages(\$mean_F,     \$stdv_F,     $K);
    FUNCS::calculate_averages(\$mean_tau,   \$stdv_tau,   $K);
    FUNCS::calculate_averages(\$mean_maxsc, \$stdv_maxsc, $K);
    FUNCS::calculate_averages(\$mean_avgid, \$stdv_avgid, $K);
    printf("#tau   %f +/- %f\n", $mean_tau,   $stdv_tau);
    printf("#F     %f +/- %f\n", $mean_F,     $stdv_F);   
    printf("#maxsc %f +/- %f [%f,%f]\n", $mean_maxsc, $stdv_maxsc, $min_maxsc, $max_maxsc);   
    printf("#avgid %f +/- %f\n", $mean_avgid, $stdv_avgid);   
    
    FUNCS::write_histogram($Nt,  $kt,  0, \@his_tau,   1.0, $taufile,   0);
    FUNCS::write_histogram($Nsc, $ksc, 0, \@his_maxsc, 1.0, $maxscfile, 0);
    open(HIS, ">>$taufile");
    print HIS "# mean_avgid $mean_avgid stdv_avgid $stdv_avgid\n";
    close(HIS);
    open(HIS, ">>$maxscfile");
    print HIS "# mean_avgid $mean_avgid stdv_avgid $stdv_avgid\n";
    close(HIS);

    my $key = "avgid=$mean_avgid +/- $stdv_avgid";
    FUNCS::gnuplot_histo($taufile,   1, 2, $taups,   $key, "tau",        "ocurrence", "$treetype", 0, 1, $min_tau-0.2, $max_tau+0.2, $Nt/2.);
    FUNCS::gnuplot_histo($maxscfile, 1, 2, $maxscps, $key, "max cov sc", "ocurrence", "$treetype", 0, 1, $min_maxsc-5, $max_maxsc+5, $Nsc/2.);
    
    if ($viewplots && $K == 1) { system("more $outfile\n"); system("open $hisfile\n"); }
    system("rm *sto\n");
    system("rm *out\n");
    system("rm *sorted\n");
    system("rm *roc\n");
    system("rm *sum\n");
    system("rm *out\n");
 }

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
    my $cmd = "time $rscape ";
    $cmd .= "-E $Eval "; 
    if ($isC2)     { $cmd .= "--C2 "; }
    if ($isC16)    { $cmd .= "--C16 "; }
    if ($nofigures) { $cmd .= "--nofigures "; }
    $cmd .= "-o $outfile $msafile";

    system("$cmd\n");
}

sub 
parse_outfile {
    my ($outfile, $Nsc, $ksc, $his_maxsc_ref, $ret_maxsc, $Nt, $kt, $his_tau_ref, $ret_tau, $ret_F, $ret_avgid) = @_;

    my $maxsc = -1;
    my $tau   = -1;
    my $F     = -1;
    my $avgid = -1;

    open (OUT, "$outfile") || die; 
    while (<OUT>) {
	if    (/\# GammaFIT:.+tau (\S+)\s*/)  { $tau = $1; }
	elsif (/\# .+avgid (\S+)\s*/)         { $avgid = $1; }
	elsif (/\#.+thresh.+\s+cov=\S+\s\[\S+\,(\S+)\]\s\[.+(\S+)\]\s*$/) { $maxsc = $1; $F = $2; }
    }
    close(OUT);

    FUNCS::fill_histo_array(1, $tau,   $Nt,  $kt,  0, $his_tau_ref);
    FUNCS::fill_histo_array(1, $maxsc, $Nsc, $ksc, 0, $his_maxsc_ref);
    $$ret_F     = $F;
    $$ret_tau   = $tau;
    $$ret_maxsc = $maxsc;
    $$ret_avgid = $avgid;
}
