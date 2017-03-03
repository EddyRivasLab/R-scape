#!/usr/bin/perl -w
#plot_evohmm.pl

use strict;
use Class::Struct;
use lib '/groups/eddy/home/rivase/projects/evohmm/scripts';
use lib '/Users/rivase/projects/evohmm/scripts';
use FUNCS;
use constant GNUPLOT => '/usr/bin/gnuplot';
#use constant GNUPLOT => '/opt/local/bin/gnuplot';

use vars qw ($opt_k $opt_I $opt_l $opt_L $opt_g $opt_n $opt_p $opt_r $opt_s $opt_t $opt_u $opt_U $opt_v $opt_x);  # required if strict used
use Getopt::Std;
getopts ('k:I:l:L:u:U:gnpr:st:vx:');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  plot_evohmm.pl [options] <srcdir> <wrkdir> \n\n";
        print "options:\n";
        print "-s    :  don't plot slopes       \n";
        print "-p    :  program is 'utest_evoprobs' [default gappenalties'       \n";
	print "-r    :  rz value       \n";
	print "-x    :  etaz value      \n";
	print "-I    :  etainf value \n";
	print "-l    :  ld \n";
	print "-L    :  ldI \n";
	print "-u    :  muE \n";
	print "-U    :  muI \n";
	print "-k    :  kappa=muE(1-rx)/(muI-rz*ldI() \n";
	print "-g    :  log scale for probabilities \n";
	print "-n    :  don't see plots\n";
 	print "-v    :  be verbose\n";
        print "-t    :  tmax [default is 3]\n";
	exit;
}
my $srcdir   = shift;
my $wrkdir   = shift;

my $tP10 = 0.114220;
my $tP30 = 0.332924;
my $tB90 = 0.903079;
my $tB62 = 1.427753;

my $GOPEN = 0; # TRUE is if you want to plot beta(1-eta)/eta
 
my $filename = $wrkdir;
if ($filename =~ /^\S+\/([^\/]+)$/) { $filename = $1; }

if ( -d $wrkdir) { print "wrkdir directory $wrkdir exist!\n"; system("rm $wrkdir/*\n"); }
else             { system("mkdir $wrkdir\n"); }

my $gappenalties = "$srcdir/src/gappenalties";
my $evoprobs     = "$srcdir/src/utest_evoprobs";

my $program = $gappenalties;
if ($opt_p) { $program = $evoprobs; }

my $plot_slope = 0;
if ($opt_s) { $plot_slope = 0; }
my $view = 1;
if ($opt_n) { $view = 0; }
my $verbose = 0;
if ($opt_v) { $verbose = 1; }

my $plot_default_only = 0;

my $tmax = -1; if ($opt_t) { $tmax = $opt_t; if ($tmax < 0) { print "tmax has to be positive\n"; die; } }

# common variables 
my $etaz       = 0.01; 
my $rz         = 0.01;
 
# variables specific for gappenalties
my $etainf  = 0.99;
my $kappa  = -1.0;
# variables specific for utest_evoprobs
my $ld  = -1.0;
my $muE = -1.0;
my $ldI = -1.0;
my $muI = -1.0;

if ($opt_x) { $etaz       = $opt_x; $plot_default_only = 1; if ( $etaz    < 0 || $etaz      > 1.0) { print "bad value of etaz=$etaz\n";       die; } }
if ($opt_r) { $rz         = $opt_r; $plot_default_only = 1; if ( $rz      < 0 || $rz        > 1.0) { print "bad value of rz=$rz\n";           die; } }
if ($opt_I) { $etainf     = $opt_I; $plot_default_only = 1; if ( $etainf  < 0                    ) { print "bad value of etainf=$etainf\n";   die; } }
if ($opt_k) { $kappa      = $opt_k; $plot_default_only = 1;                                                                                          }
if ($opt_l) { $ld         = $opt_l; $plot_default_only = 1;                                                                                          }
if ($opt_L) { $ldI        = $opt_L; $plot_default_only = 1;                                                                                          }
if ($opt_u) { $muE        = $opt_u; $plot_default_only = 1;                                                                                          }
if ($opt_U) { $muI        = $opt_U; $plot_default_only = 1;                                                                                          }


my $etaz1 = 0.0000001;
my $etaz2 = 0.001;
my $etaz3 = 0.10;
my $etaz4 = 0.15;
my $etaz5 = 0.20;

my $rz1 = 0.000001;
my $rz2 = 0.10;
my $rz3 = 0.50;
my $rz4 = 0.90;
my $rz5 = 0.9999;

my $etainf5 = 0.50;
my $etainf4 = 0.60;
my $etainf3 = 0.80;
my $etainf2 = 0.99;
my $etainf1 = 1.99;

my $kappa5 = 0.50;
my $kappa4 = 0.60;
my $kappa3 = 0.80;
my $kappa2 = 0.99;
my $kappa1 = 1.99;

my $success;


if ($plot_default_only) { 
    if ($program =~ /gappenalties/) {
	$success = run_gappenalties($program, $wrkdir, $filename, $etaz, $rz, $etainf, $kappa, $verbose);
    }
    else {
	$success = run_evoprobs($program, $wrkdir, $filename, $etaz, $rz, $ld, $muE, $ldI, $muI, $verbose);
    }
    if ($success) { plot_program_one($wrkdir, $view); }
}
else {
   
    if ($program =~ /gappenalties/) {
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz1, $rz, $etainf,  $kappa, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz2, $rz, $etainf,  $kappa, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz3, $rz, $etainf,  $kappa, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz4, $rz, $etainf,  $kappa, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz5, $rz, $etainf,  $kappa, $ld, $verbose);
	
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz1, $etainf,  $kappa, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz2, $etainf,  $kappa, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz3, $etainf,  $kappa, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz4, $etainf,  $kappa, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz5, $etainf,  $kappa, $ld, $verbose);
 
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz,  $etainf1, $kappa, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz,  $etainf2, $kappa, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz,  $etainf3, $kappa, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz,  $etainf4, $kappa, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz,  $etainf5, $kappa, $ld, $verbose);
	
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz,  $etainf,  $kappa1, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz,  $etainf,  $kappa2, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz,  $etainf,  $kappa3, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz,  $etainf,  $kappa4, $ld, $verbose);
	run_gappenalties($gappenalties, $wrkdir, $filename, $etaz, $rz,  $etainf,  $kappa5, $ld, $verbose);
    }
	
    plot_program($wrkdir, 1, $view); # etaz      variations 
    plot_program($wrkdir, 2, $view); # rz        variations 
    plot_program($wrkdir, 3, $view); # etainf  variations 
    plot_program($wrkdir, 4, $view); # kappa variations 
}

sub actual_plot {
    my ($name, $type, $tstar, $plotf1, $plotf2, $plotf3, $plotf4, $plotf5, $plotf_tstar, $view) = @_;

    my $xlabel = "substitutions per site";
    my $ylabel = "";

    my $ymin = -15;

    my $plotf = $plotf1; if ($plotf =~ /^(\S+).plot1$/) { $plotf = $1; }
    my $Pfile  = "$plotf\_prob.ps";
    my $Lfile = "$plotf\_logp.ps";

    my $tB62file = "$plotf.tB62";
    my $lB62file = "$plotf.lB62";
    
    my $tB90file = "$plotf.tB90";
    my $lB90file = "$plotf.lB90";
 
    my $tPAM30file = "$plotf.tP30";
    my $lPAM30file = "$plotf.lP30";

    my $tPAM10file = "$plotf.tP10";
    my $lPAM10file = "$plotf.lP10";

    open(T1, ">$tB62file") || die;
    print T1 "$tB62\t 0\n";
    print T1 "$tB62\t 1\n";
    close(T1);

    open(L1, ">$lB62file") || die;
    print L1 "$tB62\t$ymin\n";
    print L1 "$tB62\t0\n";
    close(L1);

    open(T1, ">$tB90file") || die;
    print T1 "$tB90\t 0\n";
    print T1 "$tB90\t 1\n";
    close(T1);

    open(L1, ">$lB90file") || die;
    print L1 "$tB90\t$ymin\n";
    print L1 "$tB90\t0\n";
    close(L1);
 
    open(T1, ">$tPAM30file") || die;
    print T1 "$tP30\t 0\n";
    print T1 "$tP30\t 1\n";
    close(T1);
    open(L1, ">$lPAM30file") || die;
    print L1 "$tP30\t$ymin\n";
    print L1 "$tP30\t0\n";
    close(L1);
    
    open(L1, ">$lPAM10file") || die;
    print L1 "$tP10\t$ymin\n";
    print L1 "$tP10\t0\n";
    close(L1);
    open(T1, ">$tPAM10file") || die;
    print T1 "$tP10\t 0\n";
    print T1 "$tP10\t 1\n";
    close(T1);
    

  my $k;

    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    $ylabel = "Probabilities";
    print GP "set output '$Pfile'\n";
    print GP "set key right top\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set title \"$name vary $type\n";
    print GP "set ylabel '$ylabel'\n";
    if ($tmax > 0) { print GP "set xrange [0:$tmax]\n"; }
    print GP "set yrange [0:1]\n";

    my $tag;
    my $tag1;
    my $tag2;
    my $tag3;
    my $tag4;
    my $tag5;
    if ($type =~ /^etaz$/) {
	$tag1 = $etaz1;
	$tag2 = $etaz2;
	$tag3 = $etaz3;
	$tag4 = $etaz4;
	$tag5 = $etaz5;
    }
    elsif ($type =~ /^rz$/) {
	$tag1 = $rz1;
	$tag2 = $rz2;
	$tag3 = $rz3;
	$tag4 = $rz4;
	$tag5 = $rz5;
    }
   elsif ($type =~ /^etainf$/) {
	$tag1 = $etainf1;
	$tag2 = $etainf2;
	$tag3 = $etainf3;
	$tag4 = $etainf4;
	$tag5 = $etainf5;
    }

    my $cmd;
    my $field;
    
    $field = 3;
    $tag = "ETA";
    $cmd   = "'$tB62file' title 'BLOSUM62' with lines ls 1, ";
    $cmd  .= "'$tB90file' title 'BLOSUM90' with lines ls 1, ";
    $cmd  .= "'$tPAM30file' title 'PAM30' with lines ls 1, ";
    $cmd  .= "'$tPAM10file' title 'PAM10' with lines ls 1, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag); 
    $cmd .= "'$plotf1' using 1:$field  title '$tag-$tag1' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style_transitions($tag."-2"); 
    $cmd .= "'$plotf2' using 1:$field  title '$tag-$tag2' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style_transitions($tag."-3"); 
    $cmd .= "'$plotf3' using 1:$field  title '$tag-$tag3' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style_transitions($tag."-4"); 
    $cmd .= "'$plotf4' using 1:$field  title '$tag-$tag4' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style_transitions($tag."-5"); 
    $cmd .= "'$plotf5' using 1:$field  title '$tag-$tag5' with lines ls $k, "; 

    if ($plot_slope) {
	$field = 4;
	$k = FUNCS::gnuplot_set_style_transitions($tag."-slope"); 
	$cmd .= "'$plotf1' using 1:$field  with lines ls $k, "; 
	if (0) {
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-2-slope"); 
	    $cmd .= "'$plotf2' using 1:$field  with lines ls $k, "; 
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-3-slope"); 
	    $cmd .= "'$plotf3' using 1:$field  with lines ls $k, "; 
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-4-slope"); 
	    $cmd .= "'$plotf4' using 1:$field   with lines ls $k, "; 
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-5-slope"); 
	    $cmd .= "'$plotf5' using 1:$field  with lines ls $k, "; 
	}
    }

    $field = 6;
    $tag = "BETA";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
    $cmd .= "'$plotf1' using 1:$field  title '$tag-$tag1' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style_transitions($tag."-2"); 
    $cmd .= "'$plotf2' using 1:$field  title '$tag-$tag2' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style_transitions($tag."-3"); 
    $cmd .= "'$plotf3' using 1:$field  title '$tag-$tag3' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style_transitions($tag."-4"); 
    $cmd .= "'$plotf4' using 1:$field  title '$tag-$tag4' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style_transitions($tag."-5"); 
    $cmd .= "'$plotf5' using 1:$field  title '$tag-$tag5' with lines ls $k, "; 

    if ($plot_slope) {
	$field = 7;
	$k = FUNCS::gnuplot_set_style_transitions($tag."-slope"); 
	$cmd .= "'$plotf1' using 1:$field  with lines ls $k, "; 
	if (0) {
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-2-slope"); 
	    $cmd .= "'$plotf2' using 1:$field  with lines ls $k, "; 
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-3-slope"); 
	    $cmd .= "'$plotf3' using 1:$field  with lines ls $k, "; 
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-4-slope"); 
	    $cmd .= "'$plotf4' using 1:$field  with lines ls $k, ";    
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-5-slope"); 
	    $cmd .= "'$plotf5' using 1:$field  with lines ls $k, "; 
	}
    }

    if ($GOPEN) {
	$field = 9;
	$tag = "BETA * (1-ETA) / ETA";
	$k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
	$cmd .= "'$plotf1' using 1:$field  title '$tag-$tag1' with lines ls $k, ";
	$k = FUNCS::gnuplot_set_style_transitions($tag."-2"); 
	$cmd .= "'$plotf2' using 1:$field  title '$tag-$tag2' with lines ls $k, ";
	$k = FUNCS::gnuplot_set_style_transitions($tag."-3"); 
	$cmd .= "'$plotf3' using 1:$field  title '$tag-$tag3' with lines ls $k, ";
	$k = FUNCS::gnuplot_set_style_transitions($tag."-4"); 
	$cmd .= "'$plotf4' using 1:$field  title '$tag-$tag4' with lines ls $k, ";
	$k = FUNCS::gnuplot_set_style_transitions($tag."-5"); 
	$cmd .= "'$plotf5' using 1:$field  title '$tag-$tag5' with lines ls $k, ";
	
	if ($plot_slope) {
	    $field = 10;
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-1-slope"); 
	    $cmd .= "'$plotf1' using 1:$field  with lines ls $k, ";
	    if (0) {
		$k = FUNCS::gnuplot_set_style_transitions($tag."-2-slope"); 
		$cmd .= "'$plotf2' using 1:$field  with lines ls $k, ";
		$k = FUNCS::gnuplot_set_style_transitions($tag."-3-slope"); 
		$cmd .= "'$plotf3' using 1:$field  with lines ls $k, ";
		$k = FUNCS::gnuplot_set_style_transitions($tag."-4-slope"); 
		$cmd .= "'$plotf4' using 1:$field  with lines ls $k, ";   
		$k = FUNCS::gnuplot_set_style_transitions($tag."-5-slope"); 
		$cmd .= "'$plotf5' using 1:$field  with lines ls $k, ";
	    }
	}
    }

    $field = 12;
    $tag = "BETA * (1-ETA)";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
    $cmd .= "'$plotf1' using 1:$field  title '$tag-$tag1' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-2"); 
    $cmd .= "'$plotf2' using 1:$field  title '$tag-$tag2' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-3"); 
    $cmd .= "'$plotf3' using 1:$field  title '$tag-$tag3' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-4"); 
    $cmd .= "'$plotf4' using 1:$field  title '$tag-$tag4' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-5"); 
    $cmd .= "'$plotf5' using 1:$field  title '$tag-$tag5' with lines ls $k, ";
    
    if ($plot_slope) {
	$field = 13;
	$k = FUNCS::gnuplot_set_style_transitions($tag."-1-slope"); 
	$cmd .= "'$plotf1' using 1:$field  with lines ls $k, ";
	if (0) {
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-2-slope"); 
	    $cmd .= "'$plotf2' using 1:$field  with lines ls $k, ";
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-3-slope"); 
	    $cmd .= "'$plotf3' using 1:$field  with lines ls $k, ";
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-4-slope"); 
	    $cmd .= "'$plotf4' using 1:$field  with lines ls $k, ";   
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-5-slope"); 
	    $cmd .= "'$plotf5' using 1:$field  with lines ls $k, ";
	}
    }

    # values at tstar
    $field = 3;
    $k = FUNCS::gnuplot_set_style_transitions("tstar"); 
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k, "; 
    $field = 5;
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k, "; 
    $field = 7;
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k, "; 
    $field = 9;
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k"; 

    $cmd .= "\n";
    print GP "plot $cmd\n";


    $ylabel = "Log Probabilities";
    print GP "set output '$Lfile'\n";
    $cmd   = "'$lB62file' title 'BLOSUM62' with lines ls 1, ";
    $cmd  .= "'$lB90file' title 'BLOSUM90' with lines ls 1, ";
    $cmd  .= "'$lPAM30file' title 'PAM30' with lines ls 1, ";
    $cmd  .= "'$lPAM10file' title 'PAM10' with lines ls 1, ";
    print GP "set yrange [$ymin:0]\n";
    print GP "set ylabel '$ylabel'\n";
    
    $field = 2;
    $tag = "ETA";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
    $cmd .= "'$plotf1' using 1:$field  title '$tag-$tag1' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-2"); 
    $cmd .= "'$plotf2' using 1:$field  title '$tag-$tag2' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-3"); 
    $cmd .= "'$plotf3' using 1:$field  title '$tag-$tag3' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-4"); 
    $cmd .= "'$plotf4' using 1:$field  title '$tag-$tag4' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-5"); 
    $cmd .= "'$plotf5' using 1:$field  title '$tag-$tag5' with lines ls $k, ";
 
    $field = 5;
    $tag = "BETA";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
    $cmd .= "'$plotf1' using 1:$field  title '$tag-$tag1' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-2"); 
    $cmd .= "'$plotf2' using 1:$field  title '$tag-$tag2' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-3"); 
    $cmd .= "'$plotf3' using 1:$field  title '$tag-$tag3' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-4"); 
    $cmd .= "'$plotf4' using 1:$field  title '$tag-$tag4' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-5"); 
    $cmd .= "'$plotf5' using 1:$field  title '$tag-$tag5' with lines ls $k, ";

    if ($GOPEN) {
	$field = 8;
	$tag = "BETA * (1-ETA) / ETA";
	$k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
	$cmd .= "'$plotf1' using 1:$field  title '$tag-$tag1' with lines ls $k, ";
	$k = FUNCS::gnuplot_set_style_transitions($tag."-2"); 
	$cmd .= "'$plotf2' using 1:$field  title '$tag-$tag2' with lines ls $k, ";
	$k = FUNCS::gnuplot_set_style_transitions($tag."-3"); 
	$cmd .= "'$plotf3' using 1:$field  title '$tag-$tag3' with lines ls $k, ";
	$k = FUNCS::gnuplot_set_style_transitions($tag."-4"); 
	$cmd .= "'$plotf4' using 1:$field  title '$tag-$tag4' with lines ls $k, ";
	$k = FUNCS::gnuplot_set_style_transitions($tag."-5"); 
	$cmd .= "'$plotf5' using 1:$field  title '$tag-$tag5' with lines ls $k, ";
    }
    
    $field = 11;
    $tag = "BETA * (1-ETA)";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
    $cmd .= "'$plotf1' using 1:$field  title '$tag-$tag1' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-2"); 
    $cmd .= "'$plotf2' using 1:$field  title '$tag-$tag2' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-3"); 
    $cmd .= "'$plotf3' using 1:$field  title '$tag-$tag3' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-4"); 
    $cmd .= "'$plotf4' using 1:$field  title '$tag-$tag4' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-5"); 
    $cmd .= "'$plotf5' using 1:$field  title '$tag-$tag5' with lines ls $k, ";
  

    # values at tstar
    $k = FUNCS::gnuplot_set_style_transitions("tstar"); 
    $field = 2;
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k, ";  
    $field = 4;
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k, "; 
    $field = 6;
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k, "; 
    $field = 8;
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";

    if ($view) {
	system("evince $Pfile&\n"); 
	system("evince $Lfile&\n"); 
    }
    
    #system("rm $t1file\n");

}

sub actual_plot_one{
    my ($name, $tstar, $plotf, $plotf_tstar, $etaz, $rz, $ld, $muE, $ldI, $muI, 
	$eta_B62, $beta_B62, $eta_B90, $beta_B90, $eta_P30, $beta_P30,  $eta_P10, $beta_P10, 
	$view) = @_;

    my $xlabel = "subtitutions per site";
    my $ylabel = "";

    my $ymin = -30;

    my $plot = $plotf; if ($plot =~ /^(\S+).plot$/) { $plot = $1; }
    my $Pfile = "$plot\_prob.ps";
    my $Lfile = "$plot\_logp.ps";

    my $tB62file = "$plot.tB62";
    my $lB62file = "$plot.lB62";

    my $tB90file = "$plot.tB90";
    my $lB90file = "$plot.lB90";
 
    my $tPAM30file = "$plot.tP30";
    my $lPAM30file = "$plot.lP30";

    my $tPAM10file = "$plot.tP10";
    my $lPAM10file = "$plot.lP10";

    open(T1, ">$tB62file") || die;
    print T1 "$tB62\t 0\n";
    print T1 "$tB62\t 1\n";
    close(T1);
    open(L1, ">$lB62file") || die;
    print L1 "$tB62\t$ymin\n";
    print L1 "$tB62\t0\n";
    close(L1);

    open(T1, ">$tB90file") || die;
    print T1 "$tB90\t 0\n";
    print T1 "$tB90\t 1\n";
    close(T1);
    open(L1, ">$lB90file") || die;
    print L1 "$tB90\t$ymin\n";
    print L1 "$tB90\t0\n";
    close(L1);
 
    open(T1, ">$tPAM30file") || die;
    print T1 "$tP30\t 0\n";
    print T1 "$tP30\t 1\n";
    close(T1);   
    open(L1, ">$lPAM30file") || die;
    print L1 "$tP30\t$ymin\n";
    print L1 "$tP30\t0\n";
    close(L1);

    open(T1, ">$tPAM10file") || die;
    print T1 "$tP10\t 0\n";
    print T1 "$tP10\t 1\n";
    close(T1);   
    open(L1, ">$lPAM10file") || die;
    print L1 "$tP10\t$ymin\n";
    print L1 "$tP10\t0\n";
    close(L1);

    my $k;
    my $eta_slope;
    my $eta_inf;
    my $beta_slope;
    my $beta_inf;
    my $kapa;

    my $eta_b62  = $eta_B62;
    my $eta_b90  = $eta_B90;   
    my $eta_p30  = $eta_P30;
    my $eta_p10  = $eta_P10;
 
    my $beta_b62 = $beta_B62;
    my $beta_b90 = $beta_B90;   
    my $beta_p30 = $beta_P30;
    my $beta_p10 = $beta_P10;

    my $gopge_b62 = $beta_B62 * (1.0-$eta_B62);
    my $gopge_b90 = $beta_B90 * (1.0-$eta_B90);
    my $gopge_p30 = $beta_P30 * (1.0-$eta_P30);
    my $gopge_p10 = $beta_P10 * (1.0-$eta_P10);

    my $logeta_b62  = ($eta_b62  > 0)? log($eta_b62)  : "-inf";
    my $logeta_b90  = ($eta_b90  > 0)? log($eta_b90)  : "-inf";
    my $logeta_p30  = ($eta_p30  > 0)? log($eta_p30)  : "-inf";
    my $logeta_p10  = ($eta_p10  > 0)? log($eta_p10)  : "-inf";

    my $logbeta_b62 = ($beta_b62 > 0)? log($beta_b62) : "-inf";
    my $logbeta_b90 = ($beta_b90 > 0)? log($beta_b90) : "-inf";
    my $logbeta_p30 = ($beta_p30 > 0)? log($beta_p30) : "-inf";
    my $logbeta_p10 = ($beta_p10 > 0)? log($beta_p10) : "-inf";

    my $loggopge_b62 = ($gopge_b62 > 0)? log($gopge_b62) : "-inf";
    my $loggopge_b90 = ($gopge_b90 > 0)? log($gopge_b90) : "-inf";
    my $loggopge_p30 = ($gopge_p30 > 0)? log($gopge_p30) : "-inf";
    my $loggopge_p10 = ($gopge_p10 > 0)? log($gopge_p10) : "-inf";

    $kapa     = $muE *(1-$rz)/($muI-$rz*$ldI);
    $eta_inf  = ($ldI<$muI)? $ldI/$muI : 1.0; 
    $beta_inf = ($ldI<$muI)? $ld/($ld+$kapa*($muI-$ldI)) : 1.0; 
    $eta_inf  = int($eta_inf  * 1000000)/1000000;
    $beta_inf = int($beta_inf * 10000000)/10000000;
    
    $eta_slope  = (1-$etaz) * ($ldI-$etaz*$muI);
    $beta_slope = $ld;
    $eta_slope  = int($eta_slope  * 10000)/10000;
    $beta_slope = int($beta_slope * 100000000)/100000000;

    $eta_b62  = int($eta_b62  * 10000)/10000;
    $eta_b90  = int($eta_b90  * 10000)/10000;
    $eta_p30  = int($eta_p30  * 10000)/10000;
    $eta_p10  = int($eta_p10  * 10000)/10000;

    $beta_b62 = int($beta_b62 * 10000)/10000;
    $beta_b90 = int($beta_b90 * 10000)/10000;
    $beta_p30 = int($beta_p30 * 10000)/10000;
    $beta_p10 = int($beta_p10 * 10000)/10000;

    $gopge_b62 = int($gopge_b62 * 10000)/10000;
    $gopge_b90 = int($gopge_b90 * 10000)/10000;
    $gopge_p30 = int($gopge_p30 * 10000)/10000;
    $gopge_p10 = int($gopge_p10 * 10000)/10000;

    $logeta_b62  = int($logeta_b62  * 100)/100;
    $logeta_b90  = int($logeta_b90  * 100)/100;
    $logeta_p30  = int($logeta_p30  * 100)/100;
    $logeta_p10  = int($logeta_p10  * 100)/100;
    
    $logbeta_b62 = int($logbeta_b62 * 100)/100;
    $logbeta_b90 = int($logbeta_b90 * 100)/100;
    $logbeta_p30 = int($logbeta_p30 * 100)/100;
    $logbeta_p10 = int($logbeta_p10 * 100)/100;

    $loggopge_b62 = int($loggopge_b62 * 100)/100;
    $loggopge_b90 = int($loggopge_b90 * 100)/100;
    $loggopge_p30 = int($loggopge_p30 * 100)/100;
    $loggopge_p10 = int($loggopge_p10 * 100)/100;
 
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    $ylabel = "Probabilities";
    print GP "set output '$Pfile'\n";
    print GP "set key right top\n";
    print GP "set xlabel '$xlabel'\n";

     
    if ($program =~ /gappenalties/) {
	if ($kappa >= 0.) {
	    print GP "set title \"$name [etaz = $etaz etainf = $etainf kappa = $kappa rz = $rz]\n";
	}
	elsif ($ld >= 0.) {
	    print GP "set title \"$name [etaz = $etaz etainf = $etainf lambda = $ld rz = $rz]\n";
	}
    }
    else {
	print GP "set title \"$name [etaz = $etaz ldI = $ldI muI = $muI ld = $ld muE = $muE rz = $rz]\n";
    }
    print GP "set ylabel '$ylabel'\n";
    if ($tmax > 0) {  print GP "set xrange  [0:$tmax]\n"; }
    print GP "set yrange [0:1]\n";
    if ($opt_g)  { 
	print GP "set logscale y\n"; 
	print GP "set yrange [$ymin:0]\n";
    }

    my $cmd;
    my $field;
    my $tag;
    
    $field = 3;
    $tag = "ETA";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
    $tag .= " ";

    $cmd   = "'$tB62file'   title 'ETA=$eta_b62 BETA=$beta_b62 GapO+GapE=$gopge_b62 || BLOSUM62' with lines ls 1, ";
    $cmd  .= "'$tB90file'   title 'ETA=$eta_b90 BETA=$beta_b90 GapO+GapE=$gopge_b90 || BLOSUM90' with lines ls 1, ";
    $cmd  .= "'$tPAM30file' title 'ETA=$eta_p30 BETA=$beta_p30 GapO+GapE=$gopge_p30 || PAM30' with lines ls 1, ";
    $cmd  .= "'$tPAM10file' title 'ETA=$eta_p10 BETA=$beta_p10 GapO+GapE=$gopge_p10 || PAM10' with lines ls 1, ";
    $tag .= " [slope $eta_slope infty-> $eta_inf]"; 
    $cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, "; 
    if ($plot_slope) {
	$field = 4;
	$k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
	$cmd .= "'$plotf' using 1:$field  with lines ls $k, "; 
    }

    $field = 6;
    $tag = "BETA";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
    $tag .= " [slope $beta_slope infty-> $beta_inf]"; 
    $cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, "; 
    if ($plot_slope) {
	$field = 7;
	$k = FUNCS::gnuplot_set_style_transitions($tag."-1-slope"); 
	$cmd .= "'$plotf' using 1:$field  with lines ls $k, "; 
   }

    if ($GOPEN) {
	$field = 9;
	$tag = "BETA * (1-ETA) / ETA";
	$k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
	$cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, ";   
	if ($plot_slope) {
	    $field = 10;
	    $k = FUNCS::gnuplot_set_style_transitions($tag."-1-slope"); 
	    $cmd .= "'$plotf' using 1:$field  with lines ls $k, ";
	}
    }

    $field = 12;
    $tag = "BETA * (1-ETA)";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
    $cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, ";   
    if ($plot_slope) {
	$field = 13;
	$k = FUNCS::gnuplot_set_style_transitions($tag."-1-slope"); 
	$cmd .= "'$plotf' using 1:$field  with lines ls $k, ";
    }

    # at time tstar
    $k = FUNCS::gnuplot_set_style_transitions("tstar"); 
    $field = 3;
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k, "; 
    $field = 5;
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k, "; 
    $field = 7;
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k", ; 
    $field = 9;
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k"; 

    $cmd .= "\n";
    print GP "plot $cmd\n";


    $ylabel = "Log Probabilities";
    print GP "set output '$Lfile'\n";
    $cmd   = "'$lB62file'   title 'logETA=$logeta_b62 logBETA=$logbeta_b62 GapO+GapE=$loggopge_b62 || BLOSUM62' with lines ls 1, ";
    $cmd  .= "'$lB90file'   title 'logETA=$logeta_b90 logBETA=$logbeta_b90 GapO+GapE=$loggopge_b90 || BLOSUM90' with lines ls 1, ";
    $cmd  .= "'$lPAM30file' title 'logETA=$logeta_p30 logBETA=$logbeta_p30 GapO+GapE=$loggopge_p30 || PAM30' with lines ls 1, ";
    $cmd  .= "'$lPAM10file' title 'logETA=$logeta_p10 logBETA=$logbeta_p10 GapO+GapE=$loggopge_p10 || PAM10' with lines ls 1, ";
    print GP "set yrange [$ymin:0]\n";
    print GP "set ylabel '$ylabel'\n";
    print GP "set key bottom right\n";

    $field = 2;
    $tag = "ETA";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
    $tag .= " ";
    $tag .= " [slope $eta_slope infty-> $eta_inf]"; 
    $cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, ";
 
    $field = 5;
    $tag = "BETA";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
    $tag .= " [slope $beta_slope infty-> $beta_inf]"; 
    $cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, ";
 
    if ($GOPEN) {
	$field = 8;
	$tag = "BETA * (1-ETA) / ETA";
	$k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
	$cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, ";
    }

    $field = 11;
    $tag = "BETA * (1-ETA)";
    $k = FUNCS::gnuplot_set_style_transitions($tag."-1"); 
    $cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, ";
 
    # at time tstar
    $k = FUNCS::gnuplot_set_style_transitions("tstar"); 
    $field = 2;
    $tag = "ETA(gap_extend=-1)";
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k, ";  
    $field = 4;
    $tag = "BETA";
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k, "; 
    if ($GOPEN) {
	$field = 6;
	$tag = "BETA * (1-ETA) / ETA (gap_open=-11)";
	$cmd .= "'$plotf_tstar' using 1:$field title '' ls $k, "; 
    }
    $field = 8;
    $tag = "BETA * (1-ETA) (gap_Open=-12)";
    $cmd .= "'$plotf_tstar' using 1:$field title '' ls $k"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";

    if ($view) { 
	system("evince $Pfile&\n"); 
	system("evince $Lfile\n"); 
    }
    #system("rm $t1file\n");

}

sub create_plots {
    my ($wrkdir, $name, $type, $file1, $file2, $file3, $file4, $file5, $view) = @_;
    
    my $plotf1      = "$wrkdir/$name-$type.plot1";
    my $plotf2      = "$wrkdir/$name-$type.plot2";
    my $plotf3      = "$wrkdir/$name-$type.plot3";
    my $plotf4      = "$wrkdir/$name-$type.plot4";
    my $plotf5      = "$wrkdir/$name-$type.plot5";
    my $plotf_tstar = "$wrkdir/$name-$type.plot_tstar";
    
    my $tstar1 = parse_outfile($plotf1, $file1, 0);
    my $tstar2 = parse_outfile($plotf2, $file2, 0);
    my $tstar3 = parse_outfile($plotf3, $file3, 0);
    my $tstar4 = parse_outfile($plotf4, $file4, 0); 
    my $tstar5 = parse_outfile($plotf5, $file5, 0); 	

    my $tstar = parse_outfile($plotf_tstar, $file5, 1); 
       
    if (abs($tstar1-$tstar) > 1e-6 ||
	abs($tstar2-$tstar) > 1e-6 || 
	abs($tstar3-$tstar) > 1e-6 || 
	abs($tstar4-$tstar) > 1e-6 || 
	abs($tstar5-$tstar) > 1e-6  ) { print "bad tstar\n"; die; }

    actual_plot($name, $type, $tstar, $plotf1, $plotf2, $plotf3, $plotf4, $plotf5, $plotf_tstar, $view);
}

sub create_plot_one {
    my ($wrkdir, $name, $file, $view) = @_;

    my $plotf       = "$wrkdir/$name.plot";
    my $plotf_tstar = "$wrkdir/$name.plot_tstar";
    
    my $etaz;
    my $rz;
    my $ld;
    my $muE;
    my $ldI;
    my $muI;

    my $eta_B62 = "";
    my $eta_B90 = "";
    my $eta_P30 = "";
    my $eta_P10 = "";
    
    my $beta_B62 = "";
    my $beta_B90 = "";
    my $beta_P30 = "";
    my $beta_P10 = "";

    my $tstar1 = parse_outfile($plotf,       $file, 0, \$etaz, \$rz, \$ld, \$muE, \$ldI, \$muI, \$eta_B62, \$beta_B62, \$eta_B90, \$beta_B90, \$eta_P30, \$beta_P30, \$eta_P10, \$beta_P10); 
    my $tstar  = parse_outfile($plotf_tstar, $file, 1); 
    if (abs($tstar1-$tstar) > 1e-6) { print "bad tstar\n"; die; }
    if ($tstar) { actual_plot_one($name, $tstar, $plotf,  $plotf_tstar, $etaz, $rz, $ld, $muE, $ldI, $muI, 
				  $eta_B62, $beta_B62, $eta_B90, $beta_B90, $eta_P30, $beta_P30, $eta_P10, $beta_P10, $view); }
}

sub parse_outfile {

    my ($plotfile, $file, $ist1, $ret_etaz, $ret_rz, $ret_ld, $ret_muE, $ret_ldI, $ret_muI, 
	$ret_eta_B62, $ret_beta_B62, $ret_eta_B90, $ret_beta_B90,
	$ret_eta_P30, $ret_beta_P30, $ret_eta_P10, $ret_beta_P10
	) = @_;

    open(OUT, ">$plotfile") || die;

    my $tstar;
    my $etaz;
    my $rz;

    my $muI;
    my $muE;
    my $ldI;
    my $ld;
    my $eta_B62;
    my $eta_B90;
    my $eta_P30;
    my $eta_P10;
    my $beta_B62;
    my $beta_B90;
    my $beta_P30;
    my $beta_P10;
 
    my $time    = -1;
    my $subsite = 0;
    my $eta     = 0;
    my $beta    = 0;
    my $gapo    = 0;
    my $gapO    = 0;
    my $logeta  = 0;
    my $logbeta = 0;
    my $loggapo = 0;
    my $loggapO = 0;

    my $eta_slope  = 0;
    my $beta_slope = 0;
    my $gapo_slope = 0;
    my $gapO_slope = 0;

    my $eta_linear  = 0;
    my $beta_linear = 0;
    my $gapo_linear = 0;
    my $gapO_linear = 0;
 
    open (FILE, "$file") || die;
    while (<FILE>) {
	
	if (/^# eta at time zero:\s+=\s+(\S+)\s+/) {
	    $etaz = $1;
	}
	elsif (/^\# r0\s+at time zero:\s+=\s+(\S+)\s+/) {
	    $rz = $1;
	}
	elsif (/^\#\s*ld\s+=\s+(\S+)/) { 
	    $ld = $1;
	}
	elsif (/^\#\s*muE\s+=\s+(\S+)/) { 
	    $muE = $1;
	}
	elsif (/^\#\s*ldI\s+=\s+(\S+)/) { 
	    $ldI = $1;
	}
	elsif (/^\#\s*muI\s+=\s+(\S+)/) { 
	    $muI = $1;
	}
	elsif (/^\#/) { next; 
	}
	elsif (/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)$/) {
	    $time     = $1;
	    $subsite  = $2;
	    $logeta   = $3;
	    $logbeta  = $4;
	    $loggapo  = $5;
	    $loggapO  = $6;
	    if ($time == 1.0) { $tstar = $subsite; }

	    $eta  = exp($logeta);
	    $beta = exp($logbeta);
	    $gapo = exp($loggapo);
	    $gapO = exp($loggapO);
	    
	    # the slopes
	    $eta_slope  = (1.-$etaz) *($ldI - $etaz*$muI);
	    $beta_slope = $ld;
	    $gapo_slope = ($etaz > 0.)? $ld * (1.-$etaz) / $etaz : 0.;
	    $gapO_slope = $ld * (1.-$etaz);

	    if ($eta_slope  < 0.) { print "bad slope for eta  $eta_slope\n";  die; }
	    if ($beta_slope < 0.) { print "bad slope for beta $beta_slope\n"; die; }
	    if ($gapo_slope < 0.) { print "bad slope for gapo $gapo_slope\n"; die; }
	    if ($gapO_slope < 0.) { print "bad slope for gapO $gapO_slope\n"; die; }

	    $eta_linear  = $time * $eta_slope + $etaz;
	    $beta_linear = $time * $beta_slope;
	    $gapo_linear = $time * $gapo_slope;
	    $gapO_linear = $time * $gapO_slope;

	    if ($subsite <= $tB62) { $eta_B62 = $eta; $beta_B62 = $beta; }
	    if ($subsite <= $tB90) { $eta_B90 = $eta; $beta_B90 = $beta; }
	    if ($subsite <= $tP30) { $eta_P30 = $eta; $beta_P30 = $beta; }
	    if ($subsite <= $tP10) { $eta_P10 = $eta; $beta_P10 = $beta; }

	    if (!$ist1) {
		#$logeta = ($beta<1&&$eta>0)?log($eta) -log(1-$beta):-1e-20;
		#$loggapO = ($beta<1&& $eta<1&&$beta>0)?log($beta)+log(1-$eta)-2*log(1-$beta):-1e-20;
		print OUT "$subsite\t$logeta\t$eta\t$eta_linear\t$logbeta\t$beta\t$beta_linear\t$loggapo\t$gapo\t$gapo_linear\t$loggapO\t$gapO\t$gapO_linear\n";
	    }
	    else {
		if ($time == 1.0) {
		    print     "$subsite\t$logeta\t$eta\t$logbeta\t$beta\t$loggapo\t$gapo\t$loggapO\t$gapO\n";
		    print OUT "$subsite\t$logeta\t$eta\t$logbeta\t$beta\t$loggapo\t$gapo\t$loggapO\t$gapO\n";
		}
	    }
	}
   }
    close (FILE);
    close (OUT);

    $$ret_etaz = $etaz;
    $$ret_rz   = $rz;
    $$ret_ld   = $ld;
    $$ret_muE  = $muE;
    $$ret_ldI  = $ldI;
    $$ret_muI  = $muI;
    $$ret_eta_B62  = $eta_B62;
    $$ret_eta_B90  = $eta_B90;
    $$ret_eta_P30  = $eta_P30;
    $$ret_eta_P10  = $eta_P10;
    $$ret_beta_B62 = $beta_B62;
    $$ret_beta_B90 = $beta_B90;
    $$ret_beta_P30 = $beta_P30;
    $$ret_beta_P10 = $beta_P10;

    return $tstar;
}

sub plot_program {
    my ($wrkdir,  $which, $view) = @_;

    my $name = $wrkdir;
    if ($name =~ /^\S+\/([^\/]+)$/) { $name = $1; }

    my $etaz_def      = 0.01; 
    my $rz_def        = 0.01; 
    my $etainf_def    = 0.99;
    my $kappa_def     = 1.0;

    my $etaz; 
    my $rz; 
    my $etainf;

    my $file1;
    my $file2;
    my $file3;
    my $file4;
    my $file5;

    if ($which == 1) {
	$rz      = $rz_def; 
	$etainf  = $etainf_def;
	$kappa = $kappa_def;
	
	$file1 = "$wrkdir/$name\_etaz$etaz1\_rz$rz\_etainf$etainf\_kappa$kappa.out";
	$file2 = "$wrkdir/$name\_etaz$etaz2\_rz$rz\_etainf$etainf\_kappa$kappa.out";
	$file3 = "$wrkdir/$name\_etaz$etaz3\_rz$rz\_etainf$etainf\_kappa$kappa.out";
	$file4 = "$wrkdir/$name\_etaz$etaz4\_rz$rz\_etainf$etainf\_kappa$kappa.out";
	$file5 = "$wrkdir/$name\_etaz$etaz5\_rz$rz\_etainf$etainf\_kappa$kappa.out";
	
	create_plots($wrkdir, $name, "etaz", $file1, $file2, $file3, $file4, $file5, $view);
    }
    elsif ($which == 2) {
	$etaz      = $etaz_def; 
	$etainf  = $etainf_def;
	$kappa = $kappa_def;
 	
	$file1 = "$wrkdir/$name\_etaz$etaz\_rz$rz1\_etainf$etainf\_kappa$kappa.out";
	$file2 = "$wrkdir/$name\_etaz$etaz\_rz$rz2\_etainf$etainf\_kappa$kappa.out";
	$file3 = "$wrkdir/$name\_etaz$etaz\_rz$rz3\_etainf$etainf\_kappa$kappa.out";
	$file4 = "$wrkdir/$name\_etaz$etaz\_rz$rz4\_etainf$etainf\_kappa$kappa.out";
	$file5 = "$wrkdir/$name\_etaz$etaz\_rz$rz5\_etainf$etainf\_kappa$kappa.out";
	
	create_plots($wrkdir, $name, "rz", $file1, $file2, $file3, $file4, $file5, $view);
    }
     elsif ($which == 3) {
	$etaz      = $etaz_def; 
	$rz        = $rz_def; 
	$kappa = $kappa_def;
	
	$file1 = "$wrkdir/$name\_etaz$etaz\_rz$rz\_etainf$etainf1\_kappa$kappa.out";
	$file2 = "$wrkdir/$name\_etaz$etaz\_rz$rz\_etainf$etainf2\_kappa$kappa.out";
	$file3 = "$wrkdir/$name\_etaz$etaz\_rz$rz\_etainf$etainf3\_kappa$kappa.out";
	$file4 = "$wrkdir/$name\_etaz$etaz\_rz$rz\_etainf$etainf4\_kappa$kappa.out";
	$file5 = "$wrkdir/$name\_etaz$etaz\_rz$rz\_etainf$etainf5\_kappa$kappa.out";
	
	create_plots($wrkdir, $name, "etainf", $file1, $file2, $file3, $file4, $file5, $view);
    }
     elsif ($which == 4) {
	$etaz     = $etaz_def; 
	$rz       = $rz_def; 
	$etainf = $etainf_def;
	
	$file1 = "$wrkdir/$name\_etaz$etaz\_rz$rz\_etainf$etainf\_kappa$kappa1.out";
	$file2 = "$wrkdir/$name\_etaz$etaz\_rz$rz\_etainf$etainf\_kappa$kappa2.out";
	$file3 = "$wrkdir/$name\_etaz$etaz\_rz$rz\_etainf$etainf\_kappa$kappa3.out";
	$file4 = "$wrkdir/$name\_etaz$etaz\_rz$rz\_etainf$etainf\_kappa$kappa4.out";
	$file5 = "$wrkdir/$name\_etaz$etaz\_rz$rz\_etainf$etainf\_kappa$kappa5.out";
	
	create_plots($wrkdir, $name, "etainf", $file1, $file2, $file3, $file4, $file5, $view);
    }

}

sub plot_program_one {
    my ($wrkdir, $view) = @_;

    my $name = $wrkdir;
    if ($name =~ /^\S+\/([^\/]+)$/) { $name = $1; }

    my $file = "$wrkdir/$name.out";

    create_plot_one($wrkdir, $name, $file, $view);
}



sub run_gappenalties {
    my ($program, $wrkdir, $filename, $etaz, $rz, $etainf, $kappa, $ld, $verbose) = @_;
    
    my $outfile = "$wrkdir/$filename.out";
    
    my $cmd;
    if    ($kappa >= 0.0) { $cmd = "$program --etaz $etaz --rz $rz --etainf $etainf --kappa  $kappa "; }
    elsif ($ld    >= 0.0) { $cmd = "$program --etaz $etaz --rz $rz --etainf $etainf --lambda $ld ";    }
    
    if ($verbose) { $cmd .= " -v "; }
    system("echo $cmd $outfile\n");
    system("$cmd $outfile\n");

    return 1;
}

sub run_evoprobs {
    my ($program, $wrkdir, $filename, $etaz, $rz, $ld, $muE, $ldI, $muI, $verbose) = @_;
    
    my $outfile = "$wrkdir/$filename.out";
    
    my $cmd = "$program --etaz $etaz --rz $rz --ld $ld --ldI $ldI --muE $muE --muI $muI  ";  
    
    if ($verbose) { $cmd .= " -v "; }
    system("echo $cmd $outfile\n");
    system("$cmd $outfile\n");

    return 1;
}
