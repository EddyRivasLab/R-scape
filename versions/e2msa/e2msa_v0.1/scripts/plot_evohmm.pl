#!/usr/bin/perl -w
#plot_evohmm.pl

use strict;
use Class::Struct;
use lib '/groups/eddy/home/rivase/projects/evohmm/scripts';
use lib '/Users/rivase/projects/evohmm/scripts';
use FUNCS;
use constant GNUPLOT => '/usr/bin/gnuplot';

use vars qw ($opt_I $opt_l $opt_r $opt_s $opt_t $opt_v $opt_x);  # required if strict used
use Getopt::Std;
getopts ('I:lr:st:vx:');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  plot_evohmm.pl [options] <result_dir> <hmm>\n\n";
        print "options:\n";
        print "-s    :  don't plot slopes       \n";
        print "-r    :  rz value       \n";
	print "-x    :  xiz value      \n";
	print "-I    :  xiIinfty value \n";
	print "-l    :  log scale for probabilities \n";
	print "-v    :  be verbose\n";
        print "-t    :  tmax [default is 20]\n";
	exit;
}
my $wrkdir   = shift;
my $hmmfile  = shift;

if ( -d $wrkdir) { die "wrkdir directory $wrkdir exist!"; }
system("mkdir $wrkdir\n");

my $hmmname = $hmmfile;
if ($hmmname =~ /^\S+\/([^\/]+)\.hmm$/) { $hmmname = $1; }
else {print "bad hmm file $hmmfile\n"; die; }

my $evohmm_utest = "/groups/eddy/home/rivase/projects/evohmm/src/evohmmer_utest";

my $plot_slope = 1;
if ($opt_s) { $plot_slope = 0; }

my $plot_default_only = 0;

my $tmax = 20; if ($opt_t) { $tmax = $opt_t; if ($tmax < 0) { print "tmax has to be positive\n"; die; } }

# variables
my $xiz      = 0.01; 
my $rz       = 0.01; 
my $xiIinfty = 0.99;

if ($opt_x) { $xiz      = $opt_x; $plot_default_only = 1; if ( $xiz < 0 || $xiz > 1.0)            { print "bad value of xiz=$xiz\n";           die; } }
if ($opt_r) { $rz       = $opt_r; $plot_default_only = 1; if ( $rz  < 0 || $rz  > 1.0)            { print "bad value of rz=$rz\n";             die; } }
if ($opt_I) { $xiIinfty = $opt_I; $plot_default_only = 1; if ( $xiIinfty  < 0 || $xiIinfty > 1.0) { print "bad value of xiIinfty=$xiIinfty\n"; die; } }

my $xiz1 = 0.0001;
my $xiz2 = 0.001;
my $xiz3 = 0.10;
my $xiz4 = 0.15;
my $xiz5 = 0.20;

my $rz1 = 0.01;
my $rz2 = 0.10;
my $rz3 = 0.50;
my $rz4 = 0.90;
my $rz5 = 0.99;

my $xiIinfty5 = 0.50;
my $xiIinfty4 = 0.60;
my $xiIinfty3 = 0.80;
my $xiIinfty2 = 0.90;
my $xiIinfty1 = 0.99;

if ($plot_default_only) {
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz, $rz, $xiIinfty);
    plot_evohmm_utest_one($wrkdir, $xiz, $rz, $xiIinfty); 
}
else {
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz1, $rz, $xiIinfty);
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz2, $rz, $xiIinfty);
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz3, $rz, $xiIinfty);
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz4, $rz, $xiIinfty);
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz5, $rz, $xiIinfty);
    
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz, $rz1, $xiIinfty);
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz, $rz2, $xiIinfty);
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz, $rz3, $xiIinfty);
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz, $rz4, $xiIinfty);
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz, $rz5, $xiIinfty);
    
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz, $rz, $xiIinfty1);
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz, $rz, $xiIinfty2);
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz, $rz, $xiIinfty3);
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz, $rz, $xiIinfty4);
    run_evohmm_utest($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz, $rz, $xiIinfty5);
    
    
    plot_evohmm_utest($wrkdir, 1); # xiz      variations 
    plot_evohmm_utest($wrkdir, 2); # rz       variations 
    plot_evohmm_utest($wrkdir, 3); # xiIinfty variations 
}

sub actual_plot {
    my ($name, $type, $m, $plotf1, $plotf2, $plotf3, $plotf4, $plotf5, $plotf_t1) = @_;

    my $xlabel = "time";
    my $ylabel = "Transition Probabilities";

    my $plotf = $plotf1; if ($plotf =~ /^(\S+).plot1$/) { $plotf = $1; }
    my $Mfile  = "$plotf\_M.ps";
    my $DIfile = "$plotf\_DI.ps";

    my $t1file = "$plotf.t1";

    open(T1, ">$t1file") || die;
    print T1 "1\t 0\n";
    print T1 "1\t 1\n";
    close(T1);

    my $k;

    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    print GP "set output '$Mfile'\n";
    print GP "set key right top\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set title \"$name position $m vary $type\n";
    print GP "set ylabel '$ylabel'\n";
    print GP "set xrange [0:$tmax]\n";
    print GP "set yrange [0:1]\n";

    my $tag1;
    my $tag2;
    my $tag3;
    my $tag4;
    my $tag5;
    if ($type =~ /^xiz$/) {
	$tag1 = $xiz1;
	$tag2 = $xiz2;
	$tag3 = $xiz3;
	$tag4 = $xiz4;
	$tag5 = $xiz5;
    }
    elsif ($type =~ /^rz$/) {
	$tag1 = $rz1;
	$tag2 = $rz2;
	$tag3 = $rz3;
	$tag4 = $rz4;
	$tag5 = $rz5;
    }
   elsif ($type =~ /^xiIinfty$/) {
	$tag1 = $xiIinfty1;
	$tag2 = $xiIinfty2;
	$tag3 = $xiIinfty3;
	$tag4 = $xiIinfty4;
	$tag5 = $xiIinfty5;
    }

    my $cmd;
    my $field;
    
    $field = 2;
    $cmd  = "'$t1file' title 'trained HMM' with lines ls 1, ";
    $k = FUNCS::gnuplot_set_style("tMM-1"); 
    $cmd .= "'$plotf1' using 1:$field  title 'tMM-$tag1' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style("tMM-2"); 
    $cmd .= "'$plotf2' using 1:$field  title 'tMM-$tag2' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style("tMM-3"); 
    $cmd .= "'$plotf3' using 1:$field  title 'tMM-$tag3' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style("tMM-4"); 
    $cmd .= "'$plotf4' using 1:$field  title 'tMM-$tag4' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style("tMM-5"); 
    $cmd .= "'$plotf5' using 1:$field  title 'tMM-$tag5' with lines ls $k, "; 

    if ($plot_slope) {
	$field = 3;
	$k = FUNCS::gnuplot_set_style("tMM-1-slope"); 
	$cmd .= "'$plotf1' using 1:$field  with lines ls $k, "; 
	if (0) {
	    $k = FUNCS::gnuplot_set_style("tMM-2-slope"); 
	    $cmd .= "'$plotf2' using 1:$field  with lines ls $k, "; 
	    $k = FUNCS::gnuplot_set_style("tMM-3-slope"); 
	    $cmd .= "'$plotf3' using 1:$field  with lines ls $k, "; 
	    $k = FUNCS::gnuplot_set_style("tMM-4-slope"); 
	    $cmd .= "'$plotf4' using 1:$field   with lines ls $k, "; 
	    $k = FUNCS::gnuplot_set_style("tMM-5-slope"); 
	    $cmd .= "'$plotf5' using 1:$field  with lines ls $k, "; 
	}
    }

    $field = 4;
    $k = FUNCS::gnuplot_set_style("tMI-1"); 
    $cmd .= "'$plotf1' using 1:$field  title 'tMI-$tag1' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style("tMI-2"); 
    $cmd .= "'$plotf2' using 1:$field  title 'tMI-$tag2' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style("tMI-3"); 
    $cmd .= "'$plotf3' using 1:$field  title 'tMI-$tag3' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style("tMI-4"); 
    $cmd .= "'$plotf4' using 1:$field  title 'tMI-$tag4' with lines ls $k, "; 
    $k = FUNCS::gnuplot_set_style("tMI-5"); 
    $cmd .= "'$plotf5' using 1:$field  title 'tMI-$tag5' with lines ls $k, "; 

    if ($plot_slope) {
	$field = 5;
	$k = FUNCS::gnuplot_set_style("tMI-1-slope"); 
	$cmd .= "'$plotf1' using 1:$field  with lines ls $k, "; 
	if (0) {
	    $k = FUNCS::gnuplot_set_style("tMI-2-slope"); 
	    $cmd .= "'$plotf2' using 1:$field  with lines ls $k, "; 
	    $k = FUNCS::gnuplot_set_style("tMI-3-slope"); 
	    $cmd .= "'$plotf3' using 1:$field  with lines ls $k, "; 
	    $k = FUNCS::gnuplot_set_style("tMI-4-slope"); 
	    $cmd .= "'$plotf4' using 1:$field  with lines ls $k, ";    
	    $k = FUNCS::gnuplot_set_style("tMI-5-slope"); 
	    $cmd .= "'$plotf5' using 1:$field  with lines ls $k, "; 
	}
    }

    $field = 6;
    $k = FUNCS::gnuplot_set_style("tMD-1"); 
    $cmd .= "'$plotf1' using 1:$field  title 'tMD-$tag1' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style("tMD-2"); 
    $cmd .= "'$plotf2' using 1:$field  title 'tMD-$tag2' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style("tMD-3"); 
    $cmd .= "'$plotf3' using 1:$field  title 'tMD-$tag3' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style("tMD-4"); 
    $cmd .= "'$plotf4' using 1:$field  title 'tMD-$tag4' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style("tMD-5"); 
    $cmd .= "'$plotf5' using 1:$field  title 'tMD-$tag5' with lines ls $k, ";
    
    if ($plot_slope) {
	$field = 7;
	$k = FUNCS::gnuplot_set_style("tMD-1-slope"); 
	$cmd .= "'$plotf1' using 1:$field  with lines ls $k, ";
	if (0) {
	    $k = FUNCS::gnuplot_set_style("tMD-2-slope"); 
	    $cmd .= "'$plotf2' using 1:$field  with lines ls $k, ";
	    $k = FUNCS::gnuplot_set_style("tMD-3-slope"); 
	    $cmd .= "'$plotf3' using 1:$field  with lines ls $k, ";
	    $k = FUNCS::gnuplot_set_style("tMD-4-slope"); 
	    $cmd .= "'$plotf4' using 1:$field  with lines ls $k, ";   
	    $k = FUNCS::gnuplot_set_style("tMD-5-slope"); 
	    $cmd .= "'$plotf5' using 1:$field  with lines ls $k, ";
	}
    }

    $field = 2;
    $k = FUNCS::gnuplot_set_style("tstar"); 
    $cmd .= "'$plotf_t1' using 1:$field  ls $k, "; 
    $field = 3;
    $cmd .= "'$plotf_t1' using 1:$field  ls $k, "; 
    $field = 4;
    $cmd .= "'$plotf_t1' using 1:$field  ls $k"; 

    $cmd .= "\n";
    print GP "plot $cmd\n";

    print GP "set output '$DIfile'\n";
    $cmd  = "'$t1file' title 'trained HMM' with lines ls 1, ";

    $field = 8;
    $k = FUNCS::gnuplot_set_style("tII-1"); 
    $cmd .= "'$plotf1' using 1:$field  title 'tII-$tag1' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style("tII-2"); 
    $cmd .= "'$plotf2' using 1:$field  title 'tII-$tag2' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style("tII-3"); 
    $cmd .= "'$plotf3' using 1:$field  title 'tII-$tag3' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style("tII-4"); 
    $cmd .= "'$plotf4' using 1:$field  title 'tII-$tag4' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style("tII-5"); 
    $cmd .= "'$plotf5' using 1:$field  title 'tII-$tag5' with lines ls $k, ";
 
    if ($plot_slope) {
	$field = 9;
	$k = FUNCS::gnuplot_set_style("tII-1-slope"); 
	$cmd .= "'$plotf1' using 1:$field with lines ls $k, ";
	if (0) {
	    $k = FUNCS::gnuplot_set_style("tII-2-slope"); 
	    $cmd .= "'$plotf2' using 1:$field with lines ls $k, ";
	    $k = FUNCS::gnuplot_set_style("tII-3-slope"); 
	    $cmd .= "'$plotf3' using 1:$field  with lines ls $k, ";
	    $k = FUNCS::gnuplot_set_style("tII-4-slope"); 
	    $cmd .= "'$plotf4' using 1:$field  with lines ls $k, ";
	    $k = FUNCS::gnuplot_set_style("tII-5-slope"); 
	    $cmd .= "'$plotf5' using 1:$field  with lines ls $k, ";
	}
    }

    $field = 10;
    $k = FUNCS::gnuplot_set_style("tDD-1"); 
    $cmd .= "'$plotf1' using 1:$field  title 'tDD-$tag1' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style("tDD-2"); 
    $cmd .= "'$plotf2' using 1:$field  title 'tDD-$tag2' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style("tDD-3"); 
    $cmd .= "'$plotf3' using 1:$field  title 'tDD-$tag3' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style("tDD-4"); 
    $cmd .= "'$plotf4' using 1:$field  title 'tDD-$tag4' with lines ls $k, ";
    $k = FUNCS::gnuplot_set_style("tDD-5"); 
    $cmd .= "'$plotf5' using 1:$field  title 'tDD-$tag5' with lines ls $k, ";
 
    if ($plot_slope) {
	$field = 11;
	$k = FUNCS::gnuplot_set_style("tDD-1-slope"); 
	$cmd .= "'$plotf1' using 1:$field  with lines ls $k, ";
	if (0) {
	    $k = FUNCS::gnuplot_set_style("tDD-2-slope"); 
	    $cmd .= "'$plotf2' using 1:$field  with lines ls $k, ";
	    $k = FUNCS::gnuplot_set_style("tDD-3-slope"); 
	    $cmd .= "'$plotf3' using 1:$field  with lines ls $k, ";
	    $k = FUNCS::gnuplot_set_style("tDD-4-slope"); 
	    $cmd .= "'$plotf4' using 1:$field  with lines ls $k, ";
	    $k = FUNCS::gnuplot_set_style("tDD-5-slope"); 
	    $cmd .= "'$plotf5' using 1:$field  with lines ls $k, ";
	}
    }

    $field = 5;
    $k = FUNCS::gnuplot_set_style("tstar"); 
    $cmd .= "'$plotf_t1' using 1:$field ls $k, ";  
    $field = 6;
    $cmd .= "'$plotf_t1' using 1:$field ls $k"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";

    #system("rm $t1file\n");

}

sub actual_plot_one{
    my ($name, $m, $plotf, $plotf_t1) = @_;

    my $xlabel = "time";
    my $ylabel = "Transition Probabilities";

    my $plot = $plotf; if ($plot =~ /^(\S+).plot$/) { $plot = $1; }
    my $Mfile  = "$plot\_M.ps";
    my $DIfile = "$plot\_DI.ps";

    my $t1file = "$plot.t1";

    open(T1, ">$t1file") || die;
    print T1 "1\t 0\n";
    print T1 "1\t 1\n";
    close(T1);

    my $k;

    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    print GP "set output '$Mfile'\n";
    print GP "set key right top\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set title \"$name position $m\n";
    print GP "set ylabel '$ylabel'\n";
    print GP "set xrange  [0:$tmax]\n";
    print GP "set yrange [0:1]\n";
    print GP "set nokey\n";
    if ($opt_l)  { 
	print GP "set logscale y\n"; 
	print GP "set yrange [-20:0]\n";
    }

    my $cmd;
    my $field;
    
    $field = 2;
    $cmd  = "'$t1file' title 'trained HMM' with lines ls 1, ";
    $k = FUNCS::gnuplot_set_style("tMM-1"); 
    $cmd .= "'$plotf' using 1:$field  title 'tMM' with lines ls $k, "; 

    if ($plot_slope) {
	$field = 3;
	$k = FUNCS::gnuplot_set_style("tMM-1"); 
	$cmd .= "'$plotf' using 1:$field  with lines ls $k, "; 
    }

    $field = 4;
    $k = FUNCS::gnuplot_set_style("tMI-1"); 
    $cmd .= "'$plotf' using 1:$field  title 'tMI' with lines ls $k, "; 

    if ($plot_slope) {
	$field = 5;
	$k = FUNCS::gnuplot_set_style("tMI-1-slope"); 
	$cmd .= "'$plotf' using 1:$field  with lines ls $k, "; 
   }

    $field = 6;
    $k = FUNCS::gnuplot_set_style("tMD-1"); 
    $cmd .= "'$plotf' using 1:$field  title 'tMD' with lines ls $k, ";
   
    if ($plot_slope) {
	$field = 7;
	$k = FUNCS::gnuplot_set_style("tMD-1-slope"); 
	$cmd .= "'$plotf' using 1:$field  with lines ls $k, ";
    }

    $field = 2;
    $k = FUNCS::gnuplot_set_style("tstar"); 
    $cmd .= "'$plotf_t1' using 1:$field  ls $k, "; 
    $field = 3;
    $cmd .= "'$plotf_t1' using 1:$field  ls $k, "; 
    $field = 4;
    $cmd .= "'$plotf_t1' using 1:$field  ls $k"; 

    $cmd .= "\n";
    print GP "plot $cmd\n";

    print GP "set output '$DIfile'\n";
    $cmd  = "'$t1file' title 'trained HMM' with lines ls 1, ";

    $field = 8;
    $k = FUNCS::gnuplot_set_style("tII-1"); 
    $cmd .= "'$plotf' using 1:$field  title 'tII' with lines ls $k, ";
 
    if ($plot_slope) {
	$field = 9;
	$k = FUNCS::gnuplot_set_style("tII-1-slope"); 
	$cmd .= "'$plotf' using 1:$field with lines ls $k, ";
    }

    $field = 10;
    $k = FUNCS::gnuplot_set_style("tDD-1"); 
    $cmd .= "'$plotf' using 1:$field  title 'tDD' with lines ls $k, ";
 
    if ($plot_slope) {
	$field = 11;
	$k = FUNCS::gnuplot_set_style("tDD-1-slope"); 
	$cmd .= "'$plotf' using 1:$field  with lines ls $k, ";
   }

    $field = 5;
    $k = FUNCS::gnuplot_set_style("tstar"); 
    $cmd .= "'$plotf_t1' using 1:$field ls $k, ";  
    $field = 6;
    $cmd .= "'$plotf_t1' using 1:$field ls $k"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";

    #system("rm $t1file\n");

}

sub create_plots {
    my ($wrkdir, $name, $type, $file1, $file2, $file3, $file4, $file5) = @_;

    my $M = get_M($file1);
 
    my $c = 0;
    while ($c < $M) {	
	my $m = $c+1;
	my $plotf1   = "$wrkdir/$type\_m$m.plot1";
	my $plotf2   = "$wrkdir/$type\_m$m.plot2";
	my $plotf3   = "$wrkdir/$type\_m$m.plot3";
	my $plotf4   = "$wrkdir/$type\_m$m.plot4";
	my $plotf5   = "$wrkdir/$type\_m$m.plot5";
	my $plotf_t1 = "$wrkdir/$type\_m$m.plot_t1";

	my $M1 = parse_outfile($m, $plotf1, $file1, 0);
	my $M2 = parse_outfile($m, $plotf2, $file2, 0); if ($M2 != $M) { print "bad file, M is $M2 but it should be $M\n"; die; }
	my $M3 = parse_outfile($m, $plotf3, $file3, 0); if ($M3 != $M) { print "bad file, M is $M3 but it should be $M\n"; die; }
	my $M4 = parse_outfile($m, $plotf4, $file4, 0); if ($M4 != $M) { print "bad file, M is $M4 but it should be $M\n"; die; }
	my $M5 = parse_outfile($m, $plotf5, $file5, 0); if ($M5 != $M) { print "bad file, M is $M5 but it should be $M\n"; die; }	

	parse_outfile($m, $plotf_t1, $file5, 1); 	

	actual_plot($name, $type, $m, $plotf1, $plotf2, $plotf3, $plotf4, $plotf5, $plotf_t1);
	$c ++;
    }
}

sub create_plot_one {
    my ($wrkdir, $name, $file) = @_;

    my $M = get_M($file);
 
    my $c = 0;
    while ($c < $M) {	
	my $m = $c+1;
	my $plotf   = "$wrkdir/m$m.plot";
	my $plotf_t1 = "$wrkdir/m$m.plot_t1";

	my $M1 = parse_outfile($m, $plotf, $file, 0); if ($M != $M) { print "bad file, M is $M1 but it should be $M\n"; die; }

	parse_outfile($m, $plotf_t1, $file, 1); 	

	actual_plot_one($name, $m, $plotf, $plotf_t1);
	$c ++;
    }
}


sub get_M {

    my ($file) = @_;

    my $M = 0;

    open (FILE, "$file") || die; 
    while (<FILE>) {
	if (/^M=(\d+)\s+muD/) { $M = $1; }
    }
    close (FILE);

    if ($M <= 0)  { print "bad file. M should be positive, but it is $M in file $file\n"; die; }
    return $M;
}

sub parse_outfile {

    my ($m, $plotfile, $file, $ist1) = @_;

    my $M = 0;

    open(OUT, ">$plotfile") || die;

    my $time = -1;
    my $tmm = 0;
    my $tmd = 0;
    my $tmi = 0;
    my $tii = 0;
    my $tdd = 0;

    my $tmm_slope = 0;
    my $tmd_slope = 0;
    my $tmi_slope = 0;
    my $tii_slope = 0;
    my $tdd_slope = 0;

    my $tmm_linear = 0;
    my $tmd_linear = 0;
    my $tmi_linear = 0;
    my $tii_linear = 0;
    my $tdd_linear = 0;
    
    open (FILE, "$file") || die;
    while (<FILE>) {
	if (/^M=(\d+)\s+muD/) { $M = $1; }
	elsif (/^TIME\s+=\s+(\S+)/) {
	    $time = $1;
	}
	#the rates
	elsif (/^$m\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)$/) {
	    my $mudS     = $1;
	    my $mudD     = $2;
	    my $muES     = $3;
	    my $muEI     = $4;
	    my $ldS      = $5;
	    my $ldI      = $6;
	    my $xiz      = $7;
	    my $rz       = $8;
	    my $xiIinfty = $9;

	    # the slopes
	    $tmm_slope = 0;
	    $tmd_slope = $mudS;
	    $tmi_slope = 0;
	    $tii_slope = 0;
	    $tdd_slope = $mudD;
	    
	    $tii_slope = ($ldI - $muEI*$xiz)*(1.0-$xiz);

	    if ( abs( $muES*(1.0-$rz)-($muEI-$rz*$ldI) )  > 0.001) {
		printf "Rates do not satisfy the necessary conditions %f", abs( $muES*(1.0-$rz)-($muEI-$rz*$ldI) ); die;
	    }
	    else {
		$tmi_slope = $ldS;
	    }

	    $tmm_slope = - $tmi_slope - $tmd_slope;
	    if ($tdd_slope < 0.) { print "bad slope for tDD $tdd_slope\n"; die; }
	    if ($tii_slope < 0.) { print "bad slope for tII $tii_slope\n"; die; }

	    if ($tmi_slope < 0.) { print "bad slope for tMI $tmi_slope\n"; die; }
	    if ($tmd_slope < 0.) { print "bad slope for tMD $tmd_slope\n"; die; }
	    if ($tmm_slope > 0.) { print "bad slope for tMM $tmm_slope\n"; die; }
	}
	#the transition probabilities
	elsif (/^$m\s+$time\s+\|\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)/) {
	    $tmm = $1;
	    $tmi = $2;
	    $tmd = $3;
	    $tii = $4;
	    $tdd = $5;

	    $tmm_linear = $time *$tmm_slope + 1.0;
	    $tmi_linear = $time *$tmi_slope;
	    $tmd_linear = $time *$tmd_slope;
	    $tii_linear = $time *$tii_slope + $xiz;
	    $tdd_linear = $time *$tdd_slope;

	    if (!$ist1) {
		print OUT "$time\t$tmm\t$tmm_linear\t$tmi\t$tmi_linear\t$tmd\t$tmd_linear\t$tii\t$tii_linear\t$tdd\t$tdd_linear\n";
	    }
	    else {
		if ($time == 1.0) {
		    print OUT "$time\t$tmm\t$tmi\t$tmd\t$tii\t$tdd\n";
		}
	    }
	}

    }
    close (FILE);
    close (OUT);
 
    return $M;
}

sub plot_evohmm_utest {
    my ($wrkdir,  $which) = @_;

    my $name = $wrkdir;
    if ($name =~ /^\S+\/([^\/]+)$/) { $name = $1; }

    my $xiz_def      = 0.01; 
    my $rz_def       = 0.01; 
    my $xiIinfty_def = 0.99;

    my $xiz; 
    my $rz; 
    my $xiIinfty;

    my $file1;
    my $file2;
    my $file3;
    my $file4;
    my $file5;

    if ($which == 1) {
	$rz       = $rz_def; 
	$xiIinfty = $xiIinfty_def;
	
	$file1 = "$wrkdir/$name\_xiz$xiz1\_rz$rz\_xiIinfty$xiIinfty.out";
	$file2 = "$wrkdir/$name\_xiz$xiz2\_rz$rz\_xiIinfty$xiIinfty.out";
	$file3 = "$wrkdir/$name\_xiz$xiz3\_rz$rz\_xiIinfty$xiIinfty.out";
	$file4 = "$wrkdir/$name\_xiz$xiz4\_rz$rz\_xiIinfty$xiIinfty.out";
	$file5 = "$wrkdir/$name\_xiz$xiz5\_rz$rz\_xiIinfty$xiIinfty.out";
	
	create_plots($wrkdir, $name, "xiz", $file1, $file2, $file3, $file4, $file5);
    }
    elsif ($which == 2) {
	$xiz      = $xiz_def; 
	$xiIinfty = $xiIinfty_def;
 	
	$file1 = "$wrkdir/$name\_xiz$xiz\_rz$rz1\_xiIinfty$xiIinfty.out";
	$file2 = "$wrkdir/$name\_xiz$xiz\_rz$rz2\_xiIinfty$xiIinfty.out";
	$file3 = "$wrkdir/$name\_xiz$xiz\_rz$rz3\_xiIinfty$xiIinfty.out";
	$file4 = "$wrkdir/$name\_xiz$xiz\_rz$rz4\_xiIinfty$xiIinfty.out";
	$file5 = "$wrkdir/$name\_xiz$xiz\_rz$rz5\_xiIinfty$xiIinfty.out";
	
	create_plots($wrkdir, $name, "rz", $file1, $file2, $file3, $file4, $file5);
    }
     elsif ($which == 3) {
	$xiz      = $xiz_def; 
	$rz       = $rz_def; 
	
	$file1 = "$wrkdir/$name\_xiz$xiz\_rz$rz\_xiIinfty$xiIinfty1.out";
	$file2 = "$wrkdir/$name\_xiz$xiz\_rz$rz\_xiIinfty$xiIinfty2.out";
	$file3 = "$wrkdir/$name\_xiz$xiz\_rz$rz\_xiIinfty$xiIinfty3.out";
	$file4 = "$wrkdir/$name\_xiz$xiz\_rz$rz\_xiIinfty$xiIinfty4.out";
	$file5 = "$wrkdir/$name\_xiz$xiz\_rz$rz\_xiIinfty$xiIinfty5.out";
	
	create_plots($wrkdir, $name, "xiIinfty", $file1, $file2, $file3, $file4, $file5);
    }

}

sub plot_evohmm_utest_one {
    my ($wrkdir, $xiz, $rz, $xiIinfty) = @_;

    my $name = $wrkdir;
    if ($name =~ /^\S+\/([^\/]+)$/) { $name = $1; }

 
    my $file = "$wrkdir/$name.out";

    create_plot_one($wrkdir, $name, $file);
}



sub run_evohmm_utest {
    my ($evohmm_utest, $wrkdir, $hmmfile, $hmmname, $xiz, $rz, $xiIinfty) = @_;
    
    my $evohmmfile = "$wrkdir/evo$hmmname\_xiz$xiz\_rz$rz\_xiIinfty$xiIinfty.hmm";
    my $outfile    = "$wrkdir/evo$hmmname\_xiz$xiz\_rz$rz\_xiIinfty$xiIinfty.out";
    
    my $cmd  = "$evohmm_utest --xiz $xiz  --rz $rz --xiIinfty $xiIinfty ";
    system("echo $cmd $hmmfile $evohmmfile\n");
    system("$cmd $hmmfile $evohmmfile > $outfile\n");
}
