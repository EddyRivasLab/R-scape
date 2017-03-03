#!/usr/bin/perl -w
#plot_evohmm.pl

use strict;
use Class::Struct;
use lib '/groups/eddy/home/rivase/projects/evohmm/scripts';
use lib '/Users/rivase/projects/evohmm/scripts';
use FUNCS;
use constant GNUPLOT => '/usr/bin/gnuplot';

use vars qw ();  # required if strict used
use Getopt::Std;
getopts ('');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  plot_evohmm.pl [options] <.plot>\n\n";
        print "options:\n";
	exit;
}

my $plotfile  = shift;
my $outfile  = "$plotfile.out";

my $time_45 = -1;
my $time_50 = -1;
my $time_62 = -1;
my $time_80 = -1;
my $time_90 = -1;
my $alpha_45 = -1;
my $alpha_50 = -1;
my $alpha_62 = -1;
my $alpha_80 = -1;
my $alpha_90 = -1;
my $beta_45 = -1;
my $beta_50 = -1;
my $beta_62 = -1;
my $beta_80 = -1;
my $beta_90 = -1;

parse_plotfile($plotfile, $outfile);

my $plot_45 = "$plotfile.45";
my $plot_50 = "$plotfile.50";
my $plot_62 = "$plotfile.62";
my $plot_80 = "$plotfile.80";
my $plot_90 = "$plotfile.90";

open(FILE, ">$plot_45") || die;
printf       "45: %f %f %f %f \n", $time_45, $alpha_45, $beta_45, $alpha_45/$beta_45;
printf FILE "%f %f \n", $time_45, 0;
printf FILE "%f %f \n", $time_45, -10000;
printf FILE "%f %f \n", $time_45, $alpha_45;
printf FILE "%f %f \n", $time_45, $beta_45;
close(FILE);

open(FILE, ">$plot_50") || die;
printf       "50: %f %f %f %f \n", $time_50, $alpha_50, $beta_50, $alpha_50/$beta_50;
printf FILE "%f %f \n", $time_50, 0;
printf FILE "%f %f \n", $time_50, -100000;
printf FILE "%f %f \n", $time_50, $alpha_50;
printf FILE "%f %f \n", $time_50, $beta_50;
close(FILE);

open(FILE, ">$plot_62") || die;
printf       "62: %f %f %f %f \n", $time_62, $alpha_62, $beta_62, $alpha_62/$beta_62;
printf FILE "%f %f \n", $time_62, 0;
printf FILE "%f %f \n", $time_62, -100000;
printf FILE "%f %f \n", $time_62, $alpha_62;
printf FILE "%f %f \n", $time_62, $beta_62;
close(FILE);

open(FILE, ">$plot_80") || die;
printf       "80: %f %f %f %f \n", $time_80, $alpha_80, $beta_80, $alpha_80/$beta_80;
printf FILE "%f %f \n", $time_80, 0;
printf FILE "%f %f \n", $time_80, -100000;
printf FILE "%f %f \n", $time_80, $alpha_80;
printf FILE "%f %f \n", $time_80, $beta_80;
close(FILE);

open(FILE, ">$plot_90") || die;
printf       "90: %f %f %f %f \n", $time_90, $alpha_90, $beta_90, $alpha_90/$beta_90;
printf FILE "%f %f \n", $time_90, 0;
printf FILE "%f %f \n", $time_90, -100000;
printf FILE "%f %f \n", $time_90, $alpha_90;
printf FILE "%f %f \n", $time_90, $beta_90;
close(FILE);


gnuplot($outfile);

sub gnuplot {
    my ($file) = @_;

    my $name = "";
    my $xlabel = "";
    my $ylabel = "";

    my $psfile = "$file.ps";
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    print GP "set output '$psfile'\n";
    print GP "set key right top\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    print GP "set xrange  [0.1:5]\n";
    print GP "set yrange [-18:0]\n";
    print GP "set nokey\n";
    my $cmd = "";
    $cmd  = "'$file' using 1:2 with lines ls 1, '$file' using 1:3 with lines ls 2, ";
    $cmd  .= "'$plot_45' using 1:2 with lines ls 1, ";
    $cmd  .= "'$plot_62' using 1:2 with lines ls 1, ";
    $cmd  .= "'$plot_90' using 1:2 with lines ls 1";
    print GP "plot $cmd\n";
   print "plot $cmd\n";
  
    close(GP);

    system("evince $psfile\n");
}

sub parse_plotfile {
    my ($file, $outfile) = @_;

    my $lambda = -1;
    my $lambdainv;
    my $lambdastar;

    open(OUT, ">$outfile") || die;

     # first pass to get lambda
    open (FILE, "$file") || die;
    while (<FILE>) {
	if (/^(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+/) {
	    my $time = $1;
	    my $tmm  = $2;
	    my $tmi  = $3;
	    my $tmd  = $4;
	    my $tii  = $5;
	    my $tdd  = $6;

	    #print "$time\t$tmm\t$tmi\t$tii\n";
	    my $beta = log($tii);
	    my $alpha;
	    if (1.0-$tii == 0.0 || $tmi == 0.0 || $tmm == 0.0) { $alpha = -100000.0; }
	    else                                               { $alpha = log(1.0-$tii) + log($tmi) - log($tmm); }


	    if ($time == 1.0) { $lambdastar = -log($tii);  }
	    if ($time > 0. && $alpha > -100000.0 && $alpha/$beta < 11.0) { $lambda = -log($tii);  $time_62 = $time; }
	}
    }
    close(FILE);

    if ($lambda == -1) { print "not reached\n"; $lambda = $lambdastar; $time_62 = 1.0; }
   if ($lambda != 0.0) { $lambdainv = 1.0/$lambda; }
    printf("time %f lambda %f lambdainv %f\n", $time_62, $lambda, $lambdainv);

    # second pass calculate alpha and beta
    open (FILE, "$file") || die;
    while (<FILE>) {
	
	if (/^(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+/) {
	    my $time = $1;
	    my $tmm  = $2;
	    my $tmi  = $3;
	    my $tmd  = $4;
	    my $tii  = $5;
	    my $tdd  = $6;

	    my $beta = log($tii);
	    my $alpha;
	    if (1.0-$tii == 0.0 || $tmi == 0.0 || $tmm == 0.0) { $alpha = -100000.0; }
	    else                                               { $alpha = log(1.0-$tii) + log($tmi) - log($tmm); }

	    $alpha *= $lambdainv;
	    $beta  *= $lambdainv;

	    printf     "%f\t%f\t%f\t%f\n", $time, $alpha, $beta, $alpha/$beta;
	    printf OUT "%f\t%f\t%f\t%f\n", $time, $alpha, $beta, $alpha/$beta;

	    if ($time < $time_62)                 {                   $alpha_62 = $alpha; $beta_62 = $beta; }
	    if ($time < $time_62 * 31.088/28.300) { $time_50 = $time; $alpha_50 = $alpha; $beta_50 = $beta; }
	    if ($time < $time_62 * 34.084/28.300) { $time_45 = $time; $alpha_45 = $alpha; $beta_45 = $beta; }
	    if ($time < $time_62 * 21.105/28.300) { $time_80 = $time; $alpha_80 = $alpha; $beta_80 = $beta; }
	    if ($time < $time_62 * 18.024/28.300) { $time_90 = $time; $alpha_90 = $alpha; $beta_90 = $beta; }
	}
    }


    close(FILE);
    close(OUT);
}
