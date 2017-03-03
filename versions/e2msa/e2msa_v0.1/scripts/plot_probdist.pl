#!/usr/bin/perl -w
#plot_probdist.pl
 
use strict;
use Class::Struct;
use lib '/groups/eddy/home/rivase/projects/evohmm/scripts';
use lib '/Users/rivase/projects/evohmm/scripts';
use FUNCS;
#use constant GNUPLOT => '/usr/bin/gnuplot';
use constant GNUPLOT => '/opt/local/bin/gnuplot';
#use constant GNUPLOT => '/sw/bin/gnuplot';

use vars qw ($opt_i $opt_I $opt_f $opt_F $opt_t);  # required if strict used
use Getopt::Std;
getopts ('i:I:f:F:t:');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  plot_prodist.pl [options] <file>  \n\n";
        print "options:\n";
	exit;
}

my $file = shift;
my $fmin = -1.0;
my $fmax = -1.0;
if ($opt_f) { $fmin = $opt_f; if ($fmin < 0.0)  { print "ilegal fmin $fmin\n"; die; } }
if ($opt_F) { $fmax = $opt_F; if ($fmax < 0.0)  { print "ilegal fmax $fmax\n"; die; } }
if ($fmin >= 0 && $fmax >= 0 && $fmin >= $fmax) { print "ilegal fmin > fmax \n"; die; }

#my $tmax = 10.;
#if ($opt_t) { $tmax = $opt_t; if ($tmax < 0.0) { print "ilegal tmax $tmax\n";  die; } }

my $minid = 0.;
my $maxid = 100.;
if ($opt_i) { $minid = $opt_i; if ($minid < 0.0) { print "ilegal minid $minid\n"; die; } }
if ($opt_I) { $maxid = $opt_I; if ($maxid < 0.0) { print "ilegal maxid $maxid\n"; die; } }

plot_probdist($file, $minid, $maxid, $fmin, $fmax);



sub plot_probdist {

    my ($file, $minid, $maxid, $fmin, $fmax, $subs) = @_;

    my $nf = 0;
    my $nft = 0;
    my $doplot = 0;
    my $psfile = "$file.ID_$minid-$maxid.ps";

    my $title  = "$file.ID[$minid,$maxid]"; if ($title =~ /\/(\S+)$/) { $title = $1; }
    my $ylabel = "LogProbability P(s1,s2 | t) / len";
    my $xlabel = "transitions divergence distance t";
 
    my $minsc = 0;
    my $maxsc = -1000000.;

    my $id;
    my $match;
    my $len;
    my $optime;
    my $optsc;

    my $tmin = 1000000.;
    my $tmax = -1.0;
 
    my $max_optime = -1.0;
    my $min_optime = 1e+10;
    my $avg_optime = 0.0;
    my $std_optime = 0.0;

    open(FILE, "$file")|| die;
    while (<FILE>) {
	if (/^\#\s+id\s+(\S+)\s+match\s+(\S+)\s+L\s+(\S+)\s+optime\s+(\S+)\s+optsc\s+(\S+)/) {
	    $id     = $1;
	    $match  = $2;
	    $len    = $3;
	    $optime = $4;
	    $optsc  = $5;
	    $doplot = 0;

	    if ($id >= $minid && $id <= $maxid) { 
		$doplot = 1; 
		$nf ++; 
		if ($optime > $max_optime) { $max_optime = $optime; }
		if ($optime < $min_optime) { $min_optime = $optime; }
		$avg_optime += $optime;
		$std_optime += $optime*$optime;
	    }
	    $nft ++;
	}
	elsif (/^(\d\S+)\s+(\S+)\s*$/) {
	    my $t  = $1;
	    my $sc = $2;
	    my $scl = $sc/$len;
	    if ($doplot) {
		if ($scl < $minsc) { $minsc = $scl; }
		if ($scl > $maxsc) { $maxsc = $scl; }
		if ($t   < $tmin && $t > 1e-8)  { $tmin = $t; }
		if ($t   > $tmax)               { $tmax = $t; }
	    }
	}
    }
    close(FILE);
    FUNCS::calculate_averages (\$avg_optime, \$std_optime, $nf);

    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    print GP "set output '$psfile'\n";
    print GP "set key right top\n";
    print GP "set title  '$title'\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    print GP "set nokey\n";
   
    printf GP "set xrange [%f:%f]\n", $tmin, $tmax+1000;
    printf GP "set yrange [%f:%f]\n", $minsc-0.1*abs($maxsc-$minsc), $maxsc+0.1*abs($maxsc-$minsc);
    print GP "set logscale x\n";

    print GP "set multiplot\n";
    my $style = 1;
    $doplot = 0;
    $nf = 0;
    open(FILE, "$file")|| die;
    while (<FILE>) {
	if (/^&$/) {
	    $style ++;
	    if ($style > 9) { $style = 0; }
	    if ($doplot) { 
		print GP "e\n";  
		printf GP "plot '-' using 1:2 ls 6666\n";
		printf GP "%f %f\n", $optime, $optsc;
		
		printf GP "e\n"; 
		#printf("%d | %f %f \n", $nf, $optime, $optsc);
	    }
	    if ($fmax < 0 || ($fmax >= $fmin && $nf >= $fmin && $nf <= $fmax)) { }
	    else {     
		close(FILE);
		close(GP);
		printf("NF $nf\n");
		system("open $psfile\n");
	    }
	}
	elsif (/^\#\s+id\s+(\S+)\s+match\s+(\S+)\s+L\s+(\S+)\s+optime\s+(\S+)\s+optsc\s+(\S+)\s*$/) {
	    $id       = int($1*10)/10.;
	    $match    = int($2*10)/10.;
	    $len      = $3;
	    $optime   = $4;
	    $optsc    = $5/$len;

	    $doplot = 0;
	    if ($id >= $minid && $id <= $maxid) { $doplot = 1; }
	    if ($doplot) { printf GP "plot '-' using 1:2 title 'id=$id' with lines ls %d\n", $style; $nf ++; }
 	}
	elsif (/^(\d\S+)\s+(\S+)\s*$/) {
	    my $t = $1;
	    my $scl = $2/$len;

	    if ($doplot) {  
		printf GP "%f %f\n", $t, $scl; 
	    }
	}
	
    }
    close(FILE);
    close(GP);

    printf("FILE = $file\n");
    printf("ID [%.4f,%.4f] NF %d/%d\n", $minid, $maxid, $nf, $nft);
    printf("optime %.4f +\= %.4f  [%.4f,%.4f]\n", $avg_optime, $std_optime, $min_optime, $max_optime);
 
    system("open $psfile&\n");

}


