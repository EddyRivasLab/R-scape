#!/usr/bin/perl -w
#YlueSimpson.pl

use strict;
use Class::Struct;
use lib '/Users/rivase/src/src/mysource/scripts';
use PDBFUNCS;
use FUNCS;
use constant GNUPLOT => '/usr/local/bin/gnuplot';

use vars qw ($opt_P $opt_D $opt_L $opt_v);  # required if strict used
use Getopt::Std;
getopts ('P:D:L:v');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  YuleSimpson.pl [options] <stofile> \n\n";
        print "options:\n";
 	exit;
}

my $stofile = shift;
my $stoname = $stofile;
if ($stoname =~ /([^\/]+)s*$/) { $stoname = $1; }
if ($stoname =~ /^(\S+).sto$/) { $stoname = $1; }

my $rscape = "/Users/rivase/src/src/mysource/bin/R-scape ";

my $seeplots = 1;

my $maxD = 8;
if ($opt_D) { $maxD = $opt_D; }

my $minL = 1;
if ($opt_L) { $minL = $opt_L; }

my $pdbfile = "";
if ($opt_P) { $pdbfile = "$opt_P"; }

my $sc_cutoff = -500000000000;

my $YS = 0;
my $outfile = "$stoname.rscape";
run_rscape($stofile, $pdbfile, $maxD, $minL, $outfile, $YS);

$YS = 1;
my $outfileYS = "$stoname.YS.rscape";
run_rscape($stofile, $pdbfile, $maxD, $minL, $outfileYS, $YS);

plot_YS($outfile, $outfileYS);

sub run_rscape {
    my ($stofile, $pdbfile, $maxD, $minL, $outfile, $YS) = @_;

    my $cmd = "$rscape --GTp  --seed 0  --naive  ";
    if ($pdbfile) { $cmd .= " --pdbfile $pdbfile --cntmaxD $maxD --cntmind $minL "; }
    if ($YS)      { $cmd .= " --YS "; }
    
    system("echo $cmd $stofile \n");
    system("     $cmd $stofile > $outfile\n");
}

sub plot_YS {
    my ($outfile, $outfileYS) = @_;

    my $plotfile1   = "$outfile.1.ys";
    my $plotfile2   = "$outfile.2.ys";
    my $plotfileYS1 = "$outfileYS.1.ys";
    my $plotfileYS2 = "$outfileYS.2.ys";
    
    my $nbps = 0;
    my $noth = 0;
    my @bps_sc;
    my @oth_sc;
    my $break = -1;
    my $minsc;
    my $maxsc;
    parse_rscapeout($outfile, \$nbps, \@bps_sc, \$noth, \@oth_sc);
    my $n = create_plotfile($plotfile1, $plotfile2, $nbps, \@bps_sc,  $noth, \@oth_sc, \$break, \$minsc, \$maxsc);
    print "break $break\n";
    
    my $nbps_YS = 0;
    my $noth_YS = 0;
    my @bps_YS_sc;
    my @oth_YS_sc;
    my $breakYS = -1;
    parse_rscapeout($outfileYS, \$nbps_YS, \@bps_YS_sc, \$noth_YS, \@oth_YS_sc);
    create_plotfile($plotfileYS1, $plotfileYS2, $nbps_YS, \@bps_YS_sc,  $noth_YS, \@oth_YS_sc, \$breakYS);
    print "nbreak $breakYS\n";
   
    YSplot($plotfile1, $plotfile2, $plotfileYS1, $plotfileYS2, $break, $minsc, $maxsc, $n+10);
}

sub parse_rscapeout {
    my ($outfile, $ret_nbps, $bps_sc_ref, $ret_noth, $oth_sc_ref) = @_;

    my $nbps = 0;
    my $noth = 0;

    open(FILE, "$outfile") || die;
    while(<FILE>) {
	if (/#/) {
	}
	elsif (/^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+\S+\s*$/) {
	    my $type = $1;
	    my $i    = $2;
	    my $j    = $3;
	    my $sc   = $4;
	    
	    if    ($type =~ /^\*$/)   {  # A WWc basepair
		$bps_sc_ref->[$nbps] = $sc;
		$nbps ++;
	    }
	    elsif ($type =~ /^\*\*$/) {  # A basepair
		$bps_sc_ref->[$nbps] = $sc;
		$nbps ++;
	    }
	    elsif ($type =~ /^c/)     {  # A contact
		$bps_sc_ref->[$nbps] = $sc;
		$nbps ++;
	    }
	    elsif ($type =~ /^\~$/)    { # Not a contact although compatible with the bpairs
		$oth_sc_ref->[$noth] = $sc;
		$noth ++;
	    }
	    else { print "type? $type\n"; die; }
	}
 	elsif (/\s+(\d+)\s+(\d+)\s+(\S+)\s+\S+\s*$/) { # Not annotated as a contact
	    my $i        = $1;
	    my $j        = $2;
	    my $sc       = $3;
	    $oth_sc_ref->[$noth] = $sc;
	    $noth ++;
	}		

    }
    close(FILE);
    $$ret_nbps = $nbps;
    $$ret_noth = $noth;
    print "$outfile: nbps $nbps noth $noth\n";
}

sub create_plotfile {
    my ($plotfile1, $plotfile2, $nbps, $bps_sc_ref,  $noth, $oth_sc_ref, $ret_break, $ret_minsc, $ret_maxsc) = @_;

    open(PLOT1, ">$plotfile1") || die;
    my $maxsc = -123456789;
    my $minsc =  123456789;
    my $n = 1;
    for (my $o = $noth-1; $o >= 0; $o--) {
	my $sc = $oth_sc_ref->[$o];
	if ($sc > $sc_cutoff) {
	    print PLOT1 "$n $sc\n";
	    if ($sc < $minsc) { $minsc = $sc; }
	    if ($sc > $maxsc) { $maxsc = $sc; }
	    $n ++;
	}
    }
    close(PLOT1);
    $$ret_break = $n;
    
    open(PLOT2, ">$plotfile2") || die;
    for (my $b = 0; $b < $nbps; $b++) {
	my $sc = $bps_sc_ref->[$b];
	if ($sc > $sc_cutoff) {
	    print PLOT2 "$n $sc\n";
	    if ($sc < $minsc) { $minsc = $sc; }
	    if ($sc > $maxsc) { $maxsc = $sc; }
	    $n ++;
	}
    }
    close(PLOT2);

    $$ret_minsc = $minsc;
    $$ret_maxsc = $maxsc;

    return $n;
}

sub YSplot {
    my ($plotfile1, $plotfile2, $plotfileYS1, $plotfileYS2, $break, $minsc, $maxsc, $max) = @_;

    my $psfile = "$plotfile1.ps";
    
    print "FILE: $psfile\n";

    my $xlabel;
    my $ylabel;
    my $title  = "$stoname";
    my $x;
    my $y;
    my $x_max, my $x_min;
    my $y_max, my $y_min;
    
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    my $gp = \*GP;

    print $gp "set terminal postscript color solid 14\n";
    print $gp "set output '$psfile'\n";    
    print $gp "set style line 1   lt 1 lc rgb 'black' pt 1 ps 0.5 lw 1\n";
    print $gp "set style line 2   lt 1 lc rgb 'brown' pt 1 ps 0.5 lw 1\n";
    print $gp "set style line 3   lt 1 lc rgb 'grey' pt 1 ps 0.5 lw 0.5\n";
    print $gp "set style line 4   lt 1 lc rgb 'cyan' pt 1 ps 0.5 lw 3\n";
    print $gp "set style line 7   lt 1 lc rgb 'red' pt 1 ps 0.5 lw 3\n";
    print $gp "set style line 5   lt 1 lc rgb 'purple' pt 1 ps 0.5 lw 3\n";
    print $gp "set style line 6   lt 1 lc rgb 'orange' pt 1 ps 0.5 lw 3\n";
    print $gp "set style line 8   lt 1 lc rgb 'blue' pt 1 ps 0.5 lw 3\n";
    print $gp "set style line 9   lt 2 lc rgb 'magenta' pt 1 ps 0.5 lw 3\n";

    print $gp "set title  '$stoname'\n";
    print $gp "set xlabel 'ordered pairs'\n";
    print $gp "set ylabel 'covariation score'\n";
   # print $gp "set xrange [1:$max]\n";
   # print $gp "set yrange [$minsc:$maxsc]\n";

    my $key   = "given alignment";
    my $keyYS = "row-shuffled alignment";
    my $m = 4;
    my $cmd = "";
    $cmd  = "'$plotfile1'   using 1:2  title '$key'   ls $m, "; $m ++; if ($m == 7) { $m = 8; }
    $cmd .= "'$plotfile2'   using 1:2  title '$key'   ls $m";   $m ++; if ($m == 7) { $m = 8; }
    print $gp "plot $cmd\n";
    
    $cmd  = "'$plotfileYS1' using 1:2  title '$keyYS' ls $m, ";   $m ++; if ($m == 7) { $m = 8; }
    $cmd .= "'$plotfileYS2' using 1:2  title '$keyYS' ls $m";     $m ++; if ($m == 7) { $m = 8; }
    print $gp "plot $cmd\n";
    
    close($gp);
    if ($seeplots) { system ("open $psfile&\n"); }

}
