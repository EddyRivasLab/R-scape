#!/usr/bin/perl -w
#plot_evosim.pl
 
use strict;
use Class::Struct;
use lib '/groups/eddy/home/rivase/projects/evohmm/scripts';
use lib '/Users/rivase/projects/evohmm/scripts';
use FUNCS;
use constant GNUPLOT => '/usr/bin/gnuplot';
#use constant GNUPLOT => '/opt/local/bin/gnuplot';
#use constant GNUPLOT => '/sw/bin/gnuplot';

use vars qw ($opt_D $opt_e $opt_E $opt_F $opt_i $opt_j $opt_I $opt_K $opt_l $opt_L $opt_n $opt_N $opt_O $opt_P $opt_t $opt_T $opt_u $opt_U $opt_v $opt_r $opt_s $opt_x $opt_y $opt_Y $opt_w $opt_W $opt_z);  # required if strict used
use Getopt::Std;
getopts ('D:e:E:Fi:jI:Kl:L:u:U:gnN:OPt:T:rs:vx:y:Y:w:W:z:');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  plot_evosim.pl [options] <srcdir> <wrkdir> \n\n";
        print "options:\n";
	print "-L    :  length of ancestral sequence [default is 1000]\n";
	print "-F    :  goodness of a geometric fit to insert lengths\n";
	print "-N    :  number of sequences generated per time [default is 1000]\n";
	print "-j    :  sample ancestral sequences instead of fix length\n";
	print "-r    :  impose reversibility\n";
	print "-x    :  rI value      \n";
	print "-l    :  ldE \n";
	print "-K    :  run TKFsim \n";
	print "-I    :  ldI \n";
	print "-u    :  muE \n";
	print "-U    :  muI \n";
	print "-D    :  muA \n";
	print "-y    :  sI \n";
	print "-w    :  sD \n";
	print "-Y    :  vI \n";
	print "-W    :  vD \n";
	print "-n    :  don't view plots\n";
 	print "-v    :  be verbose\n";
        print "-t    :  tmin [default is 0]\n";
        print "-T    :  tend [default is 3]\n";
        print "-z    :  tinterval [default is 0.2]\n";
        print "-i    :  tinc [default is 0.01]\n";
        print "-e    :  tepsilon [default is 1e-5]\n";
        print "-s    :  seed <n>\n";
        print "-O    :  plot in logarithmic scale\n";
        print "-P    :  don't run e1sim, just redo the plots\n";
	exit;
}
my $srcdir   = shift;
my $wrkdir   = shift;

my $tP10 = 0.114220;
my $tP30 = 0.332924;
my $tB90 = 0.903079;
my $tB62 = 1.427753;

my $filename = $wrkdir;
if ($filename =~ /^\S+\/([^\/]+)$/) { $filename = $1; }

if ( -d $wrkdir) { 
    print "wrkdir directory $wrkdir exist.\n";
    #system("rm $wrkdir/*\n"); 
}
else { system("mkdir $wrkdir\n"); } 
my $view = 1;
if ($opt_n) { $view = 0; }
my $verbose = 0;
if ($opt_v) { $verbose = 1; }

my $logarithmic = 0;
if ($opt_O) { $logarithmic = 1; }

# variables of the simulation
my $L         = 1000;
my $N         = 1000;
my $tmin      = 0.0;
my $tmax      = 3.0;
my $tinc      = 0.05;
my $teps      = 1e-5;
my $tinterval;

if ($opt_L) { $L         = $opt_L; if ($L         < 0) { print "L    has to be positive\n"; die; } }
if ($opt_N) { $N         = $opt_N; if ($N         < 0) { print "N    has to be positive\n"; die; } }
if ($opt_t) { $tmin      = $opt_t; if ($tmin      < 0) { print "tmin has to be positive\n"; die; } }
if ($opt_T) { $tmax      = $opt_T; if ($tmax      < 0) { print "tmax has to be positive\n"; die; } }
if ($opt_z) { $tinterval = $opt_z; if ($tinterval < 0) { print "tinterval has to be positive\n"; die; } } else { $tinterval = $tmax-$tmin; }
if ($opt_i) { $tinc      = $opt_i; if ($tinc      < 0) { print "tinc has to be positive\n"; die; } }
if ($opt_e) { $teps      = $opt_e; if ($teps      < 0) { print "teps has to be positive\n"; die; } }

# rate variables 
my $evomodel = "UN";
my $rI   = 0.01; 
my $ldE  = 0.0;
my $muE  = 0.0;
my $ldI  = 0.0;
my $muI  = 0.0;
my $muA  = 0.0;
my $sI   = 0.0;
my $sD   = 0.0;
my $vI   = 0.0;
my $vD   = 0.0;

if ($opt_E) { $evomodel   = $opt_E; }
if ($opt_x) { $rI         = $opt_x; if ( $rI  < 0 || $rI > 1.0) { print "bad value of rI=$rI\n";           die; } }
if ($opt_l) { $ldE        = $opt_l; if ( $ldE < 0)              { print "bad value of ldE=$ldE\n";         die; } }
if ($opt_I) { $ldI        = $opt_I; if ( $ldI < 0)              { print "bad value of ldI=$ldI\n";         die; } }
if ($opt_u) { $muE        = $opt_u; if ( $muE < 0)              { print "bad value of muE=$muE\n";         die; } }
if ($opt_U) { $muI        = $opt_U; if ( $muI < 0)              { print "bad value of muE=$muE\n";         die; } }
if ($opt_D) { $muA        = $opt_D; if ( $muA < 0)              { print "bad value of muA=$muA\n";         die; } }
if ($opt_y) { $sI         = $opt_y; }
if ($opt_w) { $sD         = $opt_w; }
if ($opt_Y) { $vI         = $opt_Y; }
if ($opt_W) { $vD         = $opt_W; }

my $joint = 0;
if ($opt_j) { $joint = 1; }
my $reversible = 0;
if ($opt_r) { $reversible = 1; }
my $tkfsim = 0;
if ($opt_K) { $tkfsim = 1; }

my $e1simp   = "$srcdir/src/e1sim";
my $tkfsimp  = "$srcdir/src/tkfsim";
my $program = $e1simp;
if ($opt_K) { $program = $tkfsimp; }

my $geofit = 0;
if ($opt_F) { $geofit = 1; }

if (!$opt_P) {
    if ($tkfsim) {
	run_tkfsim($program, $wrkdir, $filename, $L, $N, $tmin, $tmax, $tinterval, $tinc, $teps, $rI, $ldI, $muI, $joint, $reversible, $geofit, $verbose);
    }
    else {
	run_e1sim($program, $wrkdir, $filename, $L, $N, $tmin, $tmax, $tinterval, $tinc, $teps, $evomodel, $rI, $ldE, $muE, $ldI, $muI, $muA, $sI, $sD, $vI, $vD, $joint, $reversible, $geofit, $verbose);
    }
}

if ($opt_P) {use constant GNUPLOT => '/opt/local/bin/gnuplot'; }
plot_program($wrkdir, $reversible, $tkfsim, $view);
if (!$logarithmic) { plot_geofit($wrkdir, $reversible, $tkfsim, $view); }

 
sub actual_plot{
    my ($name, $tstar, $plotf, $plotf_tstar, $L, $N, $joint, $reversible, $tkfsim,
	$tmin, $tmax, $tinterval, $tinc, $teps, $rI, $ldE, $muE, $ldI, $muI, $muA, 
	$sE_B62, $nE_B62, $eE_B62, $sE_B90, $nE_B90, $eE_B90, 
	$sE_P30, $nE_P30, $eE_P30, $sE_P10, $nE_P10, $eE_P10, 
	$view) = @_;

    my $xlabel = "evolutionary distance t";
    my $ylabel = "";

    my $plot = $plotf; if ($plot =~ /^(\S+).plot$/) { $plot = $1; }
    my $statsfile    = "$plot\_stats.ps";
    my $Pfile        = "$plot\_prob.ps";
    my $Lfile        = "$plot\_logp.ps";
    my $pdfstatsfile = "$plot\_stats.pdf";
    my $pdfPfile     = "$plot\_prob.pdf";
    my $pdfLfile     = "$plot\_logp.pdf";

    my $subsmin = $tmin*$tB62;
    my $subsmax = $tmax*$tB62;

    my $tB62file   = "$plot.tB62";
    my $tB90file   = "$plot.tB90"; 
    my $tPAM30file = "$plot.tP30"; 
    my $tPAM10file = "$plot.tP10";
 
    my $k;
    my $eta_slope;
    my $eta_inf;
    my $beta_slope;
    my $beta_inf;

    $eta_slope  = (1-$rI) * ($ldI-$rI*$muI);
    $beta_slope = $ldE;
    $eta_slope  = int($eta_slope  * 10000)/10000;
    $beta_slope = int($beta_slope * 100000000)/100000000;
 
    $eta_inf  = ($ldI<$muI)? $ldI/$muI : 1.0; 
    $beta_inf = ($ldE+$muE > 0)? $ldE/($ldE+$muE) : 0; 
    $eta_inf  = int($eta_inf  * 1000000)/1000000;
    $beta_inf = int($beta_inf * 10000000)/10000000;
     
    my $maxL = ($eta_inf < 1.0)? $L + ($L+1.)*$beta_inf/(1-$eta_inf) : 2*$L;
    $maxL = $L;

    open(T1, ">$tB62file") || die;
    print T1 "$tB62\t 0\n";
    print T1 "$tB62\t $maxL\n";
    close(T1);

    open(T1, ">$tB90file") || die;
    print T1 "$tB90\t 0\n";
    print T1 "$tB90\t $maxL\n";
    close(T1);
 
    open(T1, ">$tPAM30file") || die;
    print T1 "$tP30\t 0\n";
    print T1 "$tP30\t $maxL\n";
    close(T1);   
 
    open(T1, ">$tPAM10file") || die;
    print T1 "$tP10\t 0\n";
    print T1 "$tP10\t $maxL\n";
    close(T1);   

    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    $ylabel = ($joint)? "ratio" : "number of residues";
    print GP "set output '$statsfile'\n";
    print GP "set key right top\n";
    print GP "set nokey\n";
    print GP "set xlabel '$xlabel'\n";

     
    if ($tkfsim) {
	if ($reversible) { print GP "set title \"TKF $name [r = $rI ld = $ldI mu = $muI] reversible \n"; }
	else             { print GP "set title \"TKF $name [r = $rI ld = $ldI mu = $muI]\n"; }
   }
    else {
	if ($reversible) { print GP "set title \"E1  $name [rI = $rI ldI = $ldI muI = $muI ldE = $ldE muE = $muE  muA = $muA] reversible\n"; }
	else             { print GP "set title \"E1  $name [rI = $rI ldI = $ldI muI = $muI ldE = $ldE muE = $muE  muA = $muA]\n"; }
    }
    print GP "set ylabel '$ylabel'\n";
    #print GP "set xrange  [$subsmin:$subsmax]\n"; 
    print GP "set xrange  [$tmin:$tmax]\n"; 
    print GP "set xrange  [0:3]\n"; 
    if ($logarithmic) { #logarithmic scale
	print GP "set log x\n";
	#print GP "set xrange [0.1:500000]\n";
    }
    if ($joint) { print GP "set yrange [0:1.1]\n"; }
    else        { 
	print GP "set yrange [0:*]\n"; 
	#print GP "set yrange [0:300000]\n"; 
	#print GP "set yrange [0:12000]\n"; 
 	#print GP "set yrange [0:40000]\n"; 
 	#print GP "set yrange [0:16000]\n"; 
	#print GP "set yrange [0:50000]\n"; 
    }
 
    my $cmd;
    my $field;
    my $field2;
    my $tag;

    $cmd   = "'$tB62file'   title 'BLOSUM62' with lines ls 1, ";
    $cmd  .= "'$tB90file'   title 'BLOSUM90' with lines ls 1, ";
    $cmd  .= "'$tPAM30file' title 'PAM30'    with lines ls 1, ";
    $cmd  .= "'$tPAM10file' title 'PAM10'    with lines ls 1, ";

    if (!$logarithmic) {
	$tag = ($joint)? "<fraction of ancestral res survive>" : "<\# Ancestral res survive>";
	$k = FUNCS::gnuplot_set_style_transitions("tMM");  
	$field  = 7;
	$field2 = 11;
	$cmd .= "'$plotf' using 1:$field:$field2  title '' with yerrorbars ls 7777, "; 
	$field  = 15;
	$field2 = 19;
	$cmd .= "'$plotf' using 1:$field:$field2 title '' with yerrorbars ls 6666, "; 
	$field = 3;
	$cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, "; 
	
	$tag = ($joint)? "<fraction of occupied Insertion Blocks>" : "<\# Insertion Blocks>";
	$k = FUNCS::gnuplot_set_style_transitions("tMD"); 
	
	$field  = 8;
	$field2 = 12;
	$cmd .= "'$plotf' using 1:$field:$field2 title '' with yerrorbars ls 7777, "; 
	$field  = 16;
	$field2 = 20;
	$cmd .= "'$plotf' using 1:$field:$field2 title '' with yerrorbars ls 6666, "; 
	$field = 4;
	$cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, "; 
	
	$tag = ($joint)? "<fraction of inserted res in alignment>" : "<\# Total Inserted res>";
	$k = FUNCS::gnuplot_set_style_transitions("tMI"); 
	$field  = 9;
	$field2 = 13;
	$cmd .= "'$plotf' using 1:$field:$field2 title '' with yerrorbars ls 7777, "; 
	$field  = 17;
	$field2 = 21;
	$cmd .= "'$plotf' using 1:$field:$field2 title '' with yerrorbars ls 6666, "; 
	$field = 5;
	$cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, "; 
	
	$tag = ($joint)? "<length of descendant seq / length of ancestral seq>" : "<\# Total descendant res>";
	$k = FUNCS::gnuplot_set_style_transitions("tDD"); 
	$field  = 10;
	$field2 = 14;
	$cmd .= "'$plotf' using 1:$field:$field2 title '' with yerrorbars ls 7777, "; 
	$field  = 18;
	$field2 = 22;
	$cmd .= "'$plotf' using 1:$field:$field2 title '' with yerrorbars ls 6666, "; 
	$field = 6;
	$cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k"; 
    }
    else {
	$tag = ($joint)? "<fraction of ancestral res survive>" : "<\# Ancestral res survive>";
	$k = FUNCS::gnuplot_set_style_transitions("tMM");  
	$field  = 7;
	$field2 = 11;
	$cmd .= "'$plotf' using 1:$field:$field2  title '' with yerrorbars ls 7777, "; 
	$field = 3;
	$cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, "; 
	
	$tag = ($joint)? "<fraction of occupied Insertion Blocks>" : "<\# Insertion Blocks>";
	$k = FUNCS::gnuplot_set_style_transitions("tMD"); 
	
	$field  = 8;
	$field2 = 12;
	$cmd .= "'$plotf' using 1:$field title '' ls 7777, "; 
	$field = 4;
	$cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, "; 
	
	$tag = ($joint)? "<fraction of inserted res in alignment>" : "<\# Total Inserted res>";
	$k = FUNCS::gnuplot_set_style_transitions("tMI"); 
	$field  = 9;
	$field2 = 13;
	$cmd .= "'$plotf' using 1:$field title '' ls 7777, "; 
	$field = 5;
	$cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, "; 
	
	$tag = ($joint)? "<length of descendant seq / length of ancestral seq>" : "<\# Total descendant res>";
	$k = FUNCS::gnuplot_set_style_transitions("tDD"); 
	$field  = 10;
	$field2 = 14;
	$cmd .= "'$plotf' using 1:$field title '' ls 7777, "; 
	$field = 6;
	$cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k"; 
    }

    $cmd .= "\n";
    print GP "plot $cmd\n";

    # elementary functions
    print GP "set output '$Pfile'\n";
    print GP "set yrange [0:1]\n";
    $ylabel = "Probability";
    print GP "set ylabel '$ylabel'\n";
    print GP "unset logscale x\n";
   
    $cmd   = "'$tB62file'   title 'BLOSUM62' with lines ls 1, ";
    $cmd  .= "'$tB90file'   title 'BLOSUM90' with lines ls 1, ";
    $cmd  .= "'$tPAM30file' title 'PAM30'    with lines ls 1, ";
    $cmd  .= "'$tPAM10file' title 'PAM10'    with lines ls 1, ";

    $tag = "GAMMA";
    $k = FUNCS::gnuplot_set_style_transitions("tMM");  
    $field  = 24;
    $field2 = 25;
    $cmd .= "'$plotf' using 1:$field:$field2  title '' with yerrorbars   ls 7777, "; 
    $field  = 26;
    $field2 = 27;
    $cmd .= "'$plotf' using 1:$field:$field2  title '' with yerrorbars   ls 6666, "; 
    $field = 23;
    $cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, "; 

    $tag = "BETA";
    $k = FUNCS::gnuplot_set_style_transitions("tMD");  
    $field  = 29;
    $field2 = 30;
    $cmd .= "'$plotf' using 1:$field:$field2  title '' with yerrorbars   ls 7777, "; 
    $field  = 31;
    $field2 = 32;
    $cmd .= "'$plotf' using 1:$field:$field2  title '' with yerrorbars   ls 6666, "; 
    $field = 28;
    $cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, "; 

    $tag = "ETA";
    $k = FUNCS::gnuplot_set_style_transitions("tMI");  
    $field  = 34;
    $field2 = 35;
    $cmd .= "'$plotf' using 1:$field:$field2  title '' with yerrorbars   ls 7777, "; 
    $field  = 36;
    $field2 = 37;
    $cmd .= "'$plotf' using 1:$field:$field2  title '' with yerrorbars   ls 6666, "; 
    $field = 33;
    $cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k"; 

    $cmd .= "\n";
    print GP "plot $cmd\n";
  

    # elementary functions -- logs
    print GP "set output '$Lfile'\n";
    $ylabel = "-Log Probability";
    print GP "set ylabel '$ylabel'\n";
    print GP "set yrange [0.0:20]\n";
    print GP "set log x\n";
    print GP "set log y\n";
    print GP "set yrange [0.0000001:10]\n";
 
    $cmd   = "'$tB62file'   title 'BLOSUM62' with lines ls 1, ";
    $cmd  .= "'$tB90file'   title 'BLOSUM90' with lines ls 1, ";
    $cmd  .= "'$tPAM30file' title 'PAM30'    with lines ls 1, ";
    $cmd  .= "'$tPAM10file' title 'PAM10'    with lines ls 1, ";

    $tag = "GAMMA";
    $k = FUNCS::gnuplot_set_style_transitions("tMM");  
    $field  = 39;
    $field2 = 40;
    $cmd .= "'$plotf' using 1:$field  title '' ls 7777, "; 
    $field  = 41;
    $field2 = 42;
    $cmd .= "'$plotf' using 1:$field  title '' ls 6666, "; 
    $field = 38;
    $cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, "; 

    if (0) {
    $tag = "ETA";
    $k = FUNCS::gnuplot_set_style_transitions("tMI");  
    $field  = 49;
    $field2 = 50;
    $cmd .= "'$plotf' using 1:$field  title '' ls 7777, "; 
    $field  = 51;
    $field2 = 52;
    $cmd .= "'$plotf' using 1:$field  title '' ls 6666, "; 
    $field = 48;
    $cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k, "; 
    }

    $tag = "BETA";
    $k = FUNCS::gnuplot_set_style_transitions("tMD");  
    $field  = 44;
    $field2 = 45;
    $cmd .= "'$plotf' using 1:$field  title '' ls 7777, "; 
    $field  = 46;
    $field2 = 47;
    $cmd .= "'$plotf' using 1:$field  title '' ls 6666, "; 
    $field = 43;
    $cmd .= "'$plotf' using 1:$field  title '$tag' with lines ls $k "; 


    $cmd .= "\n";
    print GP "plot $cmd\n";
    
    close(GP);
    
    system("ps2pdf $statsfile $pdfstatsfile\n");
    system("rm $statsfile\n");
    if (!$logarithmic) { system("ps2pdf $Pfile $pdfPfile\n"); system("rm $Pfile\n");}
    else               { system("ps2pdf $Lfile $pdfLfile\n"); system("rm $Lfile\n");}
    if ($view) { 
	system("open $pdfstatsfile&\n"); 
	if (!$logarithmic) { system("open $pdfPfile&\n"); }
	else               { system("open $pdfLfile&\n"); }
    }

}


sub create_fitplot {
    my ($wrkdir, $name, $file, $reversible, $tkfsim, $view) = @_;

    my $style = 1;
    my $plot    = "$wrkdir/$name.fitplot.ps";
    my $pdfplot = "$wrkdir/$name.fitplot.pdf";

    my $observed = 1;
    my $N = 0;
    my $maxlen = -1;
    my $maxocc = -1;

    open (FILE, "$file") || die;
    while (<FILE>) {
 	if (/^\&$/) {
	    $N ++;
	    if ($N%2==0) { $observed = 1; }
	    else         { $observed = 0; }
 	}
 	elsif (/^(\S+)\s+(\S+)$/) {
	    my $len = $1+1;
	    my $occ = $2;
	    if ($observed) { if ($len > $maxlen) { $maxlen = $len; } if ($occ > $maxocc) { $maxocc = $occ; } }
	}
    }
    close(FILE); 
    printf "MAXLEN $maxlen MAXOCC $maxocc\n";

    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);

    print GP "set output '$plot'\n";
    print GP "set key right top\n";
    print GP "set nokey\n";
    print GP "set xlabel 'insert length'\n";
    print GP "set ylabel 'number of inserts'\n";

    my $plotlen;
    my $plotocc;
    print  GP "set multiplot\n";
    print  GP "set size 1,1\n";
    print  GP "set origin 0,0\n";

    $plotlen = 350;
    $plotocc = 250000;
    
    #$plotlen = 90;
    #$plotocc = 500000;
    printf GP "set xrange [%d:%d]\n", 1, $plotlen;
    printf GP "set yrange [%d:%d]\n", 0, $plotocc;
  
    #printf GP "set xrange [%d:%d]\n", 1, $maxlen;
    #printf GP "set yrange [%d:%d]\n", 0, $maxocc;

   if (1) {# for log plots of Figure S1
	print GP "set log y\n";
	printf GP "set yrange [%f:%d]\n", 0.5, $plotocc;
 	#printf GP "set yrange [0.99:*]\n";
    }

 
    $observed = 1;
    my $n = 0;
    my $choose;
    $choose = 14; # blosum62
    $choose = 18; # blosum45
    $choose = 22; # PAM240
    #$choose = 23; # PAM250
    #$choose = 30; # t=3.0
    printf "\nHISTOGRAM for n=%d\n", $choose; 
    open (FILE, "$file") || die;
    while (<FILE>) {
	if (/^\&$/) {
	    if ($n == 2*$choose-2 || $n == 2*$choose-1) { print GP "e\n"; }
	    $n ++;
	    if ($n%2==0) { $observed = 1; }
	    else         { $observed = 0; }

	    if ($n == 2*$choose-2 || $n == 2*$choose-1) {
		if ($observed) { $style ++; printf "OBSERVED\n"; printf GP "plot '-' using 1:2 with points ls %d\n", $style; }
		else           {            printf "EXPECTED\n"; printf GP "plot '-' using 1:2 with lines  ls 888\n"; }
	    }
	}
 	elsif (/^(\S+)\s+(\S+)$/) {
	    if ($n == 2*$choose-2 || $n == 2*$choose-1) { 
		printf GP "%f %f\n", $1+2, $2; 
		printf "%f %f\n", $1+2, $2; 
	    }
	}
   }
    close(FILE); 
    close(GP); 
    
    system("ps2pdf $plot $pdfplot\n");
    system("rm $plot\n");
    if ($view) { 
	system("open $pdfplot&\n"); 
    }
}




sub create_plot {
    my ($wrkdir, $name, $file, $reversible, $tkfsim, $view) = @_;

    my $plotf       = "$wrkdir/$name.plot";
    my $plotf_tstar = "$wrkdir/$name.plot_tstar";
    
    my $N;
    my $L;
    my $joint;
    my $tmin;
    my $tmax;
    my $tinterval;
    my $tinc;
    my $teps;

    my $rI;
    my $ldE;
    my $muE;
    my $ldI;
    my $muI;
    my $muA;

    my $sE_B62 = "";
    my $nE_B62 = "";
    my $eE_B62 = ""; 
    my $sE_B90 = "";
    my $nE_B90 = "";
    my $eE_B90 = ""; 
    my $sE_P30 = "";
    my $nE_P30 = "";
    my $eE_P30 = ""; 
    my $sE_P10 = "";
    my $nE_P10 = "";
    my $eE_P10 = "";

    my $tstar1 = parse_outfile($plotf,       $file, 0, \$L, \$N, \$joint, \$tmin, \$tmax, \$tinterval, \$tinc, \$teps, \$rI, \$ldE, \$muE, \$ldI, \$muI, \$muA, 
			       \$sE_B62, \$nE_B62, \$eE_B62, \$sE_B90, \$nE_B90, \$eE_B90, 
			       \$sE_P30, \$nE_P30, \$eE_P30, \$sE_P10, \$nE_P10, \$eE_P10); 
    my $tstar  = parse_outfile($plotf_tstar, $file, 1); 

    #if (abs($tstar1-$tstar) > 1e-6) { print "bad tstar\n"; die; }
    actual_plot($name, $tstar, $plotf,  $plotf_tstar, $L, $N, $joint, $reversible, $tkfsim,
		$tmin, $tmax, $tinterval, $tinc, $teps, $rI, $ldE, $muE, $ldI, $muI, $muA, 
		$sE_B62, $nE_B62, $eE_B62, $sE_B90, $nE_B90, $eE_B90, 
		$sE_P30, $nE_P30, $eE_P30, $sE_P10, $nE_P10, $eE_P10, $view); 
}

sub parse_outfile {

    my ($plotfile, $file, $ist1, $ret_L, $ret_N, $ret_joint, $ret_tmin, $ret_tmax, $ret_tinterval, $ret_tinc, $ret_teps, 
	$ret_rI, $ret_ldE, $ret_muE, $ret_ldI, $ret_muI, $ret_muA, 
	$ret_sE_B62, $ret_nE_B62, $ret_eE_B62, $ret_sE_B90, $ret_nE_B90, $ret_eE_B90, 
	$ret_sE_P30, $ret_nE_P30, $ret_eE_P30, $ret_sE_P10, $ret_nE_P10, $ret_eE_P10
	) = @_;

    print "cannot find file $plotfile\n"; 
    open(OUT, ">$plotfile") ||die;

    my $L;
    my $N;
    my $joint;
    my $p;
    my $tmin;
    my $tmax;
    my $tinterval;
    my $tinc;
    my $teps;

    my $tstar = 0.0;
    my $rI = 0.0;

    my $muA = -1;
    my $muI = -1;
    my $muE = -1;
    my $ldI = -1;
    my $ldE  = -1;

    my $sE_B62;
    my $nE_B62;
    my $eE_B62; 
    my $sE_B90;
    my $nE_B90;
    my $eE_B90; 
    my $sE_P30;
    my $nE_P30;
    my $eE_P30; 
    my $sE_P10;
    my $nE_P10;
    my $eE_P10;

    my $time    = -1;
    my $len;

    my $sE_the;
    my $nE_the;
    my $eE_the;
    my $lE_the;
  
    my $sE_fin;
    my $nE_fin;
    my $eE_fin;
    my $lE_fin;
    my $sV_fin;
    my $nV_fin;
    my $eV_fin;
    my $lV_fin;

    my $sE_inf;
    my $nE_inf;
    my $eE_inf;
    my $lE_inf;
    my $sV_inf;
    my $nV_inf;
    my $eV_inf;
    my $lV_inf;

    my $gammaE_the;
    my $gammaE_fin;
    my $gammaV_fin;
    my $gammaE_inf;
    my $gammaV_inf;
    
    my $betaE_the;
    my $betaE_fin;
    my $betaV_fin;
    my $betaE_inf;
    my $betaV_inf;
    
    my $etaE_the;
    my $etaE_fin;
    my $etaV_fin;
    my $etaE_inf;
    my $etaV_inf;
    
    my $loggammaE_the;
    my $loggammaE_fin;
    my $loggammaV_fin;
    my $loggammaE_inf;
    my $loggammaV_inf;
    
    my $logbetaE_the;
    my $logbetaE_fin;
    my $logbetaV_fin;
    my $logbetaE_inf;
    my $logbetaV_inf;
    
    my $logetaE_the;
    my $logetaE_fin;
    my $logetaV_fin;
    my $logetaE_inf;
    my $logetaV_inf;
    
    open (FILE, "$file") || die;
    while (<FILE>) {
	
	if (/^# L\s+=\s+(\S+)\s+/) {
	    $L = $1;
	    $joint = 0;
	    $p = $L/($L+1);
	}
	elsif (/^# p\s+=\s+(\S+)\s+/) {
	    $p = $1;
	    $joint = 1;
	    $L = $p/(1.-$p);
	}
	elsif (/^# N\s+=\s+(\S+)\s+/) {
	    $N = $1;
	}
	elsif (/^# tmin\s+=\s+(\S+)\s+/) {
	    $tmin = $1;
	}
	elsif (/^# tmax\s+=\s+(\S+)\s+/) {
	    $tmax = $1;
	}
	elsif (/^# tinterval\s+=\s+(\S+)\s+/) {
	    $tinterval = $1;
	}
	elsif (/^# tinc\s+=\s+(\S+)\s+/) {
	    $tinc = $1;
	}
	elsif (/^# teps\s+=\s+(\S+)\s+/) {
	    $teps = $1;
	}
	elsif (/^# rI:\s+=\s+(\S+)\s+/) {
	    $rI = $1;
	}
	elsif (/^\#\s*ldE\s+=\s+(\S+)/) { 
	    $ldE = $1;
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
	elsif (/^\#\s*muA\s+=\s+(\S+)/) { 
	    $muA = $1;
	}
	elsif (/^\#\s*r\s+=\s+(\S+)/) { 
	    $rI = $1;  # tkfsim 
	}
	elsif (/^\#\s*ld_tkf\s+=\s+(\S+)/) { 
	    $ldI = $1;  # tkfsim 
	}
	elsif (/^\#\s*mu_tkf\s+=\s+(\S+)/) { 
	    $muI = $1;  # tkfsim 
	}
	elsif (/^\#/) { next; 
	}
	elsif (/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*/) {
	    $time     = $1;
	    $len      = $2;

	    $sE_the   = $3;
	    $nE_the   = $8;
	    $eE_the   = $13;
	    $lE_the   = $18;

	    $sE_fin   = $4;
	    $nE_fin   = $9;
	    $eE_fin   = $14;
	    $lE_fin   = $19;
	    $sV_fin   = $5;
	    $nV_fin   = $10;
	    $eV_fin   = $15;
	    $lV_fin   = $20;

	    $sE_inf   = $6;
	    $nE_inf   = $11;
	    $eE_inf   = $16;
	    $lE_inf   = $21;
	    $sV_inf   = $7;
	    $nV_inf   = $12;
	    $eV_inf   = $17;
	    $lV_inf   = $22;

	    $gammaE_the = $23;
	    $gammaE_fin = $24;
	    $gammaV_fin = $25;
	    $gammaE_inf = $26;
	    $gammaV_inf = $27;

	    $betaE_the = $28;
	    $betaE_fin = $29;
	    $betaV_fin = $30;
	    $betaE_inf = $31;
	    $betaV_inf = $32;

	    $etaE_the = $33;
	    $etaE_fin = $34;
	    $etaV_fin = $35;
	    $etaE_inf = $36;
	    $etaV_inf = $37;

	    $loggammaE_the = ($gammaE_the > 0)? -log($gammaE_the) : 50;
	    $loggammaE_fin = ($gammaE_fin > 0)? -log($gammaE_fin) : 50;
	    $loggammaV_fin = ($gammaV_fin > 0)? -log($gammaV_fin) : 50;
	    $loggammaE_inf = ($gammaE_inf > 0)? -log($gammaE_inf) : 50;
	    $loggammaV_inf = ($gammaV_inf > 0)? -log($gammaV_inf) : 50;

	    $logbetaE_the = ($betaE_the > 0)? -log($betaE_the) : 50;
	    $logbetaE_fin = ($betaE_fin > 0)? -log($betaE_fin) : 50;
	    $logbetaV_fin = ($betaV_fin > 0)? -log($betaV_fin) : 50;
	    $logbetaE_inf = ($betaE_inf > 0)? -log($betaE_inf) : 50;
	    $logbetaV_inf = ($betaV_inf > 0)? -log($betaV_inf) : 50;

	    $logetaE_the = ($etaE_the > 0)? -log($etaE_the) : 50;
	    $logetaE_fin = ($etaE_fin > 0)? -log($etaE_fin) : 50;
	    $logetaV_fin = ($etaV_fin > 0)? -log($etaV_fin) : 50;
	    $logetaE_inf = ($etaE_inf > 0)? -log($etaE_inf) : 50;
	    $logetaV_inf = ($etaV_inf > 0)? -log($etaV_inf) : 50;

	    if ($time == 1.0) { $tstar = 1.0; }

	    if ($time <= $tB62) { $sE_B62 = $sE_the; $nE_B62 = $nE_the; $eE_B62 = $eE_the; }
	    if ($time <= $tB90) { $sE_B90 = $sE_the; $nE_B90 = $nE_the; $eE_B90 = $eE_the; }
	    if ($time <= $tP30) { $sE_P30 = $sE_the; $nE_P30 = $nE_the; $eE_P30 = $eE_the; }
	    if ($time <= $tP10) { $sE_P10 = $sE_the; $nE_P10 = $nE_the; $eE_P10 = $eE_the; }

	    if (!$ist1) {
		print OUT "$time\t$len\t$sE_the\t$nE_the\t$eE_the\t$lE_the\t$sE_fin\t$nE_fin\t$eE_fin\t$lE_fin\t$sV_fin\t$nV_fin\t$eV_fin\t$lV_fin\t$sE_inf\t$nE_inf\t$eE_inf\t$lE_inf\t$sV_inf\t$nV_inf\t$eV_inf\t$lV_inf\t$gammaE_the\t$gammaE_fin\t$gammaV_fin\t$gammaE_inf\t$gammaV_inf\t$betaE_the\t$betaE_fin\t$betaV_fin\t$betaE_inf\t$betaV_inf\t$etaE_the\t$etaE_fin\t$etaV_fin\t$etaE_inf\t$etaV_inf\t$loggammaE_the\t$loggammaE_fin\t$loggammaV_fin\t$loggammaE_inf\t$loggammaV_inf\t$logbetaE_the\t$logbetaE_fin\t$logbetaV_fin\t$logbetaE_inf\t$logbetaV_inf\t$logetaE_the\t$logetaE_fin\t$logetaV_fin\t$logetaE_inf\t$logetaV_inf\n";
	    }
	    else {
		if ($time == 1.0) {
		    print OUT "$time\t$len\t$sE_the\t$nE_the\t$eE_the\n";
		}
	    }
	}
   }
    close (FILE);
    close (OUT);

    $$ret_L         = $L;
    $$ret_N         = $N;
    $$ret_joint     = $joint;
    $$ret_tmin      = $tmin;
    $$ret_tmax      = $tmax;
    $$ret_tinterval = $tinterval;
    $$ret_tinc      = $tinc;
    $$ret_teps      = $teps;

    $$ret_rI     = $rI;
    $$ret_ldE    = $ldE;
    $$ret_muE    = $muE;
    $$ret_ldI    = $ldI;
    $$ret_muI    = $muI;
    $$ret_muA    = $muA;

    $$ret_sE_B62  = $sE_B62;
    $$ret_nE_B62  = $nE_B62;
    $$ret_eE_B62  = $eE_B62;

    $$ret_sE_B90  = $sE_B90;
    $$ret_nE_B90  = $nE_B90;
    $$ret_eE_B90  = $eE_B90;

    $$ret_sE_P30  = $sE_P30;
    $$ret_nE_P30  = $nE_P30;
    $$ret_eE_P30  = $eE_P30;

    $$ret_sE_P10  = $sE_P10;
    $$ret_nE_P10  = $nE_P10;
    $$ret_eE_P10  = $eE_P10;

    return $tstar;
}


sub plot_program {
    my ($wrkdir, $reversible, $tkfsim, $view) = @_;

    my $name = $wrkdir;
    if ($name =~ /^\S+\/([^\/]+)$/) { $name = $1; }

    my $file = "$wrkdir/$name.out";

    create_plot($wrkdir, $name, $file, $reversible, $tkfsim, $view);
}

sub plot_geofit {
    my ($wrkdir, $reversible, $tkfsim, $view) = @_;

    my $name = $wrkdir;
    if ($name =~ /^\S+\/([^\/]+)$/) { $name = $1; }

    my $file = "$wrkdir/$name.fit";

    create_fitplot($wrkdir, $name, $file, $reversible, $tkfsim, $view);
}




sub run_e1sim {
    my ($program, $wrkdir, $filename, $L, $N, $tmin, $tmax, $tinterval, $tinc, $teps, $evomodel, $rI, $ldE, $muE, $ldI, $muI, $muA, $sI, $sD, $vI, $vD, $joint, $reversible, $geofit, $verbose) = @_;
    
    my $dumpfile = "$wrkdir/$filename.dump";
    my $outfile = "$wrkdir/$filename.out";
    my $fitfile = "$wrkdir/$filename.fit";
    
    my $cmd = "$program --evomodel $evomodel --rI $rI --ldEM $ldE --ldI $ldI --muEM $muE --muI $muI --muAM $muA --sI $sI --sD $sD --vI $vI --vD $vD -N $N -L $L --tmin $tmin --tmax $tmax --tinterval $tinterval --tinc $tinc --teps $teps  ";  
    
    if ($verbose)    { $cmd .= " -v "; }
    if ($joint)      { $cmd .= " --joint "; }
    if ($reversible) { $cmd .= " --rev "; }
    if ($geofit)     { $cmd .= " --geofit $fitfile "; }
    if ($opt_s) { $cmd .= " --seed $opt_s "; }
    system("echo $cmd $outfile \n");
    system("$cmd $outfile > $dumpfile\n");
}

sub run_tkfsim {
    my ($program, $wrkdir, $filename, $L, $N, $tmin, $tmax, $tinterval, $tinc, $teps, $rI, $ld, $mu, $joint, $reversible, $geofit, $verbose) = @_;
    
    my $dumpfile = "$wrkdir/$filename.dump";
    my $outfile  = "$wrkdir/$filename.out";
    my $fitfile  = "$wrkdir/$filename.fit";
    
    my $cmd = "$program --evomodel $evomodel --rI $rI --ld $ld --mu $mu --sI $sI --sD $sD --vI $vI --vD $vD -N $N -L $L --tmin $tmin --tmax $tmax --tinterval $tinterval --tinc $tinc --teps $teps  ";  
    
    if ($verbose)    { $cmd .= " -v "; }
    if ($joint)      { $cmd .= " --joint "; }
    if ($reversible) { $cmd .= " --rev "; }
    if ($geofit)     { $cmd .= " --geofit $fitfile "; }
    if ($opt_s) { $cmd .= " --seed $opt_s "; }
    system("echo $cmd $outfile\n");
    system("$cmd $outfile > $dumpfile\n");
}
