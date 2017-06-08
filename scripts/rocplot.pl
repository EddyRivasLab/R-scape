#!/usr/bin/perl -w
# rocplot.pl

use strict;
use Class::Struct;

use vars qw ($opt_D $opt_L $opt_P $opt_v);  # required if strict used
use Getopt::Std;
getopts ('D:L:P:v');


# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  rocplot.pl [options] <F> <file1>..<fileF> <afafile> <rscapedir> <gnuplotdir> \n\n";
        print "options:\n";
 	exit;
}

my $F = shift;
my @prefile;
my @prename;
my @rocfile;
my $ID = 1.0;
my $outname;
 for (my $f = 0; $f < $F; $f ++){
    $prefile[$f]  = shift;
    $prename[$f] = $prefile[$f];
    if ($prename[$f] =~ /^(\S+)\.sorted.out$/) {
	$outname = $1;
 
	$prename[$f] = "$outname.rscape";
	if ($outname     =~ /\/([^\/]+)\s*$/) { $outname = $1; }
	if ($prename[$f] =~ /\/(ID\S+)\//)    { $ID = $1; $outname .= ".$ID"; }
    }
}
my $afafile  = shift;
my $pfamname = $prefile[0];
if ($pfamname =~ /(PF[^\.]+)\./) { $pfamname = $1; }
if ($pfamname =~ /(RF[^\.]+)\./) { $pfamname = $1; }
my $oafafile = "data/$pfamname.afa";
print "original file: $oafafile\n";

my $rscapedir = shift;
my $gnuplot   = shift;
use constant GNUPLOT => '$gnuplot';
use constant GNUPLOT => '/usr/local/bin/gnuplot';
use lib '$rscapedir/scripts';
#use FUNCS;

my $currdir = $ENV{PWD};

my $maxD = 8;
if ($opt_D) { $maxD = $opt_D; }

my $minL = -1;
if ($opt_L) { $minL = $opt_L; }
for (my $f = 0; $f < $F; $f ++) {  $rocfile[$f] = ($minL>0)? "$prename[$f].minL$minL.roc":"$prename[$f].roc"; }

my $pdbfile = "";
my $ncnt = 0;
my @cnt;
if ($opt_P) { 
    $pdbfile = "$opt_P";
    contactlist_from_pdbfile($pdbfile, $afafile, $maxD, $minL, \$ncnt, \@cnt);
}

my $verbose = 0;
if ($opt_v) { $verbose = 1; }

# add a random file
my $dorandom = 0;
if ($dorandom) {
    if ($prename[0] =~ /results\/[^\/]+_([^\/]+)\//) { 
	$rocfile[$F] = "results/random_$1/$pfamname.random.minL$minL.roc";
    }
    else { 
	$rocfile[$F] = "results/random/$pfamname.random.minL$minL.roc";
    }
    $F ++;
}


for (my $f = 0; $f < $F; $f ++) {
    my$method = "";
    
    if    ($prefile[$f] =~ /results\/(\S+)_filtered\//) { $method = $1; }
    elsif ($prefile[$f] =~ /results\/([^\/]+)\//)       { $method = $1; }

    print "\n$method: $prefile[$f]\n";
    if ($method =~ /^R-scape$/) {
	create_rocfile_rscape($rocfile[$f], $prefile[$f]);
    }
    elsif ($method =~ /^random$/) {
	create_rocfile_random($rocfile[$f], $prefile[$f]);
    }
    elsif ($method =~ /^DCA$/) {
	create_rocfile_DCA($rocfile[$f], $prefile[$f]);
    }
    else { print "method $method not implemented yet\n"; die; }
    
   
}

my $xmax = 100;
my $viewplots = 1;
rocplot($pfamname, $F, \@rocfile, \@prename, $xmax, $viewplots);


sub  create_rocfile_rscape {
    my ($rocfile, $file) = @_;

    open my $fp, '>', $rocfile || die "Can't open $rocfile: $!";

    my $f   = 0;
    my $f_c = 0;
    my $f_b = 0;
    my $f_w = 0;
    my $t_c = 0;
    my $t_b = 0;
    my $t_w = 0;

    my $alen = 0;
    my $type = "";
    open(FILE, "$file") || die;
    while(<FILE>) {
	if (/^\# contacts\s+(\d+)\s+\((\d+)\s+bpairs\s+(\d+)\s+wc/) {
	    $t_c = $1;
	    $t_b = $2;
	    $t_w = $3;
	}
	elsif (/^\#.+alen\s+(\d+)\s+/) {
	    $alen = $1;
	}
	elsif (/^\#/) {
	}
	elsif (/^(\S+)\s+\d+\s+\d+\s+\S+\s+\S+\s*$/) {
	    $type = $1;
	    
	    if    ($type =~ /^\*$/)   {  # A WWc basepair
		$f_w ++;
		$f_b ++;
		$f_c ++;
	    }
	    elsif ($type =~ /^\*\*$/) {  # A basepair
		$f_b ++;
		$f_c ++;
	    }
	    elsif ($type =~ /^c/)     {  # A contact
		$f_c ++;
	    }
	    elsif ($type =~ /^\~$/)    { # Not a contact although compatible with the bpairs
	    }
	    else { print "type? $type\n"; die; }
	    
	    $f ++;
	    writeline($fp, $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $alen);
	}
 	elsif (/\s+\d+\s+\d+\s+\S+\s+\S+\s*$/) { # Not annotated as a contact
	    $f   ++;
	    writeline($fp, $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $alen);	    
	}
    }
    close(FILE);
    close($fp);
}

sub writeline {
    my ($fp, $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $alen) = @_;

    my $sen_c, my $sen_b, my $sen_w;
    my $ppv_c, my $ppv_b, my $ppv_w;
    my $F_c,   my $F_b,   my $F_w;
    
    calculateF($f_c, $t_c, $f, \$sen_c, \$ppv_c, \$F_c);
    calculateF($f_b, $t_b, $f, \$sen_b, \$ppv_b, \$F_b);
    calculateF($f_w, $t_w, $f, \$sen_w, \$ppv_w, \$F_w);
    
    printf $fp "$f\t$f_c\t$f_b\t$f_w\t\t$t_c\t$t_b\t$t_w\t$alen\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
    $f/$alen, $f_c/$alen, $f_b/$alen, $f_w/$alen, $t_c/$alen, $t_b/$alen, $t_w/$alen, $sen_c, $ppv_c, $F_c, $sen_b, $ppv_b, $F_b, $sen_w, $ppv_w, $F_w;

    #printf     "$f\t$f_c\t$f_b\t$f_w\t\t$t_c\t$t_b\t$t_w\t$alen\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
    #$f/$alen, $f_c/$alen, $f_b/$alen, $f_w/$alen, $t_c/$alen, $t_b/$alen, $t_w/$alen, $sen_c, $ppv_c, $F_c, $sen_b, $ppv_b, $F_b, $sen_w, $ppv_w, $F_w;
}

sub calculateF {
    my ($tph, $th, $fh, $ret_sen, $ret_ppv, $ret_F) = @_;
    
    my $sen;
    my $ppv;
    my $F;
    
    $sen = ($th > 0)? $tph/$th : 0.0;
    $ppv = ($fh > 0)? $tph/$fh : 0.0;
    
    if ($th+$fh > 0.0) { $F = 2.0 * $tph / ($th + $fh) }
    else               { $F = 0.0 }

    $$ret_sen = 100.*$sen;
    $$ret_ppv = 100.*$ppv;
    $$ret_F   = 100.*$F;
}

sub  create_rocfile_random {
    my ($rocfile, $prefile) = @_;
}

sub  create_rocfile_DCA {
    my ($rocfile, $prefile) = @_;

    my $method = "mfDCA";
    my $which  = "DI";
    if ($rocfile =~ /mfDCA\.MI/) {
	$which = "MI";
	$method .= "-$which";
    }
    
}


sub rocplot {
    my ($pfamname, $F, $file_ref, $prename_ref, $xmax, $seeplots) = @_;


   my $psfile = "results/$pfamname.N$F.ps";
    
    #if ($psfile =~ /\/([^\/]+)\s*$/) { $psfile = "$1"; }
    my $pdffile = $psfile;
    if ($pdffile =~ /^(\S+).ps$/) { $pdffile = "$1.pdf"; }
    print "FILE: $psfile\n";

    my $maxpp = 2;

    my $xlabel;
    my $ylabel;
    my $title  = "$pfamname";
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

    my $logscale = 0;
    $xlabel = "SEN contacts";
    $ylabel = "PPV contacts";
    $x_min = 0;
    $x_max = 100;
    $y_min = 0;
    $y_max = 100;
    $x = 16;
    $y = 17;
    oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "SEN bpairs";
    $ylabel = "PPV bpairs";
    $x_min = 0;
    $x_max = 100;
    $y_min = 0;
    $y_max = 100;
    $x = 19;
    $y = 20;
    oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);

    $logscale = 0;
    $xlabel = "number of predictions per position";
    $ylabel = "PPV contacts";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 9;
    $y = 17;
    oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "SEN contacts";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 9;
    $y = 16;
    oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "F contacts";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 9;
    $y = 18;
    oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);

    # basepairs
    $xlabel = "number of predictions per position";
    $ylabel = "PPV bpairs";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 9;
    $y = 20;
    oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "SEN bpairs";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 9;
    $y = 19;
    oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "F bpairs";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 9;
    $y = 21;
    oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
 
    $xlabel = "number of predictions";
    $ylabel = "PPV contacts";
    $x_min = 0.001;
    $x_max = $xmax;
    $y_min = 0;
    $y_max = 100;
    $x = 1;
    $y = 17;
    oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions";
    $ylabel = "SEN contacts";
    $x_min = 0.001;
    $x_max = $xmax;
    $y_min = 0;
    $y_max = 100;
    $x = 1;
    $y = 16;
    oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions";
    $ylabel = "F contacts";
    $x_min = 0.001;
    $x_max = $xmax;
    $y_min = 0;
    $y_max = 100;
    $x = 1;
    $y = 18;
    oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);


    close($gp);

    #system ("/usr/local/bin/ps2pdf $psfile $pdffile\n"); 
    #system("rm $psfile\n");
    if ($seeplots) { system ("open $psfile&\n"); }

    
}


sub oneplot {
    my ($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $xmin, $xmax, $ymin, $ymax, $logscale) = @_;
   
    my $cmd = "";
    my $m = 5;
    
    print $gp "set title  '$title'\n";
    print $gp "set xlabel '$xlabel'\n";
    print $gp "set ylabel '$ylabel'\n";
    print $gp "set xrange [$xmin:$xmax]\n";
    print $gp "set yrange [$ymin:$ymax]\n";
    if ($logscale) { print $gp "set logscale x\n"; }
    for (my $f = 0; $f < $F; $f++) {
	my $key = $prename_ref->[$f];
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using $x:$y  title ''                ls $m, " : "'$file_ref->[$f]' using $x:$y  title ''                ls $m, ";
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using $x:$y  title '$key' with lines ls $m"   : "'$file_ref->[$f]' using $x:$y  title '$key' with lines ls $m, ";
	$m ++; if ($m == 7) { $m = 8; }
    }
    print $gp "plot $cmd\n";
    if ($logscale) { print $gp "unset logscale\n"; }
}


sub contactlist_from_pdbfile 
{
    my ($pdbfile, $afafile, $maxD, $minL, $ret_ncnt, $cnt_ref) = @_;

    my $ncnt = 0;

    $$ret_ncnt = $ncnt;
}
