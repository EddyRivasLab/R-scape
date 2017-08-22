#!/usr/bin/perl -w
# rocplot_together.pl

use strict;
use Class::Struct;

# find directory where the script is installed
use FindBin;
use lib $FindBin::Bin;
use PDBFUNCS;
use FUNCS;

use vars qw ($opt_C $opt_P $opt_v);  # required if strict used
use Getopt::Std;
getopts ('CPv');


# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  rocplot_together.pl [options] <DIR> <string_name> <string_type> <string_suffix> <gnuplot> \n\n";
        print "options:\n";
 	exit;
}

my $DIR           = shift;
my $string_name   = shift;
my $string_type   = shift;
my $string_suffix = shift;
my $gnuplot       = shift;
 
my @type = split(/\s+/, $string_type);
my $M = $#type+1;
for (my $m = 0; $m < $M; $m++)
{
    $type[$m] =~ s/ //g;
}

my $famtype = "ALL";
if ($opt_P) { $famtype = "PFAM"; }
if ($opt_C) { $famtype = "CAMEO"; }

my $seeplots = 0;
my $verbose  = 0;

my $N = 2;
my $k = 50;
my $shift = 0;

my @plotfile;

for (my $m = 0; $m < $M; $m++) {
    print "\n$type[$m]\n";

    my $type = $type[$m];
    $type =~ s/\//\_/g;
    if ($type =~ /^R-scape_(\S+)$/) { $type = $1; }
    
    $plotfile[$m] = "$DIR/results/$string_name.$type.$string_suffix.plot";

    my @his_f;
    my @his_fc;
    my @his_fb;
    my @his_fw;
    my @his_tc;
    my @his_tb;
    my @his_tw;
    FUNCS::init_histo_array($N, $k, \@his_f);
    FUNCS::init_histo_array($N, $k, \@his_fc);
    FUNCS::init_histo_array($N, $k, \@his_fb);
    FUNCS::init_histo_array($N, $k, \@his_fw);
    FUNCS::init_histo_array($N, $k, \@his_tc);
    FUNCS::init_histo_array($N, $k, \@his_tb);
    FUNCS::init_histo_array($N, $k, \@his_tw);

    my $localdir = "$DIR/results/$type[$m]";
    my @family;
    FUNCS::sorted_files($localdir, \@family, $string_suffix);    
    my $F = $#family+1;
    
    my $nf = 0;
    for (my $f = 0; $f < $F; $f++)
    {
	my $rocfile = "$localdir/$family[$f]";

	my $add = ($famtype =~ /^ALL$/)? 1 : 0;
	if ($famtype =~ /^CAMEO$/ && $rocfile =~ /\/\d\S+$/)   { print "++$rocfile\n"; $add = 1; }
	if ($famtype =~ /^PFAM$/  && $rocfile =~ /\/PF\S+$/)   { $add = 1; }
	print "^^famtype $famtype roc $rocfile add $add\n";
	
	if ($add == 0) { next; }

	$nf ++;
	
	open (FILE, "$rocfile") || print "\nFILE NOT FOUND\n";
	while(<FILE>) {
	    
	    if (/\#/) {
	    }
	    elsif (/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+/) {
		my $f   = $1;
		my $fc  = $2;
		my $fb  = $3;
		my $fw  = $4;
		my $tc  = $5;
		my $tb  = $6;
		my $tw  = $7;
		my $fpp = $8; # predictions per position

		if ($fpp <= $N) {
		    FUNCS::fill_histo_array($f,       $fpp, $N, $k, $shift, \@his_f);
		    FUNCS::fill_histo_array($fc,      $fpp, $N, $k, $shift, \@his_fc);
		    FUNCS::fill_histo_array($fb,      $fpp, $N, $k, $shift, \@his_fb);
		    FUNCS::fill_histo_array($fw,      $fpp, $N, $k, $shift, \@his_fw);
		    FUNCS::fill_histo_array($tc,      $fpp, $N, $k, $shift, \@his_tc);
		    FUNCS::fill_histo_array($tb,      $fpp, $N, $k, $shift, \@his_tb);
		    FUNCS::fill_histo_array($tw,      $fpp, $N, $k, $shift, \@his_tw);
		}
		else { last; }
	    }
	}
	close(FILE);	
    }
    print "$nf/$F MSA $plotfile[$m]\n";

    open (PLOT, ">$plotfile[$m]") || die;
    for (my $i = 0; $i < $N*$k; $i++) {
	my $fpp = $i / $k;
	my $sen_c, my $ppv_c, my $F_c;
	my $sen_b, my $ppv_b, my $F_b;
	my $sen_w, my $ppv_w, my $F_w;
	FUNCS::calculateF($his_fc[$i], $his_tc[$i], $his_f[$i], \$sen_c, \$ppv_c, \$F_c);
	FUNCS::calculateF($his_fb[$i], $his_tb[$i], $his_f[$i], \$sen_b, \$ppv_b, \$F_b);
	FUNCS::calculateF($his_fw[$i], $his_tw[$i], $his_f[$i], \$sen_w, \$ppv_w, \$F_w);
	print PLOT "$fpp\t$sen_c\t$ppv_c\t$F_c\t$sen_b\t$ppv_b\t$F_b\t$sen_w\t$ppv_w\t$F_w\n";
    }
    close(PLOT);
    
}

rocplot($M, \@plotfile, \@type, $gnuplot, $seeplots);

sub rocplot {
    my ($F, $file_ref, $type_ref, $gnuplot, $seeplots) = @_;


   my $psfile = "$string_name.$string_suffix.ps";
    
    #if ($psfile =~ /\/([^\/]+)\s*$/) { $psfile = "$1"; }
    my $pdffile = $psfile;
    if ($pdffile =~ /^(\S+).ps$/) { $pdffile = "$1.pdf"; }
    print "FILE: $psfile\n";

    my $maxpp = 2.0;

    my $xlabel;
    my $ylabel;
    my $title  = "$string_name";
    my $x;
    my $y;
    my $x_max, my $x_min;
    my $y_max, my $y_min;
    
    open(GP,'|'."$gnuplot") || die "Gnuplot: $!";
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
    $x_max = 40;
    $y_min = 0;
    $y_max = 100;
    $x = 2;
    $y = 3;
    oneplot($gp, $F, $file_ref, $type_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "SEN bpairs";
    $ylabel = "PPV bpairs";
    $x_min = 0;
    $x_max = 100;
    $y_min = 0;
    $y_max = 100;
    $x = 5;
    $y = 6;
    oneplot($gp, $F, $file_ref, $type_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);

    $logscale = 0;
    $xlabel = "number of predictions per position";
    $ylabel = "PPV contacts";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 1;
    $y = 3;
    oneplot($gp, $F, $file_ref, $type_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "SEN contacts";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 1;
    $y = 2;
    oneplot($gp, $F, $file_ref, $type_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "F contacts";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 1;
    $y = 4;
    oneplot($gp, $F, $file_ref, $type_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);

    # basepairs
    $xlabel = "number of predictions per position";
    $ylabel = "PPV bpairs";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 1;
    $y = 6;
    oneplot($gp, $F, $file_ref, $type_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "SEN bpairs";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 1;
    $y = 5;
    oneplot($gp, $F, $file_ref, $type_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "F bpairs";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 1;
    $y = 7;
    oneplot($gp, $F, $file_ref, $type_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
 
    
    if ($seeplots) { system ("open $psfile&\n"); }

    
}


sub oneplot {
    my ($gp, $F, $file_ref, $type_ref, $x, $y, $xlabel, $ylabel, $title, $xmin, $xmax, $ymin, $ymax, $logscale) = @_;
   
    my $cmd = "";
    my $m = 5;
    
    print $gp "set title  '$title'\n";
    print $gp "set xlabel '$xlabel'\n";
    print $gp "set ylabel '$ylabel'\n";
    print $gp "set xrange [$xmin:$xmax]\n";
    print $gp "set yrange [$ymin:$ymax]\n";
    if ($logscale) { print $gp "set logscale x\n"; }
    for (my $f = 0; $f < $F; $f++) {
	my $key = $type_ref->[$f];
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using $x:$y  title ''                ls $m, " : "'$file_ref->[$f]' using $x:$y  title ''                ls $m, ";
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using $x:$y  title '$key' with lines ls $m"   : "'$file_ref->[$f]' using $x:$y  title '$key' with lines ls $m, ";
	$m ++; if ($m == 7) { $m = 8; }
    }
    print $gp "plot $cmd\n";
    if ($logscale) { print $gp "unset logscale\n"; }
 }
