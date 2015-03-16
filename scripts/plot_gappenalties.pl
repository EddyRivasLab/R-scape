#!/usr/bin/perl -w
#plot_evohmm.pl

use strict;
use Class::Struct;
use lib '/groups/eddy/home/rivase/projects/evohmm/scripts';
use lib '/Users/rivase/projects/evohmm/scripts';
use FUNCS;
#use constant GNUPLOT => '/usr/bin/gnuplot';
use constant GNUPLOT => '/opt/local/bin/gnuplot';

use vars qw ($opt_v );  # required if strict used
use Getopt::Std;
getopts ('k:I:l:L:u:U:gnpr:st:vx:');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  plot_evohmm.pl [options] <file1> <file1> \n\n";
        print "options:\n";
 	print "-v    :  be verbose\n";
 	exit;
}
my $file1  = shift;
my $file2  = shift;

my $sP10  = 100-90.0;
my $sP30  = 100-72.9;
my $sP70  = 100-50.2;
my $sB90  = 100-44.2;
my $sB80  = 100-39.1;
my $sB62  = 100-29.8;
my $sP120 = 100-25.1;
my $sB45  = 100-23.7;

my $tP10  = 0.074;
my $tP30  = 0.229;
my $tP70  = 0.52;
my $tB90  = 0.633;
my $tB80  = 0.741;
my $tB62  = 1.0;
my $tP120 = 1.18;
my $tB45  = 1.245;


my $verbose = 0;
if ($opt_v) { $verbose = 1; }

my $viewplot = 1;

actual_plot($file1, $file2, $viewplot);

sub actual_plot{
    my ($file1, $file2, $viewplot) = @_;

    my $xlabel;
    my $ylabel = "score in half bits";

    my $ymin = -25;
    my $xmax = 3;

    my $plot = $file1; if ($plot =~ /^(\S+).out$/) { $plot = $1; }
    my $tfile = "$plot\_time.ps";
    my $sfile = "$plot\_subs.ps";
    my $ifile = "$plot\_id.ps";
    my $pdf_tfile = "$plot\_time.pdf";
    my $pdf_sfile = "$plot\_subs.pdf";
    my $pdf_ifile = "$plot\_id.pdf";

    my $WPfile   = "$plot.WP";

    my $tB45file   = "$plot.tB45";
    my $tB62file   = "$plot.tB62";
    my $tB90file   = "$plot.tB90";
    my $tPAM30file = "$plot.tP30";
    my $tPAM10file = "$plot.tP10";

    my $sB45file   = "$plot.sB45";
    my $sB62file   = "$plot.sB62";
    my $sB90file   = "$plot.sB90";
    my $sPAM30file = "$plot.sP30";
    my $sPAM10file = "$plot.sP10";

    open(T1, ">$WPfile") || die;
    printf T1 "91.1  8.9   -2  -16\n"; #VTML10
    printf T1 "83.4  16.6  -2  -15\n"; #VTML20
    printf T1 "68.9  31.1  -1  -13\n"; #VTML40
    printf T1 "48.7  51.3  -1  -10\n"; #VTML80
    printf T1 "35.5  64.5  -1  -11\n"; #VTML120
    printf T1 "28.4  71.6  -1  -10\n"; #VTML140
    printf T1 "28.0  72.9  %f  %f\n", -1.*2./3., -12.*2./3.; #VTML160
    printf T1 "29.9  70.1  -1  -11\n";                       #BLOSUM62
    printf T1 "27.1  72.9  %f  %f\n", -10.*2./3., -2.*2./3.; #BLOSUM50 
    close(T1);

     open(T1, ">$tB45file") || die;
    print T1 "$tB45\t 0\n";
    print T1 "$tB45\t $ymin\n";
    close(T1);
    open(T1, ">$tB62file") || die;
    print T1 "$tB62\t 0\n";
    print T1 "$tB62\t $ymin\n";
    close(T1);   
    open(T1, ">$tB90file") || die;
    print T1 "$tB90\t 0\n";
    print T1 "$tB90\t $ymin\n";
    close(T1);
    open(T1, ">$tPAM30file") || die;
    print T1 "$tP30\t 0\n";
    print T1 "$tP30\t $ymin\n";
    close(T1);   
    open(T1, ">$tPAM10file") || die;
    print T1 "$tP10\t 0\n";
    print T1 "$tP10\t $ymin\n";
    close(T1);   
    
    open(T1, ">$sB45file") || die;
    print T1 "$sB45\t 0\n";
    print T1 "$sB45\t $ymin\n";
    close(T1);
    open(T1, ">$sB62file") || die;
    print T1 "$sB62\t 0\n";
    print T1 "$sB62\t $ymin\n";
    close(T1);   
    open(T1, ">$sB90file") || die;
    print T1 "$sB90\t 0\n";
    print T1 "$sB90\t $ymin\n";
    close(T1);
    open(T1, ">$sPAM30file") || die;
    print T1 "$sP30\t 0\n";
    print T1 "$sP30\t $ymin\n";
    close(T1);       
    open(T1, ">$sPAM10file") || die;
    print T1 "$sP10\t 0\n";
    print T1 "$sP10\t $ymin\n";
    close(T1);   
    
 
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    print GP "set key right top\n";
    print GP "set yrange [$ymin:0.0]\n";
    
    my $cmd;
    my $field;
    my $tag;
    
    #print GP "set logscale y 10\n";

    print GP "set output '$tfile'\n";
    $field = 1;
    $xlabel = "divergence time t";
    print GP "set xlabel '$xlabel'\n";
    print GP "set xrange [0:$xmax]\n";
    $tag = "";
    $cmd  = "'$tB62file'   title '' with lines ls 3, ";
    $cmd .= "'$tB90file'   title '' with lines ls 3, ";
    $cmd .= "'$tPAM30file' title '' with lines ls 3, ";
    $cmd .= "'$tPAM10file' title '' with lines ls 3, ";
    $cmd .= "'$tB45file'   title '' with lines ls 3, ";
    $cmd .= "'$file2' using $field:4  title '' with lines ls 1115, "; 
    $cmd .= "'$file2' using $field:5  title '' with lines ls 1112, "; 
    $cmd .= "'$file2' using $field:6  title ''    with lines ls 1111, "; 
    #$cmd .= "'$file2' using $field:7  title 'go1' with lines ls 3, "; 
    #$cmd .= "'$file2' using $field:8  title 'go2' with lines ls 4, "; 
    #$cmd .= "'$file2' using $field:9  title 'go3' with lines ls 5, "; 
    #$cmd .= "'$file2' using $field:10 title 'go4' with lines ls 2, "; 
    #$cmd .= "'$file2' using $field:11 title 'go5' with lines ls 7, "; 
    #$cmd .= "'$file2' using $field:12 title 'go6' with lines ls 8, "; 
    $cmd .= "'$file1' using $field:4  title '' ls 1115, "; 
    $cmd .= "'$file1' using $field:5  title '' ls 1112, "; 
    $cmd .= "'$file1' using $field:6  title '' ls 1111"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";

    print GP "set output '$sfile'\n";
    $field = 3;    
    $xlabel = "% substitutions per site";
    print GP "set xlabel '$xlabel'\n";
    print GP "set xrange [0:95]\n";
    $tag = "";
    $cmd  = "'$sB62file'                  title '' with lines ls 3, ";
    $cmd .= "'$sB90file'                  title '' with lines ls 3, ";
    $cmd .= "'$sPAM30file'                title '' with lines ls 3, ";
    $cmd .= "'$sPAM10file'                title '' with lines ls 3, ";
    $cmd .= "'$sB45file'                  title '' with lines ls 3, ";
    $cmd .= "'$file1'       using $field:5 title '' ls 1112, "; 
    $cmd .= "'$file1'       using $field:6 title '' ls 1111, "; 
    $cmd .= "'$file2'       using $field:5 title '' with lines ls 1112, "; 
    $cmd .= "'$file2'       using $field:6 title '' with lines ls 1111, "; 
    $cmd .= "'$WPfile'     using 2:3      title '' ls 88, ";
    $cmd .= "'$WPfile'     using 2:4      title '' ls 88";
    $cmd .= "\n";
    print GP "plot $cmd\n";

    print GP "set output '$ifile'\n";
    $field = 2;    
    $xlabel = "% identities per conserved site";
    print GP "set xlabel '$xlabel'\n";
    print GP "set xrange [100:5]\n";
    $tag = "";
    $cmd   = "'$sB62file'                 title '' with lines ls 3, ";
    $cmd  .= "'$sB90file'                 title '' with lines ls 3, ";
    $cmd  .= "'$sPAM30file'               title '' with lines ls 3, ";
    $cmd  .= "'$sPAM10file'               title '' with lines ls 3, ";
    $cmd  .= "'$sB45file'                 title '' with lines ls 3, ";
    $cmd .= "'$file1'       using $field:5 title '' ls 1112, "; 
    $cmd .= "'$file1'       using $field:6 title '' ls 1111, "; 
    $cmd .= "'$file2'       using $field:5 title '' with lines ls 1112, "; 
    $cmd .= "'$file2'       using $field:6 title '' with lines ls 1111, "; 
    $cmd .= "'$WPfile'     using 1:3      title '' ls 88, ";
    $cmd .= "'$WPfile'     using 1:4      title '' ls 88";
    $cmd .= "\n";
    print GP "plot $cmd\n";

    close(GP);

    system("ps2pdf $tfile\n");
    system("ps2pdf $sfile\n");
    system("ps2pdf $ifile\n");
    
    system("rm $tfile\n");
    system("rm $sfile\n");
    system("rm $ifile\n");
    
    if ($viewplot) {
	system("open $pdf_tfile&\n"); 
	system("open $pdf_sfile&\n"); 
	#system("open $pdf_ifile&\n"); 
    }
}
