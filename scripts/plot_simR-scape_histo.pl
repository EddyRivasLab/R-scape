#!/usr/bin/perl -w
#plot_simR-scape_histo.pl

use strict;
use Class::Struct;
use lib '/Users/rivase/projects/R-scape/scripts';
use FUNCS;
#use constant GNUPLOT => '/usr/bin/gnuplot';
use constant GNUPLOT => '/opt/local/bin/gnuplot';

use vars qw ($opt_x $opt_X $opt_v $opt_V);  # required if strict used
use Getopt::Std;
getopts ('x:X:vV');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  plot_simR-scape_histo.pl [options] <outpsfile> <N> <histo_1>..<histo_N> \n\n";
        print "options:\n";
 	print "-x <s>    : xlabel\n";
 	print "-X <x>    : xmax\n";
 	print "-v        : be verbose\n";
 	print "-V        : viewplots\n";
 	exit;
}
my $psfile = shift;
my $N = shift;
my @hisfile;
for (my $n = 0; $n < $N; $n ++) {
    $hisfile[$n]  = shift;
}

my $xfield = 1;
my $yfield = 2;
my $title = "$psfile";
my $xlabel = "max cov score"; if ($opt_x) { $xlabel = "$opt_x"; }
my $ylabel = "ocurrence";
my $iscum = 0;
my $seeplots = 1;
my $xleft = "";
my $xright = ""; if ($opt_X) { $xright = "$opt_X"; }
my $ymax = 20;
gnuplot_N_histo($N, \@hisfile, $xfield, $yfield, $psfile, $title, $xlabel, $ylabel, $iscum, $seeplots, $xleft, $xright, $ymax);


sub gnuplot_N_histo {

    my ($N, $hfile_ref, $xfield, $yfield, $psfile, $title, $xlabel, $ylabel, $iscum, $seeplots, $xleft, $xright, $ymax) = @_;

    my $pdffile = $psfile;
    if ($pdffile =~ /^(\S+).ps$/) { $pdffile = "$1.pdf"; }
    
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    
    print GP "set terminal postscript color solid 14\n";
    FUNCS::gnuplot_define_styles (*GP);

    print GP "set output '$psfile'\n";
    #print GP "set nokey\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set xrange [$xleft:$xright]\n";
    print GP "set yrange [0:$ymax]\n";
    print GP "set style histogram cluster gap 1\n";
    print GP "set style fill solid noborder \n";

    #print GP "set title \"$title\\n\\n$key\"\n";
    print GP "set title '$title'\n";

    print GP "set ylabel '$ylabel'\n";

    my $m = 2;
    my $cmd = "";
    for (my $n = 0; $n < $N; $n ++) {
	my $mean_avgid;
	my $stdv_avgid;
	parse_hisfile($hfile_ref->[$n], \$mean_avgid, \$stdv_avgid);
	my $key = "$hfile_ref->[$n]";
	if ($key =~ /([^\/]+)$/) { $key = "$1\_$mean_avgid"; }

	if ($iscum) {
	    $cmd .= ($n == $N-1)? 
		"'$hfile_ref->[$n]' using $xfield:$yfield  with lines title '$key' ls $m" : 
		"'$hfile_ref->[$n]' using $xfield:$yfield  with lines title '$key' ls $m, ";
	}
	else {
	    $cmd .= ($n == $N-1)? 
		"'$hfile_ref->[$n]' using $xfield:$yfield  with boxes title '$key' ls $m" : 
		"'$hfile_ref->[$n]' using $xfield:$yfield  with boxes title '$key' ls $m, ";
	} 
	$m ++;
    }

    print  "plot $cmd\n";
    print GP "plot $cmd\n";

    close (GP);

    system ("ps2pdf $psfile\n"); 
    system("rm $psfile\n");
    if ($seeplots) { system ("open $pdffile&\n"); }
}

sub parse_hisfile {
    my ($hisfile, $ret_mean_avgid, $ret_stdv_avgid)= @_;

    my $mean_avgid = -1;
    my $stdv_avgid = -1;

    open (HIS, "$hisfile") || die; 
    while (<HIS>) {
	if    (/\# mean_avgid (\S+) stdv_avgid (\S+)/)  { $mean_avgid = $1; $stdv_avgid = $2; }
    }
    close(HIS);
    
    $$ret_mean_avgid = int($mean_avgid * 1000)/1000;
    $$ret_stdv_avgid = int($stdv_avgid * 1000)/1000;
}
