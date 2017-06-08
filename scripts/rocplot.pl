#!/usr/bin/perl -w
# rocplot.pl

use strict;
use Class::Struct;

use vars qw ($opt_D $opt_L $opt_P $opt_R $opt_W $opt_v);  # required if strict used
use Getopt::Std;
getopts ('D:L:P:RW:v');


# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  rocplot.pl [options] <F> <file1>..<fileF> <stofile> <rscapebin> <gnuplotdir> \n\n";
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
my $stofile  = shift;
my $pfamname = $prefile[0];
if ($pfamname =~ /(PF[^\.]+)\./) { $pfamname = $1; }
if ($pfamname =~ /(RF[^\.]+)\./) { $pfamname = $1; }
my $oafafile = "data/$pfamname.afa";
print "original file: $oafafile\n";

my $rscapebin = shift;
my $gnuplot   = shift;
#use constant GNUPLOT => '$gnuplot';
use constant GNUPLOT => '/usr/local/bin/gnuplot';

my $lib;
#BEGIN { $lib = "$rscapebin/../scripts" };
#use lib '$lib';
use lib '/Users/rivase/src/src/mysource/scripts';
use PDBFUNCS;
use FUNCS;

my $currdir = $ENV{PWD};

my $maxD = 8;
if ($opt_D) { $maxD = $opt_D; }

my $minL = -1;
if ($opt_L) { $minL = $opt_L; }
for (my $f = 0; $f < $F; $f ++) {  $rocfile[$f] = ($minL>0)? "$prename[$f].minL$minL.roc":"$prename[$f].roc"; }

my $dornaview = 0;
if ($opt_R) { $dornaview = 1; }
my $which = "MIN"; #options: CA C MIN AVG NOH / C1' (for RNA suggested by Westhof)
if ($opt_W) { $which = "$opt_W"; }
my $seeplots = 0;
my $ncnt = 0;
my $nbp  = 0;
my $nwc  = 0;
my @cnt;
my $pdbfile = "";
if ($opt_P) { 
    $pdbfile = "$opt_P";
    PDBFUNCS::contacts_from_pdbfile ($gnuplot, $rscapebin, $pdbfile, $stofile, \$ncnt, \@cnt, $maxD, $minL, $which, $dornaview, "", "", $seeplots);
    PDBFUNCS::print_contactlist(\*STDOUT, $ncnt, \@cnt);
    PDBFUNCS::bpinfo_contactlist($ncnt, \@cnt, \$nbp, \$nwc);
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

my $alen = -1;
my $alenDCA = -1;
my @mapDCA;
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
    elsif ($method =~ /^mfDCA$/) {
	$alenDCA = create_rocfile_mfDCA ($rocfile[$f], $prefile[$f], $stofile, $minL, $ncnt, $nbp, $nwc, \@cnt,  \$alen, \$alenDCA, \@mapDCA);
    }
    elsif ($method =~ /^plmDCA$/) {
	$alenDCA = create_rocfile_plmDCA($rocfile[$f], $prefile[$f], $stofile, $minL, $ncnt, $nbp, $nwc, \@cnt, \$alen, \$alenDCA, \@mapDCA);
    }
    else { print "method $method not implemented yet\n"; die; }
    
   
}

my $xmax = 100;
my $viewplots = 1;
rocplot($pfamname, $F, \@rocfile, \@prename, $xmax, $viewplots);






####################### routines

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


sub  create_rocfile_mfDCA {
    my ($rocfile, $prefile, $stofile, $minL, $ncnt, $nbp, $nwc, $cnt_ref, $ret_alen, $ret_alenDCA, $mapDCA_ref, $which) = @_;

    my $method = "mfDCA";
    my $which  = "DI";
    if ($rocfile =~ /mfDCA\.MI/) {
	$which = "MI";
	$method .= "-$which";
    }

    my $alenDCA = $$ret_alenDCA;
    if ($alenDCA < 0) {
	mapDCA2MSA($stofile, $mapDCA_ref, $ret_alen, \$alenDCA);
    }

    parse_mfDCA($rocfile, $prefile, $minL, $ncnt, $nbp, $nwc, $cnt_ref, $mapDCA_ref, $$ret_alen, $alenDCA, $which);
    
    $$ret_alenDCA = $alenDCA;
    
}

sub  create_rocfile_plmDCA {
    my ($rocfile, $prefile, $stofile, $minL, $ncnt, $nbp, $nwc, $cnt_ref, $ret_alen, $ret_alenDCA, $mapDCA_ref) = @_;

    my $method = "plmDCA";
  
    my $alenDCA = $$ret_alenDCA;
    if ($alenDCA < 0) {
	mapDCA2MSA($stofile, $mapDCA_ref, $ret_alen, \$alenDCA);
    }

    parse_plmDCA($rocfile, $prefile, $minL, $ncnt, $nbp, $nwc, $cnt_ref, $mapDCA_ref, $$ret_alen, $alenDCA);

    $$ret_alenDCA = $alenDCA;
  
}

sub  create_rocfile_random {
    my ($rocfile, $prefile) = @_;
}



# map the coordenates in DCA output files to those of the input alignment
#
# mlDCA and plmDCA trim internally the alignmen by removing coluns with . or lower case
# we need mapDCA() because mlDCA and plmDCA report coordinates relative
# to this filtered alignment, not relative to the input alignment
#
# alen    length of the input alignment
# alenDCA length of the columns used by DCA (no . no lower case)
# mapDCA[1..alenDCA] valued in 1..alen
sub mapDCA2MSA {
    my ($stofile, $mapDCA_ref, $ret_alen, $ret_alenDCA) = @_;
    my $alen     = 0;
    my $asq      = "";
    my $alenDCA  = 0;
    my $asqDCA   = "";

    my $afafile = "$stofile.afa";
    my $reformat = "$rscapebin/../lib/hmmer/easel/miniapps/esl-reformat afa ";
    system("$reformat $stofile > $afafile\n");
    
    # grab the first sequence
    my $n = 0;
    open (FA, "$afafile") || die;
    while (<FA>) {
	if (/^\>/) { $n ++; if ($n > 1) { last; } }
	elsif ($n == 1 && /^(\S+)\s*$/) {
	    $asq .= $1;
	}       
    }
    close (FA);
    
    $alen = length($asq);
    
    my $aux = $asq;
    my $char;
    my $apos = 0;
    my $pos = 0;
    while ($aux) {
	$apos ++;
	$aux =~ s/^(\S)//; $char = $1;
	if ($char =~ /^\.$/ || $char =~ /^[a-z]$/) {
	}
	else { 
	    $asqDCA .= "$char"; 
	    $pos ++;     
	    $mapDCA_ref->[$pos] = $apos;
	}
    }
    $alenDCA = length($asqDCA);
    
    print "alen    $alen\n$asq\n";
    print "alenDCA $alenDCA\n$asqDCA\n";

    system("rm $afafile\n");
    $$ret_alenDCA = $alenDCA;
    $$ret_alen    = $alen;
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


sub parse_mfDCA {
    my ($rocfile, $file, $minL, $ncnt, $nbp, $nwc, $cnt_ref, $mapDCA_ref, $alen, $alenDCA, $which) = @_;
    my $npre = 0;

    my $sortfile = sort_mfDCA($file, $which);
    
    open my $fp, '>', $rocfile || die "Can't open $rocfile: $!";

    my $f   = 0;
    my $f_c = 0;
    my $f_b = 0;
    my $f_w = 0;
    my $t_c = $ncnt;
    my $t_b = $nbp;
    my $t_w = $nwc;

    my $type;
    open(FILE, "$sortfile") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+\S+\s*$/) {
	    my $idca = $1;
	    my $jdca = $2;
	    my $i    = $mapDCA_ref->[$idca];
	    my $j    = $mapDCA_ref->[$jdca];
	    if (PDBFUNCS::found_in_contactlist($i, $j, $minL, $ncnt, $cnt_ref, \$type)) {
		if    ($type ==  0) { $f_w ++; }
		elsif ($type <  12) { $f_b ++; }
		$f_c ++;
	    }
	    $f ++;
	    writeline($fp,      $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $alen);
	    writeline(\*STDOUT, $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $alen);
	}
    }
    close(FILE);
    close($fp);

    system("rm $sortfile\n");
}

sub parse_plmDCA {
    my ($rocfile, $file, $minL, $ncnt, $nbp, $nwc, $cnt_ref, $mapDCA_ref, $alen, $alenDCA) = @_;
    my $npre = 0;

    my $sortfile = sort_plmDCA($file);

    open(FILE, "$sortfile") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+\S+\s*$/) {
	    my $i = $1;
	    my $j = $2;

	}
    }
    close(FILE);

   system("rm $sortfile\n");
}

sub parse_gremplin {
    my ($file, $minL) = @_;
    my $npre = 0;

    my $sortfile = sort_gremlin($file);

    open(FILE, "$sortfile") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+\S+\s*$/) {
	    my $i = $1;
	    my $j = $2;
	}
    }
    close(FILE);

    system("rm $sortfile\n");
}

sub sort_mfDCA {
    my ($file, $which) = @_;

    my $sortfile = "$file.sort";
    my $n = 0;
    my %i;
    my %j;
    my %sc;
    my $di;
    my $mi;
    open(SORT, ">$sortfile") || die;
    open(FILE, "$file") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s*$/) {
	    $i{$n} = $1;
	    $j{$n} = $2;
	    $mi    = $3;
	    $di    = $4;
	    if    ($which =~ /^DI$/) { $sc{$n} = $di; }
	    elsif ($which =~ /^MI$/) { $sc{$n} = $mi; }
	    $n ++;
	}
    }
    close(FILE);
    
    my @newkey = sort { $sc{$b} <=> $sc{$a} } keys(%sc);
    my @newi   = @i{@newkey};
    my @newj   = @j{@newkey};
    my @newsc  = @sc{@newkey};

   for (my $x = 0; $x < $n; $x ++) {
	print SORT "$newi[$x] $newj[$x] $newsc[$x]\n";
    }
    close(SORT);
   
    return $sortfile;
}

sub sort_plmDCA {
    my ($file) = @_;

    my $sortfile = "$file.sort";

    my $n = 0;
    my %i;
    my %j;
    my %sc;
    my $di;
    my $mi;
    my %newi;
    my %newj;
    my %newsc;
    open(SORT, ">$sortfile") || die;
    open(FILE, "$file") || print "FILE DOES NOT EXIST\n";
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+),(\d+),(\S+)\s*$/) {
	    $i{$n}  = $1;
	    $j{$n}  = $2;
	    $sc{$n} = $3;
	    $n ++;
	}
    }
    close(FILE);

    my @newkey = sort { $sc{$b} <=> $sc{$a} } keys(%sc);
    my @newi   = @i{@newkey};
    my @newj   = @j{@newkey};
    my @newsc  = @sc{@newkey};

    for(my $x = 0; $x < $n; $x ++) {
	print SORT "$newi[$x] $newj[$x] $newsc[$x]\n";
    }
    close(SORT);
    
    return $sortfile;
}

sub sort_gremlin {
    my ($file) = @_;

    my $sortfile = "$file.sort";

    my $n = 0;
    my %i;
    my %j;
    my %sc;
    my $di;
    my $mi;
    my %newi;
    my %newj;
    my %newsc;
    
    my $row = 0;
    my $dim;
    open(SORT, ">$sortfile") || die;
    print "file:$file\n";
    open(FILE, "$file") || print "FILE DOES NOT EXIST\n";
    while(<FILE>) {
	my $line = $_;
	$line =~ s/\n//g;
	my @vec = split(/ /, $line);

	$dim = $#vec+1;
	for (my $col = $row+1; $col < $dim; $col ++) {
	    $i{$n}  = $row+1;
	    $j{$n}  = $col+1;
	    $sc{$n} = $vec[$col];
	    $n ++;
	}
	$row ++;
    }
    close(FILE);
     if ($row != $dim) { print "sort_gremlin() error\n"; }
    
    my @newkey = sort { $sc{$b} <=> $sc{$a} } keys(%sc);
    my @newi   = @i{@newkey};
    my @newj   = @j{@newkey};
    my @newsc  = @sc{@newkey};

    for(my $x = 0; $x < $n; $x ++) {
	print SORT "$newi[$x] $newj[$x] $newsc[$x]\n";
    }
    close(SORT);
    
    return $sortfile;
}
