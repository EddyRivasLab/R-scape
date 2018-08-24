#!/usr/bin/perl -w
# 
#  Rank a R-scape output by the F measure
#

use strict;
use Class::Struct;

use constant GNUPLOT => '/usr/local/bin/gnuplot';
use lib '/Users/erivas/src/src/mysource/scripts';
use FUNCS;

use vars qw ($opt_v $opt_F);  # required if strict used
use Getopt::Std;
getopts ('vF:');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  rank_by_F.pl [options] <R-scape.out> \n\n";
        print "options:\n";
 	print "-v    :  be verbose\n";
 	exit;
}
my $file    = shift;
my $outfile = "$file";
if ($outfile =~ /^(\S+).rscape/) { $outfile = "$1.rank"; }
else                             { $outfile .= ".rank";  }

my $file3d;
if ($opt_F) { $file3d = "$opt_F"; }

my $n3d = 0;
my @fam3d;
my @name3d;
if ($file3d) { 
    parse_3dfile($file3d, \$n3d, \@fam3d, \@name3d);
    print "3D fam $n3d\n";
    for (my $f = 0; $f < $n3d; $f ++) {
	printf "%d |$fam3d[$f]| |$name3d[$f]|\n", $f+1;
    }
}

# Method Target_E-val [cov_min,conv_max] [FP | TP True Found | Sen PPV F] 
# GTp    0.05         [-9.85,91.47]     [0 | 5 6 5 | 83.33 100.00 90.91] 

my @fam;
my %fam_id;
my %fam_alen;
my %fam_nseq;
my %fam_avgsub;

my %fam_tp;
my %fam_true;
my %fam_found;
my %fam_bpexp;
my %fam_power;
my %fam_S;
my %fam_P;
my %fam_F;

my %fam_tp_cyk;
my %fam_true_cyk;
my %fam_found_cyk;
my %fam_bpexp_cyk;
my %fam_power_cyk;
my %fam_S_cyk;
my %fam_P_cyk;
my %fam_F_cyk;
my $nf = 0;
my $fam;
my $usefam;
my $acc;
my $iscyk;

open (OUT, ">$outfile") || die;
open (FILE, "$file")    || die;

while(<FILE>) {
    # MSA RF00001_5S_rRNA nseq 712 (712) alen 119 (230) avgid 56.09 (55.86) nbpairs 34 (34)
    if (/^# MSA\s+(\S+)\s+nseq\s+(\S+)\s+\(.+alen\s+(\S+)\s+\(.+avgid\s+(\S+)\s+/) {
	$fam      = $1;
	my $nseq  = $2;
	my $alen  = $3;
	my $avgid = $4;
	$acc = $fam;
	if ($acc =~ /^(RF\d\d\d\d\d)/) { $acc = $1; }
	
	$usefam = 0;
	$iscyk  = 0;
	if ($n3d == 0 || ($n3d > 0 && is_3dfam($acc, $n3d, \@fam3d)) ) { $usefam = 1; }	    

 
	if ($usefam) {
	    $fam[$nf]       = $fam;
	    $fam_nseq{$fam} = $nseq;
	    $fam_alen{$fam} = $alen;
	    $fam_id{$fam}   = $avgid;
	    
	    $fam_tp{$fam}    = 0.0;
	    $fam_true{$fam}  = 0.0;
	    $fam_found{$fam} = 0.0;
	    
	    $fam_S{$fam}     = 0.0;
	    $fam_P{$fam}     = 0.0;
	    $fam_F{$fam}     = 0.0;
	    $fam_power{$fam} = 0.0;
	    
	    $fam_tp_cyk{$fam}    = 0.0;
	    $fam_true_cyk{$fam}  = 0.0;
	    $fam_found_cyk{$fam} = 0.0;
	    
	    $fam_S_cyk{$fam}     = 0.0;
	    $fam_P_cyk{$fam}     = 0.0;
	    $fam_F_cyk{$fam}     = 0.0;
	    $fam_power_cyk{$fam} = 0.0;
	}
    }
    # Method Target_E-val [cov_min,conv_max] [FP | TP True Found | Sen PPV F] 
    # GTp    0.05         [-9.35,1389.53]     [8 | 22 34 30 | 64.71 73.33 68.75] 
    elsif (/^#\s+\S+\s+\S+\s+\[\S+\]\s+\[\d+\s+\|\s+(\d+)\s+(\d+)\s+(\d+)\s+\|\s+(\S+)\s+(\S+)\s+(\S+)\]/) {
	my $tp      = $1;
	my $true    = $2;
	my $found   = $3;
	
	my $sen     = $4;
	my $ppv     = $5;
	my $F       = $6;
	
	if ($usefam) {	    
	    #printf "%d FAM $fam[$nf] sen $sen ppv $ppv F $F\n", $nf+1;
	    if ($iscyk) {
		$fam_tp_cyk{$fam}    = $tp;
		$fam_true_cyk{$fam}  = $true;
		$fam_found_cyk{$fam} = $found;
		
		$fam_S_cyk{$fam}     = $sen;
		$fam_P_cyk{$fam}     = $ppv;
		$fam_F_cyk{$fam}     = $F;
		$fam_power_cyk{$fam} = ($true > 0)? int(10000*$fam_bpexp{$fam}/$true)/100 : 0.0;
	    }
	    else {
		$fam_tp{$fam}    = $tp;
		$fam_true{$fam}  = $true;
		$fam_found{$fam} = $found;
		
		$fam_S{$fam}     = $sen;
		$fam_P{$fam}     = $ppv;
		$fam_F{$fam}     = $F;
		$fam_power{$fam} = ($true > 0)? int(10000*$fam_bpexp{$fam}/$true)/100 : 0.0;
		$nf ++;
	    }
	}
    }
    # BPAIRS 20
    # avg substitutions per BP 15.8
    # BPAIRS expected covary 7.1
    elsif (/^#\s+avg substitutions per BP\s+(\S+)/) {
	if ($usefam) {
	    my $avgsub = $1;
	    $fam_avgsub{$fam} = $avgsub;
	}
    }
    elsif (/^#\s+BPAIRS expected covary\s+(\S+)/) {
	if ($usefam) {
	    my $bpexp = $1;
	    $fam_bpexp{$fam} = $bpexp;
	}
    }
    elsif (/^# The predicted cyk-cov structure/) {
	$iscyk = 1;
    }
    
}
close (FILE);

create_rank_plot($outfile);

my $m = 0;
for (my $f = 0; $f < $nf; $f ++) {
    my $fam = $fam[$f];
    my $sen = $fam_S{$fam};
    my $power = $fam_power{$fam};
    if ($sen < 2 && $power > 10) {
	$m ++;
	print "outlier $m $fam sen $sen power $power\n";
    }
}
$m = 0;
for (my $f = 0; $f < $nf; $f ++) {
    my $fam = $fam[$f];
    my $sen = $fam_S{$fam};
    my $power = $fam_power{$fam};
    if ($power < 2) {
	$m ++;
	print "no power $m $fam sen $sen power $power\n";
    }
}

if ($iscyk) {
    my $n = 0;
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam     = $fam[$f];
	my $sen     = $fam_S{$fam};
	my $sen_cyk = $fam_S_cyk{$fam};
	if ($sen_cyk > $sen) {
	    $n ++;
	    print "better ss $n $fam sen $sen sen_cyk $sen_cyk\n";
	}
    }
    
    $n = 0;
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam           = $fam[$f];
	my $tp            = $fam_tp{$fam};
	my $tp_cyk        = $fam_tp_cyk{$fam};
	my $sen           = $fam_S{$fam};
	my $sen_cyk       = $fam_S_cyk{$fam};
	if ($sen_cyk > $sen  && $tp_cyk > $tp+2) {
	    $n ++;
	    print "better ss2 $n $fam sen $sen sen_cyk $sen_cyk | tp $tp tp_cyk $tp_cyk\n";
	}
    }
    $n = 0;
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam           = $fam[$f];
	my $tp            = $fam_tp{$fam};
	my $tp_cyk        = $fam_tp_cyk{$fam};
	my $sen           = $fam_S{$fam};
	my $sen_cyk       = $fam_S_cyk{$fam};
	my $F             = $fam_F{$fam};
	my $F_cyk         = $fam_F_cyk{$fam};
	if ($tp_cyk < $tp ) {
	    $n ++;
	    print "worse ss $n $fam sen $sen sen_cyk $sen_cyk | tp $tp tp_cyk $tp_cyk\n";
	}
    }
}


sub create_rankF_extra {
    my ($outfile) = @_;
    
    my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} } keys(%fam_S);
    my @P_order = sort { $fam_P{$b} <=> $fam_P{$a} } keys(%fam_P);
    my @F_order = sort { $fam_F{$b} <=> $fam_F{$a} } keys(%fam_F);

    my $n = 0;
    printf OUT "\n# nfam = $nf\n";
    printf OUT "# family\t\t\t\t\t\tF\tSEN\tPPV\tPOWER\tTRUE\tFOUND\tTP\tEXP\tavgsub\tavgid\talen\tnseq\n";
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam = $F_order[$f];
	if ($fam_true{$fam} > 0 ) {
	    $n ++;
	    printf OUT "%d\t%-40s\t$fam_F{$fam}\t$fam_S{$fam}\t$fam_P{$fam}\t$fam_power{$fam}\t$fam_true{$fam}\t$fam_found{$fam}\t$fam_tp{$fam}\t$fam_F_cyk{$fam}\t$fam_S_cyk{$fam}\t$fam_P_cyk{$fam}\t$fam_power_cyk{$fam}\t$fam_true_cyk{$fam}\t$fam_found_cyk{$fam}\t$fam_tp_cyk{$fam}\t$fam_bpexp{$fam}\t$fam_avgsub{$fam}\t$fam_id{$fam}\t$fam_alen{$fam}\t$fam_nseq{$fam}\n", $n, $fam;
	}
    }
    close(OUT);

    plot_rank($outfile);
    
    if ($iscyk) { plot($outfile, 6, 4, 13, 11); }
    else        { plot($outfile, 6, 4, -1, -1); }
}

sub create_rank_plot {
    my ($outfile) = @_;
    
    my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} } keys(%fam_S);
    my @P_order = sort { $fam_P{$b} <=> $fam_P{$a} } keys(%fam_P);
    my @F_order = sort { $fam_F{$b} <=> $fam_F{$a} } keys(%fam_F);

    my $n = 0;
    printf OUT "\n# nfam = $nf\n";
    printf OUT "# family\t\t\t\t\t\tF\tSEN\tPPV\tPOWER\tTRUE\tFOUND\tTP\tEXP\tavgsub\tavgid\talen\tnseq\n";
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam = $S_order[$f];
	if ($fam_true{$fam} >= 0 ) {
	    $n ++;
	    printf OUT "%d\t%-40s\t$fam_F{$fam}\t$fam_S{$fam}\t$fam_P{$fam}\t$fam_power{$fam}\t$fam_true{$fam}\t$fam_found{$fam}\t$fam_tp{$fam}\t$fam_bpexp{$fam}\t$fam_avgsub{$fam}\t$fam_id{$fam}\t$fam_alen{$fam}\t$fam_nseq{$fam}\n", $n, $fam;
	}
    }
    close(OUT);

    plot_rank($outfile);
    
    plot($outfile, 6, 4, -1, -1); 
}


sub parse_3dfile {
    my ($file3d, $ret_n3d, $fam3d_ref, $name3d_ref) = @_;

    my $n3d = 0;

    open(FILE, "$file3d");
    while (<FILE>) {
	if (/^(.+)(RF\d+)\s*\n/) {
	    $name3d_ref->[$n3d] = $1;
	    $fam3d_ref->[$n3d]  = $2;
	    $n3d ++;
	}
    }
    close(FILE);
    
    $$ret_n3d = $n3d;
}

sub is_3dfam {
    my ($acc, $n3d, $fam3d_ref) = @_;

    for (my $f = 0; $f < $n3d; $f ++) {
	if ($acc =~ /^$fam3d_ref->[$f]$/) { return 1; }
    }

    return 0;
}

sub plot {
    my ($file, $xfield, $yfield, $xfield2, $yfield2) = @_;

    my $pdf = "$file.ps";
    my $xlabel = "Expected Covaring  (% basepairs)";
    my $ylabel = "Significantly Covaring (% basepairs)";
    
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    print GP "set output '$pdf'\n";
    print GP "set key right top\n";
    print GP "set nokey\n";
    print GP "set size square\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    print GP "set xrange [-0.2:10]\n";
    print GP "set yrange [-0.2:10]\n";
    #print GP "set logscale xy\n";

    my $cmd;
    $cmd .= "'$file' using $xfield:$yfield  title '' with points ls 2, "; 
    if ($xfield2 > 0) { $cmd .= "'$file' using $xfield2:$yfield2  title '' with points ls 8"; }
    $cmd .= "\n";
    print GP "plot $cmd\n";
    
    print GP "set xrange [-1:80]\n";
    print GP "set yrange [-0.2:10]\n";
    print GP "set size nosquare\n";

    $cmd = "";
    $cmd .= "'$file' using $xfield:$yfield  title '' with points ls 2, "; 
    if ($xfield2 > 0) { $cmd .= "'$file' using $xfield2:$yfield2  title '' with points ls 8"; }
    $cmd .= "\n";
    print GP "plot $cmd\n";
    
    print GP "set xrange [-1:100]\n";
    print GP "set yrange [-1:100]\n";
    print GP "set size square\n";

    $cmd = "";
    $cmd .= "'$file' using $xfield:$yfield  title '' with points ls 2, "; 
    if ($xfield2 > 0) { $cmd .= "'$file' using $xfield2:$yfield2  title '' with points ls 8"; }
    $cmd .= "\n";
    print GP "plot $cmd\n";
    
    
    system("open $pdf\n");

}

sub plot_rank {
    my ($file) = @_;

    my $pdf = "$file.rank.ps";
    my $xlabel = "RNA familiy";
    my $ylabel;
    my $cmd;
    
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    print GP "set output '$pdf'\n";
    print GP "set key right top\n";
    print GP "set nokey\n";
    print GP "set xlabel '$xlabel'\n";

    $ylabel = "Sensitivity (%)";
    print GP "set ylabel '$ylabel'\n";
    $cmd  = "";
    $cmd .= "'$file' using 1:4   title '' with boxes ls 1"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
    
    $ylabel = "F measure (%)";
    print GP "set ylabel '$ylabel'\n";
    $cmd .= "'$file' using 1:3   title '' with boxes ls 1"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";

    
    system("open $pdf\n");

}

