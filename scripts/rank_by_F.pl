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
my $outname = "$file";
if ($outname =~ /^(\S+).rscape/) { $outname = "$1"; }

my $outfile_rank         = "$outname.rank";        # file with all families ranked by Sensitivity/power 
my $outfile_allfam       = "$outname.allfam";      # file with all families ranked by Sensitivity to plot
my $outfile_nopower      = "$outname.nopower";     # families with no power (power = 0) 
my $outfile_withpower    = "$outname.withpower";   # families with power   (power > 0)
my $outfile_nocov        = "$outname.nocov";       # families with no significant covariations   (sen = 0)
my $outfile_withcov      = "$outname.withcov";     # families with    significant covariations   (sen > 0)
my $outfile_greyzone     = "$outname.greyzone";
my $outfile_outlierp     = "$outname.outlierp";
my $outfile_outlierc     = "$outname.outlierc";
my $outfile_goodpower    = "$outname.goodpower";
my $outfile_betterss     = "$outname.betterss";
my $outfile_muchbettss   = "$outname.muchbetterss";
my $outfile_worsess      = "$outname.worsess";
my $outfile_equalss      = "$outname.equalss";

my $file3d;
my $n3d = 0;
my @fam3d;
my @name3d;
if ($opt_F) { 
    $file3d = "$opt_F"; 
    parse_3dfile($file3d, \$n3d, \@fam3d, \@name3d);
    print "3D fam $n3d\n";
    for (my $f = 0; $f < $n3d; $f ++) {
	printf "%d |$fam3d[$f]| |$name3d[$f]|\n", $f+1;
    }
}

# Method Target_E-val [cov_min,conv_max] [FP | TP True Found | Sen PPV F] 
# GTp    0.05         [-9.85,91.47]     [0 | 5 6 5 | 83.33 100.00 90.91] 

my @fam;
my %fam_idx;
my %fam_id;
my %fam_alen;
my %fam_nseq;

my %fam_tp;
my %fam_true;
my %fam_found;
my %fam_tpexp;
my %fam_avgsub;

my %fam_S;
my %fam_P;
my %fam_F;
my %fam_Spower;

my %fam_tp_cyk;
my %fam_true_cyk;
my %fam_found_cyk;
my %fam_tpexp_cyk;
my %fam_avgsub_cyk;

my %fam_S_cyk;
my %fam_P_cyk;
my %fam_F_cyk;
my %fam_Spower_cyk;

my %fam_all;

my $nf = 0;
my $fam;
my $usefam;
my $acc;
my $iscyk;

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
	    $fam_idx{$fam}  = $nf+1;
	    $fam_nseq{$fam} = $nseq;
	    $fam_alen{$fam} = $alen;
	    $fam_id{$fam}   = $avgid;
	    
	    $fam_tp{$fam}         = 0.0;
	    $fam_true{$fam}       = 0.0;
	    $fam_found{$fam}      = 0.0;
	    $fam_tpexp{$fam}      = 0.0;
	    
	    $fam_S{$fam}          = 0.0;
	    $fam_P{$fam}          = 0.0;
	    $fam_F{$fam}          = 0.0;
	    
	    $fam_Spower{$fam}     = 0.0;
	    
	    $fam_tp_cyk{$fam}     = 0.0;
	    $fam_true_cyk{$fam}   = 0.0;
	    $fam_found_cyk{$fam}  = 0.0;
	    $fam_tpexp_cyk{$fam}  = 0.0;
	    
	    $fam_S_cyk{$fam}      = 0.0;
	    $fam_P_cyk{$fam}      = 0.0;
	    $fam_F_cyk{$fam}      = 0.0;
	    
	    $fam_Spower_cyk{$fam} = 0.0;
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
	    }
	    else {
		$fam_tp{$fam}    = $tp;
		$fam_true{$fam}  = $true;
		$fam_found{$fam} = $found;
		$nf ++;
	    }
	}
    }
    # BPAIRS 20
    # avg substitutions per BP 15.8
    # BPAIRS expected covary 7.1
    elsif (/^#\s+avg substitutions per BP\s+(\S+)/) {
	my $avgsub = $1;
	if ($usefam) {
	    if ($iscyk) { $fam_avgsub_cyk{$fam} = $avgsub; }
	    else        { $fam_avgsub{$fam}     = $avgsub; }
	}
    }
    elsif (/^#\s+BPAIRS expected covary\s+(\S+)/) {
	my $tpexp = $1;
	if ($usefam) { 
	    if ($iscyk) { $fam_tpexp_cyk{$fam} = $tpexp; }
	    else        { $fam_tpexp{$fam}     = $tpexp; }
	}
    }
    elsif (/^# The predicted cyk-cov structure/) {
	$iscyk = 1;
    }
    # add compatible pairs in the cyk structure
    #~	               320	     321	194.17596	0.0153764	127	0.87
    #~	 *	        98	       106	121.80433	3.80688e-10	12	0.35
    elsif (/^~\s+.+\d+\s+\d+\s+\S+\s+\S+\s+\d+\s+(\S+)$/) {
	my $pp = $1;
	if ($iscyk) {
	    $fam_tp_cyk{$fam}    ++; 
	    $fam_true_cyk{$fam}  ++; 
	    $fam_tpexp_cyk{$fam} += $pp; 
	}
   }
}
close (FILE);

# fields:
#
# idx       1     (cyk)
# fam       2
#
# S         3      12
# PPV       4      13
# F         5      14
# Spower    6      15
#
# true      7      16
# found     8      17
# tp        9      18
# tpexp    10      19
# avgsub   11      20
#
# avgid    21
# alen     22
# nseq     23
#
my $idx         = 1;
my $ifam        = 2;

my $iS          = 3;
my $iP          = 4;
my $iF          = 5;
my $iSpower     = 6;

my $itrue       = 7;
my $ifound      = 8;
my $itp         = 9;
my $itpexp      = 10;
my $avgsub      = 11;

my $iS_c        = 12;
my $iP_c        = 13;
my $iF_c        = 14;
my $iSpower_c   = 15;

my $itrue_c     = 16;
my $ifound_c    = 17;
my $itp_c       = 18;
my $itpexp_c    = 19;
my $avgsub_c    = 20;

my $iavgid      = 21;
my $ialen       = 22;
my $inseq       = 23;

for (my $f = 0; $f < $nf; $f ++) {
    my $fam = $fam[$f];

    # recalculate in case we added "compatible" covarying pairs ~
    FUNCS::calculateF  ($fam_tp{$fam},     $fam_true{$fam},     $fam_found{$fam},     \$fam_S{$fam},     \$fam_P{$fam},     \$fam_F{$fam});
    FUNCS::calculateF  ($fam_tp_cyk{$fam}, $fam_true_cyk{$fam}, $fam_found_cyk{$fam}, \$fam_S_cyk{$fam}, \$fam_P_cyk{$fam}, \$fam_F_cyk{$fam});
    
    FUNCS::calculateSEN($fam_tpexp{$fam},     $fam_true{$fam},     \$fam_Spower{$fam});
    FUNCS::calculateSEN($fam_tpexp_cyk{$fam}, $fam_true_cyk{$fam}, \$fam_Spower_cyk{$fam});
    
    my $name = $fam; while (length($name) < 40) { $name .= " "; }
    $fam_all{$fam}  = "$name";
    
    $fam_all{$fam} .= "\t$fam_S{$fam}\t$fam_P{$fam}\t$fam_F{$fam}";
    $fam_all{$fam} .= "\t$fam_Spower{$fam}";    
    $fam_all{$fam} .= "\t$fam_true{$fam}\t$fam_found{$fam}\t$fam_tp{$fam}\t$fam_tpexp{$fam}\t$fam_avgsub{$fam}";
    
    $fam_all{$fam} .= "\t$fam_S_cyk{$fam}\t$fam_P_cyk{$fam}\t$fam_F_cyk{$fam}";
    $fam_all{$fam} .= "\t$fam_Spower_cyk{$fam}";    
    $fam_all{$fam} .= "\t$fam_true_cyk{$fam}\t$fam_found_cyk{$fam}\t$fam_tp_cyk{$fam}\t$fam_tpexp_cyk{$fam}\t$fam_avgsub_cyk{$fam}";
    
    $fam_all{$fam} .= "\t$fam_id{$fam}\t$fam_alen{$fam}\t$fam_nseq{$fam}";
}

outfile_rank       ($outfile_rank, $outfile_allfam);
outfile_nopower    ($outfile_nopower);
outfile_withpower  ($outfile_withpower);
outfile_nocov      ($outfile_nocov);
outfile_withcov    ($outfile_withcov);
outfile_outlierp   ($outfile_outlierp);
outfile_outlierc   ($outfile_outlierc);
outfile_greyzone   ($outfile_greyzone);
outfile_betterss   ($outfile_betterss);
outfile_muchbettss ($outfile_muchbettss);
outfile_worsess    ($outfile_worsess);
outfile_equalss    ($outfile_equalss);


my $plotallfamfile = "$outname.plot.allfam";
plot_allfam($plotallfamfile);
 
my $plotallvsallfile = "$outname.plot.allvsall";
plot_allvsall($plotallvsallfile);

my $plotssfile = "$outname.plot.ss";
plot_ss($plotssfile); 


sub plot_allfam {
    my ($file) = @_;

    my $pdf    = "$file.ps";
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
    print GP "set xrange [1:$nf]\n";

    $ylabel = "fraction of basepairs that covary (%)";
    print GP "set ylabel '$ylabel'\n";
    $cmd  = "";
    $cmd .= "'$outfile_allfam'      using $idx:$iSpower:(0.7) with boxes  title '' ls 1111, "; 
    $cmd .= "'$outfile_outlierp'    using $idx:$iSpower:(0.7) with boxes  title '' ls 1117, "; 
    $cmd .= "'$outfile_allfam'      using $idx:$iS            with points title '' ls 1116 "; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
   
    $ylabel = "fraction of covarying pairs in the structure (%)";
    print GP "set ylabel '$ylabel'\n";
    $cmd  = "";
    $cmd .= "'$outfile_allfam'      using $idx:$iSpower:(0.7) with boxes  title '' ls 1111, "; 
    $cmd .= "'$outfile_outlierp'    using $idx:$iSpower:(0.7) with boxes  title '' ls 1117, "; 
    $cmd .= "'$outfile_allfam'      using $idx:$iP            with points title '' ls 1118 "; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
    
   print GP "set xrange [1:1500]\n";
    $ylabel = "fraction of basepairs that covary (%)";
    print GP "set ylabel '$ylabel'\n";
    $cmd  = "";
    $cmd .= "'$outfile_allfam'      using $idx:$iSpower:(0.7) with boxes  title '' ls 1111, "; 
    $cmd .= "'$outfile_outlierp'    using $idx:$iSpower:(0.7) with boxes  title '' ls 1117, "; 
    $cmd .= "'$outfile_allfam'      using $idx:$iS            with points title '' ls 1116 "; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
   
    $ylabel = "fraction of covarying pairs in the structure (%)";
    print GP "set ylabel '$ylabel'\n";
    $cmd  = "";
    $cmd .= "'$outfile_allfam'      using $idx:$iSpower:(0.7) with boxes  title '' ls 1111, "; 
    $cmd .= "'$outfile_outlierp'    using $idx:$iSpower:(0.7) with boxes  title '' ls 1117, "; 
    $cmd .= "'$outfile_allfam'      using $idx:$iP            with points title '' ls 1116 "; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
    
    system("open $pdf\n");
    
}
sub plot_ss {
    my ($file) = @_;
    
    my $pdf    = "$file.ps";
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
    print GP "set xrange [1:$nf]\n";

    $ylabel = "fraction of covarying pairs in the structure (%)";
    print GP "set ylabel '$ylabel'\n";
    $cmd  = "";
    $cmd .= "'$outfile_allfam'     using $idx:$iSpower:(0.7) with boxes  title '' ls 1114, "; 
    #$cmd .= "'$outfile_worsess'    using $idx:$iP            with points title '' ls 1116, "; 
    #$cmd .= "'$outfile_betterss'   using $idx:$iP            with points title '' ls 1116, "; 
    $cmd .= "'$outfile_worsess'    using $idx:$iP_c          with points title '' ls 1113, "; 
    $cmd .= "'$outfile_betterss'   using $idx:$iP_c          with points title '' ls 1115 "; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
    
    $ylabel = "fraction of basepairs that covary (%)";
    print GP "set ylabel '$ylabel'\n";
    $cmd  = "";
    $cmd .= "'$outfile_allfam'     using $idx:$iSpower:(0.7) with boxes  title '' ls 1114, "; 
    #$cmd .= "'$outfile_worsess'    using $idx:$iS            with points title '' ls 1116, "; 
    #$cmd .= "'$outfile_betterss'   using $idx:$iS            with points title '' ls 1116, "; 
    $cmd .= "'$outfile_worsess'    using $idx:$iS_c          with points title '' ls 1113, "; 
    $cmd .= "'$outfile_betterss'   using $idx:$iS_c          with points title '' ls 1115 "; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
   
    system("open $pdf\n");

}

sub plot_allvsall {
    my ($file) = @_;

    my $pdf    = "$file.ps";
    my $xlabel = "Expected Covaring  (% basepairs)";
    my $ylabel = "Significantly Covaring (% basepairs)";
    my $cmd;
    
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    print GP "set output '$pdf'\n";
    print GP "set key right top\n";
    print GP "set nokey\n";
    print GP "set size square\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    
    print GP "set xrange [-0.5:100]\n";
    print GP "set yrange [-0.5:100]\n";
    $cmd = "";
    $cmd     .= "'$outfile_allfam' using $iSpower:$iS      title '' with points ls 1112"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
    
    print GP "set xrange [-0.5:20]\n";
    print GP "set yrange [-0.5:20]\n";
    $cmd = "";
    $cmd     .= "'$outfile_allfam' using $iSpower:$iS      title '' with points ls 1112"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
    
    
    system("open $pdf\n");
}

sub outfile_rank {
    my ($outfile_rank, $outfile_allfam) = @_;
    
    my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} or $fam_Spower{$b} <=> $fam_Spower{$a} } keys(%fam_S);
    #my @S_order = sort { $fam_Spower{$b} <=> $fam_Spower{$a}  } keys(%fam_S);
    
    open (OUT1, ">$outfile_rank")   || die;
    open (OUT2, ">$outfile_allfam") || die;
	
    my $cmd = "";
    $cmd .= "# family\t\t\t\t";
    $cmd .= "\tSEN\tPPV\tF";
    $cmd .= "\tSENpower";
    $cmd .= "\tTRUE\tFOUND\tTP\tTPexp\tavgsub";
    $cmd .= "\tSENcyk\tPPVcyk\tFcyk";
    $cmd .= "\tSENpowercyk";
    $cmd .= "\tTRUEcyk\tFOUNDcyk\tTPcyk\tTPexpcyk\tavgsubcyk";
    $cmd .= "\tavgid\talen\tnseq";
    
    printf OUT1 "$cmd\n";
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam = $S_order[$f];
	my $n = $f+1;
	printf OUT1 "$fam_all{$fam}\n";
	$fam_all{$fam} = $n." ".$fam_all{$fam};
	printf OUT2 "$fam_all{$fam}\n";
    }
    
    close(OUT1);  
    close(OUT2);      
}

sub outfile_nopower{
    my ($outfile_nopower) = @_;

    open (OUT, ">$outfile_nopower") || die;    
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $fam[$f];
	my $power = $fam_Spower{$fam};
	if ($power == 0) {
	    $m ++;
	    printf     "nopower $m $fam power $power\n";
	    printf OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}

sub outfile_withpower{
    my ($outfile_withpower) = @_;

    open (OUT, ">$outfile_withpower") || die;    
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $fam[$f];
	my $power = $fam_Spower{$fam};
	if ($power > 0) {
	    $m ++;
	    print     "withpower $m $fam power $power\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}

sub outfile_nocov{
    my ($outfile_nocov) = @_;

    open (OUT, ">$outfile_nocov") || die;    
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam = $fam[$f];
	my $sen = $fam_S{$fam};
	if ($sen == 0) {
	    $m ++;
	    print     "nocov $m $fam sen $sen\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}

sub outfile_withcov{
    my ($outfile_withcov) = @_;

    open (OUT, ">$outfile_withcov") || die;    
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $fam[$f];
	my $sen = $fam_S{$fam};
	if ($sen > 0) {
	    $m ++;
	    print     "cov $m $fam sen $sen\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}

sub outfile_greyzone{
    my ($outfile_greyzone) = @_;

    open (OUT, ">$outfile_greyzone") || die;    
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $fam[$f];
	my $sen   = $fam_S{$fam};
	my $power = $fam_Spower{$fam};
	if ($sen == 0 && $power > 0 && $power < 4) {
	    $m ++;
	    print     "greyzone $m $fam sen $sen power $power\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}

sub outfile_outlierp{
    my ($outfile_outlierp) = @_;

    open (OUT, ">$outfile_outlierp") || die;
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $fam[$f];
	my $sen   = $fam_S{$fam};
	my $power = $fam_Spower{$fam};
	if ($sen < 2 && $power > 10) {
	    $m ++;
	    print     "outlierp $m $fam sen $sen power $power\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}

sub outfile_outlierc{
    my ($outfile_outlierc) = @_;

    open (OUT, ">$outfile_outlierc") || die;
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $fam[$f];
	my $sen   = $fam_S{$fam};
	my $power = $fam_Spower{$fam};
	if ($sen > 10 && $power < 2) {
	    $m ++;
	    print     "outlierc $m $fam sen $sen power $power\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}


sub outfile_betterss{
    my ($outfile_betterss) = @_;

    if (!$iscyk) { return; }
    
    open (OUT, ">$outfile_betterss") || die;    
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam     = $fam[$f];
	my $sen     = $fam_S{$fam};
	my $sen_cyk = $fam_S_cyk{$fam};
	my $ppv     = $fam_P{$fam};
	my $ppv_cyk = $fam_P_cyk{$fam};
	my $tp      = $fam_tp{$fam};
	my $tp_cyk  = $fam_tp_cyk{$fam};
	my $power = $fam_Spower{$fam};
	if ($ppv_cyk > $ppv) {
	    $m ++;
	    print     "better ss $m $fam sen $sen sen_cyk $sen_cyk ppv $ppv ppv_cyk $ppv_cyk tp $tp tp_cyk $tp_cyk\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}
 
sub outfile_muchbettss{
    my ($outfile_muchbettss) = @_;

    if (!$iscyk) { return; }
    
    open (OUT, ">$outfile_muchbettss") || die;    
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam     = $fam[$f];
	my $sen     = $fam_S{$fam};
	my $sen_cyk = $fam_S_cyk{$fam};
	my $ppv     = $fam_P{$fam};
	my $ppv_cyk = $fam_P_cyk{$fam};
	my $tp      = $fam_tp{$fam};
	my $tp_cyk  = $fam_tp_cyk{$fam};
	my $power = $fam_Spower{$fam};
	if ($tp_cyk > $tp+2) {
	    $m ++;
	    print     "much better ss $m $fam sen $sen sen_cyk $sen_cyk ppv $ppv ppv_cyk $ppv_cyk tp $tp tp_cyk $tp_cyk\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}

sub outfile_worsess{
    my ($outfile_worsess) = @_;

    if (!$iscyk) { return; }
    
    open (OUT, ">$outfile_worsess") || die;    
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam     = $fam[$f];
	my $sen     = $fam_S{$fam};
	my $sen_cyk = $fam_S_cyk{$fam};
	my $ppv     = $fam_P{$fam};
	my $ppv_cyk = $fam_P_cyk{$fam};
	my $tp      = $fam_tp{$fam};
	my $tp_cyk  = $fam_tp_cyk{$fam};
	my $power = $fam_Spower{$fam};
	if ($ppv_cyk < $ppv) {
	    $m ++;
	    print "worse ss $m $fam sen $sen sen_cyk $sen_cyk ppv $ppv ppv_cyk $ppv_cyk tp $tp tp_cyk $tp_cyk\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}
 
sub outfile_equalss{
    my ($outfile_equalss) = @_;

    if (!$iscyk) { return; }
    
    open (OUT, ">$outfile_equalss") || die;    
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam     = $fam[$f];
	my $sen     = $fam_S{$fam};
	my $sen_cyk = $fam_S_cyk{$fam};
	my $ppv     = $fam_P{$fam};
	my $ppv_cyk = $fam_P_cyk{$fam};
	my $tp      = $fam_tp{$fam};
	my $tp_cyk  = $fam_tp_cyk{$fam};
	my $power = $fam_Spower{$fam};
	if ($ppv_cyk == $ppv) {
	    $m ++;
	    print "equal ss $m $fam sen $sen sen_cyk $sen_cyk ppv $ppv ppv_cyk $ppv_cyk tp $tp tp_cyk $tp_cyk\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
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
