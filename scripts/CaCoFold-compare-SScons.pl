#!/usr/bin/perl -w
# 
#  Compare the CaCoFold structure to the consensus structure provided with the alignment
#
use strict;
use Class::Struct;
use POSIX;
    
use constant GNUPLOT => '/usr/local/bin/gnuplot';
use lib '/Users/erivas/src/Mysrc/R-scape/scripts';
use FUNCS;

use vars qw ($opt_b $opt_p $opt_v $opt_f);  # required if strict used
use Getopt::Std;
getopts ('b:p:vf:');
# Print a helpful message if the user provides no input file.
if (!@ARGV) {
    print "usage:  CaCoFold-compare-SS_cons.pl [options] <R-scape.out> \n\n";
    print "options:\n";
    print "-v    :  be verbose\n";
    exit;
}
my $file    = shift;
my $filename;
my $outdir;
my $outdirname = "compare";
if    ($file =~ /^(\S+)\/([^\/]+).rscape/) { $outdir = "$1/$outdirname"; $filename = "$2"; }
elsif ($file =~ /^([^\/]+).rscape/)        { $outdir = "$outdirname";    $filename = "$1"; }
system("mkdir $outdir\n");
my $outname = "$outdir/$filename";
print "outdir $outdir\noutname $outname\n";

# treshholds
my $nopower      = 0;  # (%)
my $cov_plus     = 2;  # $tp_fold > $tp+$cov_plus
my $thresh_tpexp = 4;  # greyzone: no cov <= $thresh_cov and tpexp > $thresh_tpexp
my $thresh_cov   = 0;  # greyzone: no cov <= $thresh_cov and tpexp > $thresh_tpexp

# nbps thresholds for reporting total # families
my $min_nbps  = 1;
if ($opt_b) { $min_nbps = $opt_b; }


my $outfile_rank         = "$outname.rank";        # file with all families ranked by Sensitivity/power 
my $outfile_allfam       = "$outname.allfam";      # file with all families ranked by Sensitivity to plot
my $outfile_nopower      = "$outname.nopower";     # families with no power (power = 0) 
my $outfile_withpower    = "$outname.withpower";   # families with power   (power > 0)
my $outfile_nocov        = "$outname.nocov";       # families with no significant covariations   (sen = 0)
my $outfile_withcov      = "$outname.withcov";     # families with    significant covariations   (sen > 0)
my $outfile_greyzone     = "$outname.greyzone";
my $outfile_outlierp     = "$outname.outlierp";
my $outfile_outlierc     = "$outname.outlierc";
my $outfile_nopnoc       = "$outname.nopnoc";
my $outfile_nopc         = "$outname.nopc";
my $outfile_pnoc         = "$outname.pnoc";
my $outfile_goodpower    = "$outname.goodpower";
my $outfile_betterss     = "$outname.betterss";
my $outfile_muchbettss   = "$outname.muchbetterss";
my $outfile_worsess      = "$outname.worsess";
my $outfile_equalss      = "$outname.equalss";


my @allfam;
my %allfam_idx;
my %allfam_id;
my %allfam_alen;
my %allfam_nseq;

my %allfam_tp;
my %allfam_true;
my %allfam_found;
my %allfam_tpexp;
my %allfam_avgsub;

my %allfam_S;
my %allfam_P;
my %allfam_F;
my %allfam_Spower;
 
my %allfam_tp_fold;
my %allfam_true_fold;
my %allfam_found_fold;
my %allfam_tpexp_fold;
my %allfam_avgsub_fold;

my %allfam_S_fold;
my %allfam_P_fold;
my %allfam_F_fold;
my %allfam_Spower_fold;

my %allfam_all;
my %allfam_table;
my %allfam_tabless;


# fields:
#
# idx       1     (fold)
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

my $nfall = parse_rscapeout($file);
print "NALLFAM $nfall\n";

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

my %fam_tp_fold;
my %fam_true_fold;
my %fam_found_fold;
my %fam_tpexp_fold;
my %fam_avgsub_fold;

my %fam_S_fold;
my %fam_P_fold;
my %fam_F_fold;
my %fam_Spower_fold;
my %fam_all;
my %fam_table;
my %fam_tabless;

my $nofilter = 1;
my $nf = filter_fam($nofilter);

outfile_rank       ($outfile_rank, $outfile_allfam);
outfile_nopower    ($nopower, $outfile_nopower);
outfile_withpower  ($nopower, $outfile_withpower);
outfile_nocov      ($outfile_nocov);
outfile_withcov    ($outfile_withcov);
outfile_outlierp   ($thresh_tpexp, $outfile_outlierp);
outfile_outlierc   ($thresh_cov, $thresh_tpexp, $outfile_outlierc);
outfile_nopnoc     ($outfile_nopnoc);
outfile_nopc       ($outfile_nopc);
outfile_pnoc       ($outfile_pnoc);
outfile_greyzone   ($thresh_tpexp, $outfile_greyzone);
outfile_betterss   ($outfile_betterss);
outfile_muchbettss ($cov_plus, $outfile_muchbettss);
outfile_worsess    ($outfile_worsess);
outfile_equalss    ($outfile_equalss);


my $plotallvsallfile = "$outname.plot.allvsall";
plot_allvsall($plotallvsallfile);



sub parse_rscapeout {
    my  ($file) = @_;

    my $usefam;
    my $isfold;
    my $fam;

    my $tp      = $1;
    my $true    = $2;
    my $found   = $3;
    
    my $sen     = $4;
    my $ppv     = $5;
    my $F       = $6;

    my $nmode;
    my $nf  = 0;
    open (FILE, "$file")    || die;
    while(<FILE>) {
	# MSA RF00001_5S_rRNA nseq 712 (712) alen 119 (230) avgid 56.09 (55.86) nbpairs 34 (34)
	if (/^# MSA\s+(\S+)\s+nseq\s+(\S+)\s+\(.+alen\s+(\S+)\s+\(.+avgid\s+(\S+)\s+/) {
	    $fam      = $1;
	    my $nseq  = $2;
	    my $alen  = $3;
	    my $avgid = $4;
	    my $acc   = $fam;
	    if ($acc =~ /^(RF\d\d\d\d\d)/) { $acc = $1; }
	    #if ($fam =~ /^RF\d\d\d\d\d\_(\S+)$/) { $fam = $1; }
	    
	    $usefam = 1;
	    
	    $isfold  = 0;
	    $nmode   = 0;
	    
	    if ($usefam) {
		$allfam[$nf]       = $fam; 
		$allfam_idx{$fam}  = $nf+1;
		$allfam_nseq{$fam} = $nseq;
		$allfam_alen{$fam} = $alen;
		$allfam_id{$fam}   = $avgid;
		
		$allfam_tp{$fam}         = 0.0;
		$allfam_true{$fam}       = 0.0;
		$allfam_found{$fam}      = 0.0;
		$allfam_tpexp{$fam}      = 0.0;
		
		$allfam_S{$fam}          = 0.0;
		$allfam_P{$fam}          = 0.0;
		$allfam_F{$fam}          = 0.0;
		
		$allfam_Spower{$fam}     = 0.0;
		
		$allfam_tp_fold{$fam}     = 0.0;
		$allfam_true_fold{$fam}   = 0.0;
  		$allfam_found_fold{$fam}  = 0.0;
		$allfam_tpexp_fold{$fam}  = 0.0;
		
		$allfam_S_fold{$fam}      = 0.0;
		$allfam_P_fold{$fam}      = 0.0;
		$allfam_F_fold{$fam}      = 0.0;

		$allfam_avgsub{$fam}      = 0.0;
		$allfam_avgsub_fold{$fam} = 0.0;
		$allfam_tpexp{$fam}       = 0.0;
		$allfam_tpexp_fold{$fam}  = 0.0;
		
		$allfam_Spower_fold{$fam} = 0.0;

	    }
	}
	# Method Target_E-val [cov_min,conv_max] [FP | TP True Found | Sen PPV F]
	elsif (/^# Method/) {
	    $nmode ++;
	    if ($nmode == 2) { $isfold = 1; }
	}
	# GTp    0.05         [-9.35,1389.53]     [8 | 22 34 30 | 64.71 73.33 68.75] 
	elsif (/^#\s+\S+\s+\S+\s+\[\S+\]\s+\[\d+\s+\|\s+(\d+)\s+(\d+)\s+(\d+)\s+\|\s+(\S+)\s+(\S+)\s+(\S+)\]/) {
	    $tp      = $1;
	    $true    = $2;
	    $found   = $3;
	    
	    $sen     = $4;
	    $ppv     = $5;
	    $F       = $6;
	}
	# BPAIRS 20
	# avg substitutions per BP 15.8
	# BPAIRS expected to covary 7.1
	elsif (/^#\s+avg substitutions per BP\s+(\S+)/) {
	    my $avgsub = $1;
	    if ($usefam) {
		if ($isfold) { $allfam_avgsub_fold{$fam} = $avgsub; }
		else         { $allfam_avgsub{$fam}      = $avgsub; }
	    }
	}
	elsif (/^#\s+BPAIRS expected to covary\s+(\S+)/) {
	    my $tpexp = $1;
	    if ($usefam) { 
		if ($isfold) { $allfam_tpexp_fold{$fam} = $tpexp; }
		else         { $allfam_tpexp{$fam}      = $tpexp; }
	    }

	    if ($usefam) {	    
		if ($nmode == 2) {
		    printf "%d FAM CaCoFold $allfam[$nf-1] bp $true bp&cov $tp cov $found exp_cov $tpexp\n", $nf;
		    $allfam_tp_fold{$fam}    = $tp;
		    $allfam_true_fold{$fam}  = $true;
		    $allfam_found_fold{$fam} = $found;
		}
		else {
		    printf "\n%d FAM SScons   $allfam[$nf] bp $true bp&cov $tp cov $found exp_cov $tpexp\n", $nf+1;
		    $allfam_tp{$fam}    = $tp;
		    $allfam_true{$fam}  = $true;
		    $allfam_found{$fam} = $found;
		    $nf ++;
		}
	    }
	}
    }
    close (FILE);
    print "NFAMILIES parsed $nf\n";
       
    # check
    for (my $f = 0; $f < $nf-1; $f ++) {
	my $fam = $allfam[$f];
	for (my $ff = $f+1; $ff < $nf; $ff ++) {
	    my $ffam = $allfam[$ff];
	    if ($ffam =~ /^$fam$/) { printf("^^AH $f $fam $ff $ffam\n"); die; }
	}
    }
    
    # cummulatives
    my $tp_tot         = 0;
    my $true_tot       = 0;
    my $found_tot      = 0;
    my $tp_fold_tot    = 0;
    my $true_fold_tot  = 0;
    my $found_fold_tot = 0;
    my $nf_used = 0;
   
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam = $allfam[$f];

	if ($allfam_true{$fam}  >= $min_nbps) 
	{
	    $tp_tot         += $allfam_tp{$fam};
	    $true_tot       += $allfam_true{$fam};
	    $found_tot      += $allfam_found{$fam};
	    $tp_fold_tot    += $allfam_tp_fold{$fam};
	    $true_fold_tot  += $allfam_true_fold{$fam};
	    $found_fold_tot += $allfam_found_fold{$fam};
	    $nf_used ++;
	}
    }
    printf("Totals for %d/$nf families (min_bp = $min_nbps)\n               SScons    CaCoFOld\n", $nf_used, $nf);
    printf("cov_bps        %d %d\n", $tp_tot,    $tp_fold_tot);
    printf("all_bps        %d %d\n", $true_tot,  $true_fold_tot);
    
 
	
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam = $allfam[$f];
	
	# recalculate in case we added "compatible" covarying pairs ~
	FUNCS::calculateF($allfam_tp{$fam},      $allfam_true{$fam},      $allfam_found{$fam},      \$allfam_S{$fam},      \$allfam_P{$fam},      \$allfam_F{$fam});
	FUNCS::calculateF($allfam_tp_fold{$fam}, $allfam_true_fold{$fam}, $allfam_found_fold{$fam}, \$allfam_S_fold{$fam}, \$allfam_P_fold{$fam}, \$allfam_F_fold{$fam});
	
	FUNCS::calculateSEN($allfam_tpexp{$fam},      $allfam_true{$fam},      \$allfam_Spower{$fam});
	FUNCS::calculateSEN($allfam_tpexp_fold{$fam}, $allfam_true_fold{$fam}, \$allfam_Spower_fold{$fam});

	$allfam_Spower{$fam}  = sprintf("%.1f", $allfam_Spower{$fam});
	$allfam_S{$fam}       = sprintf("%.1f", $allfam_S{$fam});
	$allfam_P{$fam}       = sprintf("%.1f", $allfam_P{$fam});
 	$allfam_F{$fam}       = sprintf("%.1f", $allfam_F{$fam});
	$allfam_S_fold{$fam}  = sprintf("%.1f", $allfam_S_fold{$fam});
	$allfam_P_fold{$fam}  = sprintf("%.1f", $allfam_P_fold{$fam});
	$allfam_F_fold{$fam}  = sprintf("%.1f", $allfam_F_fold{$fam});
	$allfam_id{$fam}      = sprintf("%.1f", $allfam_id{$fam});
	    
	my $name = $fam;
	$name =~ s/\_/ /g;
	
	my $longname = $fam; $longname =~ s/\_/-/g;
	while (length($longname) < 40) { $longname .= " "; }
	
	$allfam_all{$fam}  = "$longname";	
	$allfam_all{$fam} .= "\t$allfam_S{$fam}\t$allfam_P{$fam}\t$allfam_F{$fam}";
	$allfam_all{$fam} .= "\t$allfam_Spower{$fam}";    
	$allfam_all{$fam} .= "\t$allfam_true{$fam}\t$allfam_found{$fam}\t$allfam_tp{$fam}\t$allfam_tpexp{$fam}\t$allfam_avgsub{$fam}";
	
	$allfam_all{$fam} .= "\t$allfam_S_fold{$fam}\t$allfam_P_fold{$fam}\t$allfam_F_fold{$fam}";
	$allfam_all{$fam} .= "\t$allfam_Spower_fold{$fam}";    
	$allfam_all{$fam} .= "\t$allfam_true_fold{$fam}\t$allfam_found_fold{$fam}\t$allfam_tp_fold{$fam}\t$allfam_tpexp_fold{$fam}\t$allfam_avgsub_fold{$fam}";
	
	$allfam_all{$fam} .= "\t$allfam_id{$fam}\t$allfam_alen{$fam}\t$allfam_nseq{$fam}";

	$allfam_table{$fam}  = "$name & ";
	$allfam_table{$fam} .= "$allfam_S{$fam} ($allfam_tp{$fam}/$allfam_true{$fam}) & ";    
	$allfam_table{$fam} .= "$allfam_Spower{$fam} & ";    
	$allfam_table{$fam} .= "$allfam_P{$fam} ($allfam_tp{$fam}/$allfam_found{$fam}) & ";    
	$allfam_table{$fam} .= "$allfam_avgsub{$fam} & ";    
	$allfam_table{$fam} .= "$allfam_id{$fam} & ";    
	$allfam_table{$fam} .= "$allfam_nseq{$fam}\n\\\\"; 
   
	$allfam_tabless{$fam}  = "$name & ";
	$allfam_tabless{$fam} .= "$allfam_tp_fold{$fam}    & $allfam_tp{$fam}    & ";    
	$allfam_tabless{$fam} .= "$allfam_true_fold{$fam}  & $allfam_true{$fam}  & ";    
	$allfam_tabless{$fam} .= "$allfam_S_fold{$fam}     & $allfam_S{$fam}     & ";    
	$allfam_tabless{$fam} .= "$allfam_found{$fam} & ";    
	$allfam_tabless{$fam} .= "$allfam_P_fold{$fam}     & $allfam_P{$fam}     & ";    
	$allfam_tabless{$fam} .= "$allfam_tpexp_fold{$fam} & $allfam_tpexp{$fam} ";    
	$allfam_tabless{$fam} .= "\n\\\\";    
    }

    return $nf;
}

sub filter_fam {
    my ($nofilter) = @_;
    
    my $nf = 0;

    for (my $f = 0; $f < $nfall; $f ++) {

	my $fam = $allfam[$f];

	if ($nofilter || (!$nofilter && $allfam_true{$fam}  > 10)) {
	    
	    $fam[$nf] = $fam;
	    
	    $fam_idx{$fam}  = $nf + 1;
	    $fam_id{$fam}   = $allfam_id{$fam};
	    $fam_alen{$fam} = $allfam_alen{$fam};
	    $fam_nseq{$fam} = $allfam_nseq{$fam};
	    
	    $fam_tp{$fam}     = $allfam_tp{$fam};
	    $fam_true{$fam}   = $allfam_true{$fam};
	    $fam_found{$fam}  = $allfam_found{$fam};
	    $fam_tpexp{$fam}  = $allfam_tpexp{$fam};
	    $fam_avgsub{$fam} = $allfam_avgsub{$fam};
	    
	    $fam_S{$fam}      = $allfam_S{$fam};
	    $fam_P{$fam}      = $allfam_P{$fam};
	    $fam_F{$fam}      = $allfam_F{$fam};
	    $fam_Spower{$fam} = $allfam_Spower{$fam};
	    
	    $fam_tp_fold{$fam}     = $allfam_tp_fold{$fam};
	    $fam_true_fold{$fam}   = $allfam_true_fold{$fam};
	    $fam_found_fold{$fam}  = $allfam_found_fold{$fam};
	    $fam_tpexp_fold{$fam}  = $allfam_tpexp_fold{$fam};
	    $fam_avgsub_fold{$fam} = $allfam_avgsub_fold{$fam};
	    
	    $fam_S_fold{$fam}      = $allfam_S_fold{$fam};
	    $fam_P_fold{$fam}      = $allfam_P_fold{$fam};
	    $fam_F_fold{$fam}      = $allfam_F_fold{$fam};
	    $fam_Spower_fold{$fam} = $allfam_Spower_fold{$fam};
	    
	    $fam_all{$fam}     = $allfam_all{$fam};
	    $fam_table{$fam}   = $allfam_table{$fam};
	    $fam_tabless{$fam} = $allfam_tabless{$fam};
	    
	    $nf ++;
	}
	
    }
    print "$nf FILTERED FAM\n";
    return $nf;
}


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
    #$cmd .= "'$outfile_allfam'      using $idx:$iP            with points title '' ls 1118,  "; 
    $cmd .= "'$outfile_allfam'      using $idx:$iS            with points title '' ls 1114 "; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
   
    print GP "set xrange [1:1031]\n";
    $ylabel = "fraction of covarying pairs in the structure (%)";
    print GP "set ylabel '$ylabel'\n";
    $cmd  = "";
    $cmd .= "'$outfile_allfam'      using $idx:$iSpower:(0.7) with boxes  title '' ls 1111, "; 
    $cmd .= "'$outfile_outlierp'    using $idx:$iSpower:(0.7) with boxes  title '' ls 1117, "; 
    $cmd .= "'$outfile_allfam'      using $idx:$iS            with points title '' ls 1114, "; 
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

sub plot_allvsall {
    my ($file) = @_;

    my $pdf    = "$file.ps";
    my $xlabel = "base pairs expected to covary";
    my $ylabel = "base pairs observed to covary";
    my $cmd;
    
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    print GP "set output '$pdf'\n";
    print GP "set key right top\n";
    #print GP "set nokey\n";
    print GP "set size square\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    
    print GP "set xrange [-0.5:*]\n";
    print GP "set yrange [-0.5:*]\n";
    $cmd  = "";
    $cmd .= "'$outfile_allfam' using $itpexp:$itp    title 'SScons' with points ls 1114"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
    
    print GP "set xrange [-0.5:*]\n";
    print GP "set yrange [-0.5:*]\n";
    $cmd  = "";
    $cmd .= "'$outfile_allfam' using $itpexp_c:$itp_c    title 'CaCoFold' with points ls 1112"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
    
    print GP "set xrange [-0.5:*]\n";
    print GP "set yrange [-0.5:*]\n";
    $cmd  = "";
    $cmd .= "'$outfile_allfam' using $itpexp_c:$itp_c:$ifam      title ''         with labels ls 1112, "; 
    $cmd .= "'$outfile_allfam' using $itpexp_c:$itp_c            title 'CaCoFold' with points ls 1112"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";

    print GP "set xrange [-0.5:50]\n";
    print GP "set yrange [-0.5:50]\n";
    $cmd  = "";
    $cmd .= "'$outfile_allfam' using $itpexp_c:$itp_c      title 'CaCoFold' with points ls 1112"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
    
    print GP "set xrange [-0.5:50]\n";
    print GP "set yrange [-0.5:50]\n";
    $cmd  = "";
    $cmd .= "'$outfile_allfam' using $itpexp_c:$itp_c:$ifam      title ''         with labels ls 1112, "; 
    $cmd .= "'$outfile_allfam' using $itpexp_c:$itp_c            title 'CaCoFold' with points ls 1112"; 
    $cmd .= "\n";
    print GP "plot $cmd\n";
    
    
    system("open $pdf\n");
}



sub outfile_rank {
    my ($outfile_rank, $outfile_allfam) = @_;
    
    #my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} } keys(%fam_S);
    my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} or $fam_Spower{$b} <=> $fam_Spower{$a} } keys(%fam_S);
   
    open (OUT1,  ">$outfile_rank")      || die;
    open (OUT2,  ">$outfile_allfam")    || die;
 
    my $cmd = "";
    $cmd .= "# family\t\t\t\t";
    $cmd .= "\tSEN\tPPV\tF";
    $cmd .= "\tSENpower";
    $cmd .= "\tTRUE\tFOUND\tTP\tTPexp\tavgsub";
    $cmd .= "\tSENfold\tPPVfold\tFfold";
    $cmd .= "\tSENpowerfold";
    $cmd .= "\tTRUEfold\tFOUNDfold\tTPfold\tTPexpfold\tavgsubfold";
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
    my ($nopower, $outfile_nopower) = @_;

    open (OUT, ">$outfile_nopower") || die;
    printf OUT "\nFamilies with power <= $nopower\n";
    printf     "\nFamilies with power <= $nopower\n";
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $fam[$f];
	my $power = $fam_Spower{$fam};
	my $tp_fold  = $fam_tp_fold{$fam};

	if ($power <= $nopower) {
	    $m ++;
	    printf     "nopower $m $fam power $power cov pairs $tp_fold\n";
	    printf OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}

sub outfile_withpower{
    my ($nopower, $outfile_withpower) = @_;

    open (OUT, ">$outfile_withpower") || die;    
    printf OUT "\nFamilies with power > $nopower\n";
    printf     "\nFamilies with power > $nopower\n";
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $fam[$f];
	my $power = $fam_Spower{$fam};
	my $sen   = $fam_S{$fam};
	if ($power > $nopower) {
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
    printf OUT "\nFamilies without covariation\n";
    printf     "\nFamilies without covariation\n";
 
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
    printf OUT "\nFamilies with covariation\n";
    printf     "\nFamilies with covariation\n";
  
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
    my ($thresh_tpexp, $outfile_greyzone) = @_;

    open (OUT, ">$outfile_greyzone") || die; 
    printf     "greyzone no cov and tpexp <= %d\n", $thresh_tpexp;
    printf OUT "greyzone no cov and tpexp <= %d\n", $thresh_tpexp;
   
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $fam[$f];
	my $sen   = $fam_S{$fam};
	my $tpexp = $fam_tpexp{$fam};
	if ($sen == 0 && $tpexp <= $thresh_tpexp) {
	    $m ++;
	    print     "greyzone $m $fam sen $sen tpexp $tpexp\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}

sub outfile_outlierp{
    my ($thresh_tpexp, $outfile_outlierp) = @_;

    my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} or $fam_Spower{$b} <=> $fam_Spower{$a} } keys(%fam_S);
 
    open (OUT, ">$outfile_outlierp") || die;
    printf OUT "\nFamilies without cov but with tpexp > %d\n", $thresh_tpexp;
    printf     "\nFamilies without cov but with tpexp > %d\n", $thresh_tpexp;

    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $S_order[$f];
	my $sen   = $fam_S{$fam};
	my $power = $fam_Spower{$fam};
	my $tpexp = $fam_tpexp{$fam};
	if ($sen == 0  && $tpexp > $thresh_tpexp) {
	    $m ++;
	    print     "outlierp $m $fam sen $sen tpexp $tpexp\n";
	    printf OUT "$fam_table{$fam} %f\n", $power*$fam_true{$fam}/100;
	}
    }
    close(OUT);
}

sub outfile_outlierc{
    my ($thresh_cov, $thresh_tpexp, $outfile_outlierc) = @_;

    my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} or $fam_Spower{$b} <=> $fam_Spower{$a} } keys(%fam_S);

    open (OUT, ">$outfile_outlierc") || die;
    printf OUT "\nFamilies with cov > %d but tpexp <= %d\n", $thresh_cov, $thresh_tpexp;
    printf     "\nFamilies with cov > %d but tpexp <= %d\n", $thresh_cov, $thresh_tpexp;

    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $S_order[$f];
	my $sen   = $fam_S{$fam};
	my $power = $fam_Spower{$fam};
	my $tpexp = $fam_tpexp{$fam};
	
	if ($sen > $thresh_cov && $tpexp <= $thresh_tpexp) {
	    $m ++;
	    print     "outlierc $m $fam sen $sen tpexp $tpexp\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}
sub outfile_nopnoc{
    my ($outfile_nopnoc) = @_;

    my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} or $fam_Spower{$b} <=> $fam_Spower{$a} } keys(%fam_S);

    open (OUT, ">$outfile_nopnoc") || die;
    printf OUT "\nFamilies without cov and power\n";
    printf     "\nFamilies without cov and power \n";

    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $S_order[$f];
	my $sen   = $fam_S{$fam};
	my $power = $fam_Spower{$fam};
	if ($sen == 0 && $power == 0) {
	    $m ++;
	    print     "nopnoc $m $fam sen $sen power $power\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}
sub outfile_nopc{
    my ($outfile_nopc) = @_;

    my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} or $fam_Spower{$b} <=> $fam_Spower{$a} } keys(%fam_S);

    open (OUT, ">$outfile_nopc") || die;
    printf OUT "\nFamilies with cov > 0 but power = 0\n";
    printf     "\nFamilies with cov > 0 but power = 0\n";

    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $S_order[$f];
	my $sen   = $fam_S{$fam};
	my $power = $fam_Spower{$fam};
	if ($sen > 0 && $power == 0) {
	    $m ++;
	    print     "nopc $m $fam sen $sen power $power\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}

sub outfile_pnoc{
    my ($outfile_pnoc) = @_;

    my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} or $fam_Spower{$b} <=> $fam_Spower{$a} } keys(%fam_S);

    open (OUT, ">$outfile_pnoc") || die;
    printf OUT "\nFamilies with cov = 0 but power > 0\n";
    printf     "\nFamilies with cov = 0 but power > 0\n";

    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam   = $S_order[$f];
	my $sen   = $fam_S{$fam};
	my $power = $fam_Spower{$fam};
	if ($sen == 0 && $power > 0) {
	    $m ++;
	    print     "pnoc $m $fam sen $sen power $power\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}


sub outfile_betterss{
    my ($outfile_betterss) = @_;
    
    my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} or $fam_Spower{$b} <=> $fam_Spower{$a} } keys(%fam_S);

    open (OUT, ">$outfile_betterss") || die;
    printf OUT "\nFamilies with better ss: CaCoFold ss has more covarying pairs than Rfam SScons\n";
    printf     "\nFamilies with better ss: CaCoFold ss has more covarying pairs than Rfam SScons\n";

    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam      = $S_order[$f];
	my $which    = $fam_idx{$fam};
	my $sen      = $fam_S{$fam};
	my $sen_fold = $fam_S_fold{$fam};
	my $ppv      = $fam_P{$fam};
	my $ppv_fold = $fam_P_fold{$fam};
	my $tp       = $fam_tp{$fam};
	my $tp_fold  = $fam_tp_fold{$fam};
	my $power    = $fam_Spower{$fam};
	if ($tp_fold > $tp) {
	    $m ++;
	    print     "better_ss $m $which $fam sen $sen sen_fold $sen_fold ppv $ppv ppv_fold $ppv_fold tp $tp tp_fold $tp_fold\n";
	    print OUT "$fam_tabless{$fam}\n";
	}
    }
    close(OUT);
}
 
sub outfile_muchbettss{
    my ($cov_plus, $outfile_muchbettss) = @_;

    my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} or $fam_Spower{$b} <=> $fam_Spower{$a} } keys(%fam_S);

    open (OUT, ">$outfile_muchbettss") || die;
    printf OUT "\nFamilies with much_better ss: CaCoFold ss has at least %d covarying pairs than Rfam SScons\n", $cov_plus;
    printf     "\nFamilies with much_better ss: CaCoFold ss has at least %d covarying pairs than Rfam SScons\n", $cov_plus;

    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam     = $S_order[$f];
	my $which   = $fam_idx{$fam};
	my $sen     = $fam_S{$fam};
	my $sen_fold = $fam_S_fold{$fam};
	my $ppv     = $fam_P{$fam};
	my $ppv_fold = $fam_P_fold{$fam};
	my $tp      = $fam_tp{$fam};
	my $tp_fold  = $fam_tp_fold{$fam};
	my $power   = $fam_Spower{$fam};
	if ($tp_fold > $tp+$cov_plus) {
	    $m ++;
	    print     "much_better_ss $m $which $fam sen $sen sen_fold $sen_fold ppv $ppv ppv_fold $ppv_fold tp $tp tp_fold $tp_fold\n";
	    print OUT "$fam_tabless{$fam}\n";
	}
    }
    close(OUT);
}

sub outfile_worsess{
    my ($outfile_worsess) = @_;
    
    my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} or $fam_Spower{$b} <=> $fam_Spower{$a} } keys(%fam_S);

    open (OUT, ">$outfile_worsess") || die;    
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam      = $S_order[$f];
	my $which    = $fam_idx{$fam};
	my $sen      = $fam_S{$fam};
	my $sen_fold = $fam_S_fold{$fam};
	my $ppv      = $fam_P{$fam};
	my $ppv_fold = $fam_P_fold{$fam};
	my $tp       = $fam_tp{$fam};
	my $tp_fold  = $fam_tp_fold{$fam};
	my $power    = $fam_Spower{$fam};
	if ($ppv_fold < $ppv) {
	    $m ++;
	    print "worse_ss $m $which $fam sen $sen sen_fold $sen_fold ppv $ppv ppv_fold $ppv_fold tp $tp tp_fold $tp_fold\n";
	    print OUT "$fam_tabless{$fam}\n";
	}
    }
    close(OUT);
}
 
sub outfile_equalss{
    my ($outfile_equalss) = @_;
    
    open (OUT, ">$outfile_equalss") || die;  
    printf OUT "\nFamilies with equal ss: CaCoFold ss has same covarying pairs than Rfam SScons\n";
    #printf     "\nFamilies with equal ss: CaCoFold ss has same covarying pairs than Rfam SScons\n";
  
    my $m = 0;    
    for (my $f = 0; $f < $nf; $f ++) {
	my $fam      = $fam[$f];
	my $which    = $fam_idx{$fam};
	my $sen      = $fam_S{$fam};
	my $sen_fold = $fam_S_fold{$fam};
	my $ppv      = $fam_P{$fam};
	my $ppv_fold = $fam_P_fold{$fam};
	my $tp       = $fam_tp{$fam};
	my $tp_fold  = $fam_tp_fold{$fam};
	my $power = $fam_Spower{$fam};
	if ($ppv_fold == $ppv) {
	    $m ++;
	    #print "equal_ss $m $which $fam sen $sen sen_fold $sen_fold ppv $ppv ppv_fold $ppv_fold tp $tp tp_fold $tp_fold\n";
	    print OUT "$fam_all{$fam}\n";
	}
    }
    close(OUT);
}
 



