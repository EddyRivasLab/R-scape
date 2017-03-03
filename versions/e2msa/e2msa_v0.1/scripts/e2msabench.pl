#!/usr/bin/perl -w
#e2msabench.pl
 
use strict;
use Class::Struct;
use lib '/groups/eddy/home/rivase/projects/evohmm/scripts';
use lib '/Users/rivase/projects/evohmm/scripts';
use FUNCS;
#use constant GNUPLOT => '/usr/bin/gnuplot';
use constant GNUPLOT => '/opt/local/bin/gnuplot';
#use constant GNUPLOT => '/sw/bin/gnuplot';

use vars qw ($opt_b $opt_D $opt_E $opt_i $opt_I $opt_p $opt_s $opt_t);  # required if strict used
use Getopt::Std;
getopts ('b:D:Ei:I:pts:');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  e2msabench.pl [options] N <benchfile1>..<benchfileN>  \n\n";
        print "options:\n";
	exit;
}


my $N = shift;
my @file;
my $name = "";
my $filename;
my @filename;
for (my $n = 0; $n < $N-1; $n ++) {
    $file[$n] = shift;
    $filename = $file[$n];
    if ($filename =~ /^(fix\S+)\/([^\/]+).bench$/)         { $filename = "$1_$2"; }
    if ($filename =~ /^(slower)\/([^\/]+).bench\d*$/)      { $filename = "$1_$2"; }
    if ($filename =~ /^(pam\S+)\/([^\/]+).bench\d*$/)      { $filename = "$1_$2"; }
    if ($filename =~ /^(blosum\S+)\/([^\/]+).bench\d*$/)   { $filename = "$1_$2"; }
    if ($filename =~ /^(super\S+)\/([^\/]+).bench\d*$/)    { $filename = "$1_$2"; }
    if ($filename =~ /^(msaprobs)\/([^\/]+).bench$/)       { $filename = "$1_$2"; }
    if ($filename =~ /^(muscle)\/([^\/]+).bench$/)         { $filename = "$1_$2"; }
    if ($filename =~ /^(\S*blast)\/([^\/]+).bench$/)       { $filename = "$1_$2"; }
    if ($filename =~ /^(\S*blast-ws4)\/([^\/]+).bench$/)   { $filename = "$1_$2"; }
    if ($filename =~ /^(optimal)\/([^\/]+).bench$/)        { $filename = "$1_$2"; }
    if ($filename =~ /^(phmmer)\/([^\/]+).bench$/)         { $filename = "$1_$2"; }
    if ($filename =~ /^(phmmer3)\/([^\/]+).bench$/)        { $filename = "$1_$2"; }
    if ($filename =~ /^(ephmmer)\/([^\/]+).bench$/)        { $filename = "$1_$2"; }
    if ($filename =~ /^(phmmer-max)\/([^\/]+).bench$/)     { $filename = "$1_$2"; }
    if ($filename =~ /^(phmmer3-max)\/([^\/]+).bench$/)    { $filename = "$1_$2"; }
    if ($filename =~ /^(ephmmer-max\S*)\/([^\/]+).bench$/) { $filename = "$1_$2"; }
    if ($filename =~ /^(ssearch\S*)\/([^\/]+).bench$/)     { $filename = "$1_$2"; }
    if ($filename =~ /\/([^\/]+)\/([^\/]+).bench$/)        { $filename = "$1_$2"; }
    $name .= "$filename-";
    $filename[$n] = $filename;
    print "$filename\n";
}
$file[$N-1] = shift;
$filename = $file[$N-1];
if ($filename =~ /^(fix\S+)\/([^\/]+).bench$/)         { $filename = "$1_$2"; }
if ($filename =~ /^(slower)\/([^\/]+).bench\d*$/)      { $filename = "$1_$2"; }
if ($filename =~ /^(pam\S+)\/([^\/]+).bench\d*$/)      { $filename = "$1_$2"; }
if ($filename =~ /^(blosum\S+)\/([^\/]+).bench\d*$/)   { $filename = "$1_$2"; }
if ($filename =~ /^(super\S+)\/([^\/]+).bench\d*$/)    { $filename = "$1_$2"; }
if ($filename =~ /^(msaprobs)\/([^\/]+).bench$/)       { $filename = "$1_$2"; }
if ($filename =~ /^(muscle)\/([^\/]+).bench$/)         { $filename = "$1_$2"; }
if ($filename =~ /^(\S*blast)\/([^\/]+).bench$/)       { $filename = "$1_$2"; }
if ($filename =~ /^(\S*blast-ws4)\/([^\/]+).bench$/)   { $filename = "$1_$2"; }
if ($filename =~ /^(optimal)\/([^\/]+).bench$/)        { $filename = "$1_$2"; }
if ($filename =~ /^(phmmer)\/([^\/]+).bench$/)         { $filename = "$1_$2"; }
if ($filename =~ /^(phmmer3)\/([^\/]+).bench$/)        { $filename = "$1_$2"; }
if ($filename =~ /^(ephmmer)\/([^\/]+).bench$/)        { $filename = "$1_$2"; }
if ($filename =~ /^(phmmer-max)\/([^\/]+).bench$/)     { $filename = "$1_$2"; }
if ($filename =~ /^(phmmer3-max)\/([^\/]+).bench$/)    { $filename = "$1_$2"; }
if ($filename =~ /^(ephmmer-max\S*)\/([^\/]+).bench$/) { $filename = "$1_$2"; }
if ($filename =~ /^(ssearch\S*)\/([^\/]+).bench$/)     { $filename = "$1_$2"; }
if ($filename =~ /\/([^\/]+)\/([^\/]+).bench$/)        { $filename = "$1_$2"; }
 $name .= "$filename";
$filename[$N-1] = $filename;
    print "$filename\n";

my $minid = 0.0;
my $maxid = 100.0;
my $binsize = 10.0;
if ($opt_b) { $binsize = $opt_b; if ($binsize <= 0.0) { print "ilegal binsize $binsize\n"; } }
if ($opt_i) { $minid   = $opt_i; if ($minid   <  0.0) { print "ilegal minid $minid\n"; } }
if ($opt_I) { $maxid   = $opt_I; if ($maxid   <  0.0) { print "ilegal maxid $maxid\n"; } }
my $shift = 0.0;
if ($opt_s) { $shift = $opt_s; if ($shift <= 0.0) { print "ilegal shift $shift\n"; } }

my $min_nmsa = 1; 
my $summaryfile;
my $NS;
my $ndom;
my @name;
my @d1pid;
my @d2pid;
my @d1n;
my @d2n;
my @len;
my $ave_dn = 0;
my $std_dn = 0;

my $lentot = 0;
my $homtot = 0;
my $homfrac;
my $ave_homfrac = 0;
my $std_homfrac = 0;
my $ntot = 0;
if ($opt_D) { 
    $summaryfile = "$opt_D"; 
    $ndom = parse_summary($summaryfile, \$NS, \@name, \@d1pid, \@d2pid, \@d1n, \@d2n, \@len);
    print "\nNS=$NS ndomains = $ndom\n";

    for (my $n = 0; $n < $NS; $n ++) {
	$ave_dn += $d1n[$n];
	$std_dn += $d1n[$n] * $d1n[$n];

	$lentot += $len[$n];

	$homtot += $d1n[$n];
	$homfrac = ($len[$n] > 0)? $d1n[$n] / $len[$n] : 0;

	$ave_homfrac += $homfrac;
	$std_homfrac += $homfrac * $homfrac;
	$ntot ++;

	if ($ndom == 2) { 
	    $ave_dn += $d2n[$n]; 
	    $std_dn += $d2n[$n] * $d2n[$n]; 

	    $homtot += $d1n[$n];
	    $homfrac = ($len[$n] > 0)? $d1n[$n] / $len[$n] : 0;
	    
	    $ave_homfrac += $homfrac;
	    $std_homfrac += $homfrac * $homfrac;
	    
	    $ntot ++;
	}
    }
    FUNCS::calculate_averages(\$ave_dn, \$std_dn, $ntot);
    print "domain length $ave_dn +/- $std_dn\n"; 

    my $homfractot = ($lentot > 0)? $homtot / $lentot : 0;

    FUNCS::calculate_averages(\$ave_homfrac, \$std_homfrac, $ntot);
    printf "domain coverage %.2f (%.2f +/- %.2f)\n\n", 100 * $homfractot, 100*$ave_homfrac, 100*$std_homfrac; 
}

my $efficiency = 0;
if ($opt_E) { $efficiency = 1; }

my $viewplots = 1;
if ($opt_p) { $viewplots = 0; }

my $blast    = "NCBIBLAST";
my $muscle   = "MUSCLE";
my $msaprobs = "MSAProbs";
my $phmmer   = "PHMMER";
my $ephmmer  = "EPHMMER";

my $together = 0;
if ($opt_t) { $together = 1; }

# calculate efficiency
# N >= 2
# we assume that the last file is the one with the "optimal" scores
#
if ($efficiency && $N >1) {
    e2msabench_efficiency($N, \@file, \@filename, $name, $viewplots);
}

if (!$efficiency) {
    if (!$together) {
	for (my $n = 0; $n < $N; $n ++) {
	    e2msabench_onefile($file[$n], $minid, $maxid, $binsize, $shift, $viewplots);
	}
    }
    
# compare
    $viewplots = 1;
    if ($N > 0) {
	my @statfile;
	for (my $n = 0; $n < $N; $n ++) {
	    if    ($filename[$n] =~ /^msaprobs/) { $statfile[$n] = "$file[$n].msaprob.out.stat"; }
	    elsif ($filename[$n] =~ /^muscle/)   { $statfile[$n] = "$file[$n].muscle.out.stat"; }
	    elsif ($filename[$n] =~ /^blast/)    { $statfile[$n] = "$file[$n].blast.out.stat"; }
	    else                                 { $statfile[$n] = "$file[$n].out.stat"; }
	}
	e2msabench_together($N, \@statfile, \@filename, $name, $viewplots);
    }
}

sub e2msabench_efficiency{
    my ($Nf, $file_ref, $filename_ref, $name, $viewplots) = @_;
    
    my @effifile;
    
    #histogram
    my $N     = 100;
    my $k     = 1.0/$binsize;
    my $shift = 0;
    

    my $nfref;
    my @nameref;
    my @scref;
    my @timeref;
    benchfile_getsc($file_ref->[$Nf-1], \$nfref, \@nameref, \@scref, \@timeref);
    printf("nmsa %d\n", $nfref);

    for (my $n = 0; $n < $Nf-1; $n ++) {
	$effifile[$n] = "$file_ref->[$n].E.n$n";
	benchfile_getEff($effifile[$n], $N, $k, $shift, $file_ref->[$n], $nfref, \@nameref, \@scref, \@timeref);
    }
    
    plot_efficiency($Nf, \@effifile, "E");
    plot_efficiency($Nf, \@effifile, "scE");
    
    for (my $n = 0; $n <= $Nf-2; $n ++) { system("rm $effifile[$n]\n"); }
}

sub plot_efficiency
{

    my ($Nf, $effifile_ref, $which) = @_;

    my $psfile  = "$name.$which.ps";
    my $pdffile = "$name.$which.pdf";
    my $xlabel  = "\% ID of aligned region";
    my $ylabel;
    
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
    
    print GP "set output '$psfile'\n";
    print GP "set key right top\n";
      
    my $ymin;
    my $field;
    if ($which =~ /^E$/) {
	$ylabel  = "Score Efficiency (\%)";
	$field = 2;
	$ymin = 0.;
    }
    elsif ($which =~ /^scE$/) {
	$ylabel  = "Score Efficiency per position(\%)";
	$field = 4;
	$ymin = 70.;
    }
    my $field_std = $field + 1;
    
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    print GP "set yrange  [0:100]\n";
    print GP "set xrange  [100:0]\n";
    print GP "set key bottom\n";

    my $cmd = "";
    my $x = 1111;
    for (my $n = 0; $n < $Nf-2; $n ++) {
	$cmd  .= "'$effifile_ref->[$n]'   using 1:$field:$field_std with yerrorbars title '' ls $x, '$effifile_ref->[$n]'   using 1:$field with linespoints title '' ls $x, ";
	$x ++; if ($x > 1120) { $x = 1111; }
    }
    $cmd  .= "'$effifile_ref->[$Nf-2]'   using 1:$field:$field_std with yerrorbars title '' ls $x, '$effifile_ref->[$Nf-2]'   using 1:$field  with linespoints title '' ls $x";
    print GP "plot $cmd\n";
    close(GP);


    system("ps2pdf $psfile\n");
    system("rm $psfile\n");
    if ($viewplots) { system("open $pdffile&\n"); }

}

sub benchfile_getsc{
   my ($benchfile, $n_ref, $name_ref, $sc_ref, $time_ref) = @_;
    
   my $n = 0;
   open(BENCH,   "$benchfile")   || die;
   while (<BENCH>) {
       if    (/\#/) { next; }
       elsif (/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s*$/) {
	   my $name        = $1;
	   my $treeavgt    = $11;
	   my $method      = $13;
	   my $sc          = $24;
	   if ($method =~ /^$muscle$/ || $method =~ /^$msaprobs$/ || $method =~ /^$blast$/ || $method =~ /^$phmmer/) { next; }
	   else {
	       $name_ref->[$n] = $name;
	       $sc_ref->[$n]   = $sc;
	       $time_ref->[$n] = $treeavgt;
	       print "$n $name $sc $treeavgt\n";

	       $n ++;
	   }
       }
   }
   close(BENCH);

   $$n_ref = $n;
}

sub benchfile_getEff{
   my ($effifile, $N, $k, $shift, $benchfile, $nref, $name_ref, $scref_ref, $time_ref) = @_;
   
   my @hisID;
   my @ave_EffID;
   my @std_EffID;
   my @ave_EffscaledID;
   my @std_EffscaledID;

   my $nt = 0;
   my $avgEff   = 0;
   my $stdEff   = 0;
   my $avgEffsc = 0;
   my $stdEffsc = 0;

   my $nt50 = 0;
   my $avgEff50   = 0;
   my $stdEff50   = 0;

   my $nt60 = 0;
   my $avgEff60   = 0;
   my $stdEff60   = 0;

   FUNCS::init_histo_array($N, $k, \@hisID);
   FUNCS::init_histo_array($N, $k, \@ave_EffID);
   FUNCS::init_histo_array($N, $k, \@std_EffID);
   FUNCS::init_histo_array($N, $k, \@ave_EffscaledID);
   FUNCS::init_histo_array($N, $k, \@std_EffscaledID);
   
   my $idx = 0;
   open(BENCH,   "$benchfile")   || die;
   while (<BENCH>) {
       if    (/\#/) { next; }
       elsif (/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s*$/) {
	   my $name       = $1;
	   my $tph        = $2;
	   my $th         = $3;
	   my $fh         = $4;
	   my $nhe        = $5; # non homologous inferred
	   my $nhr        = $6; # non-homologous reference
	   my $rid        = $7; # the avgpid of the reference alignment 
	   my $eid        = $8;
	   my $rmatch     = $9;
	   my $ematch     = $10;
	   my $evodist    = $11;
	   my $time       = $12;
	   my $method     = $13;
	   my $ralen      = $14;
	   my $ealen      = $15;
	   my $avglen     = $16;
	   my $sc         = $24;

	    if ($tph == 0 && $fh == 0 && $nhe == 0 && $th > 0) { #this indicates the method did not detect this homology
		next;
	    } 

	   my $ID = $rid;
	   if ($opt_D) { $ID = pid_from_domain($idx, $name, $NS, \@name, \@d1pid, \@d2pid, \@d1n, \@d2n, \@len); }
	   $idx ++;
	   
	   if ($method =~ /^$muscle$/ || $method =~ /^$msaprobs$/ || $method =~ /^$blast$/ || $method =~ /^$phmmer/) { next; }
	   else {
	       my $eff;
	       my $effsc;
	       my $n;


	       for ($n = 0; $n < $nref; $n ++) {
		   if ($name =~ /^$name_ref->[$n]$/) {

		       my $optsc = $scref_ref->[$n];

		       if ($optsc > -1000000) {
			   $effsc = 100 * exp($sc/$avglen - $optsc/$avglen);
			   $eff   = ($sc >= $scref_ref->[$n])? 100. : 100 * exp($sc - $optsc);
			   printf("%s> %s nameref %s | %f %f | %f %f %f | avglen %f ralen %f | $tph $th $fh\n", $effifile, $name, $name_ref->[$n], $ID, $time_ref->[$n], $eff, $sc, $optsc, $avglen, $ralen);
			   
			   $nt ++;
			   $avgEff   += $eff;
			   $stdEff   += $eff*$eff;
			   $avgEffsc += $effsc;
			   $stdEffsc += $effsc*$effsc;
			   
			   if ($ID < 55.0 && $ID > 45.0) { $avgEff50 += $eff; $stdEff50 += $eff*$eff; $nt50 ++; }
			   if ($ID > 60.0)               { $avgEff60 += $eff; $stdEff60 += $eff*$eff; $nt60 ++; }
			   
			   FUNCS::fill_histo_array(1,             $ID, $N, $k, $shift, \@hisID);
			   FUNCS::fill_histo_array($eff,          $ID, $N, $k, $shift, \@ave_EffID);
			   FUNCS::fill_histo_array($eff*$eff,     $ID, $N, $k, $shift, \@std_EffID);    	    
			   FUNCS::fill_histo_array($effsc,        $ID, $N, $k, $shift, \@ave_EffscaledID);
			   FUNCS::fill_histo_array($effsc*$effsc, $ID, $N, $k, $shift, \@std_EffscaledID);   
		       } 	    
		       last;
		   }
	       }
	       if ($n == $nref) { last; }
	   }   
       }
   }
   close(BENCH);
   
   FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_EffID,       \@std_EffID,       , 0);
   FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_EffscaledID, \@std_EffscaledID, , 0);
   
   FUNCS::calculate_averages(\$avgEff, \$stdEff, $nt);
   printf("Eff         %f +\- %f\n", $avgEff, $stdEff);
   FUNCS::calculate_averages(\$avgEff50, \$stdEff50, $nt50);
   printf("Eff 45-55   %f +\- %f\n", $avgEff50, $stdEff50);
   FUNCS::calculate_averages(\$avgEff60, \$stdEff60, $nt60);
   printf("Eff >60     %f +\- %f\n", $avgEff60, $stdEff60);

   FUNCS::calculate_averages(\$avgEffsc, \$stdEffsc, $nt);
   printf("Eff  scaled %f +\- %f\n", $avgEffsc, $stdEffsc);

   open(EFILE, ">$effifile") || die;
   my $dim = $N * $k;
   
   for (my $i=0; $i<=$dim; $i++) { 
       my $idl = $i/$k     - $shift;
       my $idh = (($i+1)/$k < $N)? ($i+1)/$k - $shift : $N-$shift;
       
      printf EFILE "%f %f %f %f %f\n", $idl+0.5*($idh-$idl), $ave_EffID[$i], $std_EffID[$i], $ave_EffscaledID[$i], $std_EffscaledID[$i]; 
   }
   close(EFILE);
}

sub e2msabench_together{
    my ($N, $statfile_ref, $filename_ref, $name, $viewplots) = @_;

    my $which     = "totalF";
    FUNCS::plot_id2totalF($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "totalSEN";
    FUNCS::plot_id2totalF($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "totalPPV";
    FUNCS::plot_id2totalF($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "totalSPE";
    FUNCS::plot_id2totalF($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);

    $which     = "totalCF";
    FUNCS::plot_id2totalF($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "totalCSEN";
    FUNCS::plot_id2totalF($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "totalCPPV";
    FUNCS::plot_id2totalF($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);

    $viewplots = 0;
    $which     = "F";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "SEN";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "PPV";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "SPE";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
  
    $which     = "CF";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "CSEN";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "CPPV";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
   
    if (0) {
    $which     = "SCL";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "MATCH";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "OPENG";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    }
}

sub e2msabench_onefile{
    my ($benchfile, $minid, $maxid, $binsize, $shift, $viewplots) = @_;

    my $outf    = "$benchfile.out";
    my $outf_pr = "$benchfile.msaprob.out";
    my $outf_mu = "$benchfile.muscle.out";
    my $outf_bl = "$benchfile.blast.out";
    
    parse_benchfile($benchfile, $outf, $outf_pr, $outf_mu, $outf_bl, $minid, $maxid);
    
    my $xlabel = "Tree mean distance";
    plot_id2divergence($outf, 3, $xlabel, $viewplots);
     
    plot_id2bench_histo($outf,    $binsize, $shift, $viewplots);
    
    $viewplots = 0;
    plot_id2bench_histo($outf_pr, $binsize, $shift, $viewplots);
    plot_id2bench_histo($outf_mu, $binsize, $shift, $viewplots);
    plot_id2bench_histo($outf_bl, $binsize, $shift, $viewplots);
    
}

sub final_stats_file {
    my ($outfile, $N, $k, $shift, 
	$hisID_ref, $hisMATCH_ref,
	$tph_histo_ref, $th_histo_ref, $fh_histo_ref, 
	$nhe_histo_ref, $nhr_histo_ref,  
	$hr_histo_ref, $he_histo_ref, $hre_histo_ref, 
	$ave_senhisID_ref,   $std_senhisID_ref, 
	$ave_ppvhisID_ref,   $std_ppvhisID_ref, 
	$ave_FhisID_ref,     $std_FhisID_ref, 
	$ave_SPEhisID_ref,   $std_SPEhisID_ref, 
	$ave_CsenhisID_ref,  $std_CsenhisID_ref, 
	$ave_CppvhisID_ref,  $std_CppvhisID_ref, 
	$ave_CFhisID_ref,    $std_CFhisID_ref, 
	$ave_sclhisID_ref,   $std_sclhisID_ref, 
	$ave_matchhisID_ref, $std_matchhisID_ref, 
	$ave_openghisID_ref, $std_openghisID_ref) = @_;
    
    open(OUTF, ">$outfile") || die;
    printf "$outfile\n";

    my $dim = $N * $k;
    
    my $AUCsen  = 0.0;
    my $AUCppv  = 0.0;
    my $AUCF    = 0.0;
    my $AUCSPE  = 0.0;

    my $AUCCsen = 0.0;
    my $AUCCppv = 0.0;
    my $AUCCF   = 0.0;
 
    my $AUCavgsen = 0.0;
    my $AUCavgppv = 0.0;
    my $AUCavgF   = 0.0;
    my $AUCavgSPE = 0.0;

    my $AUCCavgsen = 0.0;
    my $AUCCavgppv = 0.0;
    my $AUCCavgF   = 0.0;

    my $num = 0;

    printf("\n");
    for (my $i=0; $i<=$dim; $i++) { 
	my $idl = $i/$k     - $shift;
	my $idh = (($i+1)/$k < $N)? ($i+1)/$k - $shift : $N-$shift;
	
	my $sen;
	my $ppv;
	my $F;
	FUNCS::calculateF($tph_histo_ref->[$i], $th_histo_ref->[$i], $fh_histo_ref->[$i], \$sen, \$ppv, \$F);
	
	my $SPE = ($nhr_histo_ref->[$i])? 100.*$nhe_histo_ref->[$i]/$nhr_histo_ref->[$i]: 0.0;
	
	my $Csen;
	my $Cppv;
	my $CF;
	FUNCS::calculateF($hre_histo_ref->[$i], $hr_histo_ref->[$i], $he_histo_ref->[$i], \$Csen, \$Cppv, \$CF);

	my $nmsa = $hisID_ref->[$i];
	
	if ($nmsa >= $min_nmsa) {
	    $num ++;
	    $AUCsen  += $sen;
	    $AUCppv  += $ppv;
	    $AUCF    += $F;
	    $AUCSPE  += $SPE;
	    
	    $AUCCsen += $Csen;
	    $AUCCppv += $Cppv;
	    $AUCCF   += $CF;

	    $AUCavgsen  += $ave_senhisID_ref->[$i];
	    $AUCavgppv  += $ave_ppvhisID_ref->[$i];
	    $AUCavgF    += $ave_FhisID_ref->[$i];
	    $AUCavgSPE  += $ave_SPEhisID_ref->[$i];

	    $AUCCavgsen += $ave_CsenhisID_ref->[$i];
	    $AUCCavgppv += $ave_CppvhisID_ref->[$i];
	    $AUCCavgF   += $ave_CFhisID_ref->[$i];
	}

	if (0) {
	    printf("nmsa %10d nmsaMATCH %10d id [%3.2f,%3.2f) ", $nmsa, $hisMATCH_ref->[$i], $idl, $idh);
	    printf("sen %3.2f [%3.2f+/-%3.2f]\tppv %3.2f [%3.2f+/-%3.2f]\tF %3.2f [%3.2f+/-%3.2f] \tSPE %3.2f [%3.2f+/-%3.2f]\tCsen %3.2f [%3.2f+/-%3.2f]\tCppv %3.2f [%3.2f+/-%3.2f]\tCF %3.2f [%3.2f+/-%3.2f]  [%3.2f+/-%3.2f] [%3.2f+/-%3.2f] [%3.2f+/-%3.2f] \n", 
		   $sen, $ave_senhisID_ref->[$i], $std_senhisID_ref->[$i], 
		   $ppv, $ave_ppvhisID_ref->[$i], $std_ppvhisID_ref->[$i], 
		   $F,   $ave_FhisID_ref->[$i],   $std_FhisID_ref->[$i],
		   $SPE, $ave_SPEhisID_ref->[$i], $std_SPEhisID_ref->[$i],
		   $Csen, $ave_CsenhisID_ref->[$i], $std_CsenhisID_ref->[$i], 
		   $Cppv, $ave_CppvhisID_ref->[$i], $std_CppvhisID_ref->[$i], 
		   $CF,   $ave_CFhisID_ref->[$i],   $std_CFhisID_ref->[$i],
		   $ave_sclhisID_ref->[$i],   $std_sclhisID_ref->[$i],
		   $ave_matchhisID_ref->[$i], $std_matchhisID_ref->[$i],
		   $ave_openghisID_ref->[$i], $std_openghisID_ref->[$i]);
	}
	printf OUTF "%d %d %f %f ", $nmsa, $hisMATCH_ref->[$i], $idl, $idl + 0.5*($idh-$idl); 
	if ($nmsa >= $min_nmsa) {
	    printf OUTF "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
	    $sen,  $ave_senhisID_ref->[$i],  $std_senhisID_ref->[$i], 
	    $ppv,  $ave_ppvhisID_ref->[$i],  $std_ppvhisID_ref->[$i], 
	    $F,    $ave_FhisID_ref->[$i],    $std_FhisID_ref->[$i],
	    $SPE,  $ave_SPEhisID_ref->[$i],  $std_SPEhisID_ref->[$i],
	    $Csen, $ave_CsenhisID_ref->[$i], $std_CsenhisID_ref->[$i], 
	    $Cppv, $ave_CppvhisID_ref->[$i], $std_CppvhisID_ref->[$i], 
	    $CF,   $ave_CFhisID_ref->[$i],   $std_CFhisID_ref->[$i],
	    $ave_sclhisID_ref->[$i],   $std_sclhisID_ref->[$i],
	    $ave_matchhisID_ref->[$i], $std_matchhisID_ref->[$i],
	    $ave_openghisID_ref->[$i], $std_openghisID_ref->[$i]; 
	}
	printf OUTF "\n"; 
	
    }

    my $denom = $num;
    $denom = ($denom>0)? 1.0/$denom : 1.0;

    $AUCsen  *= $denom;
    $AUCppv  *= $denom;
    $AUCF    *= $denom;
    $AUCSPE  *= $denom;

    $AUCCsen *= $denom;
    $AUCCppv *= $denom;
    $AUCCF   *= $denom;
    
    $AUCavgsen  *= $denom;
    $AUCavgppv  *= $denom;
    $AUCavgF    *= $denom;
    $AUCavgSPE  *= $denom;
    
    $AUCCavgsen *= $denom;
    $AUCCavgppv *= $denom;
    $AUCCavgF   *= $denom;
    
    printf OUTF "#AUC %.1f %.1f %.1f %.1f | %.1f %.1f %.1f %.1f\n", $AUCsen, $AUCppv, $AUCF, $AUCSPE, $AUCavgsen, $AUCavgppv, $AUCavgF, $AUCavgSPE;  
    printf      "#AUC %.1f %.1f %.1f %.1f | %.1f %.1f %.1f %.1f\n", $AUCsen, $AUCppv, $AUCF, $AUCSPE, $AUCavgsen, $AUCavgppv, $AUCavgF, $AUCavgSPE; 
    
    printf OUTF "#C_AUC %.1f %.1f %.1f  | %.1f %.1f %.1f \n", $AUCCsen, $AUCCppv, $AUCCF, $AUCCavgsen, $AUCCavgppv, $AUCCavgF;  
    printf      "#C_AUC %.1f %.1f %.1f  | %.1f %.1f %.1f \n", $AUCCsen, $AUCCppv, $AUCCF, $AUCCavgsen, $AUCCavgppv, $AUCCavgF; 

    close (OUTF);
    
}


sub parse_benchfile{
    my ($benchfile, $outfile, $outfile_pr, $outfile_mu, $outfile_bl, $minid, $maxid) = @_;
    
    my $ttph = 0;
    my $tth  = 0;
    my $tfh  = 0;
    my $tnhr = 0;
    my $tnhe = 0;
    my $ttph_bl = 0;
    my $tth_bl  = 0;
    my $tfh_bl  = 0;
    my $tnhe_bl = 0;
    my $tnhr_bl = 0;
    my $ttph_mu = 0;
    my $tth_mu  = 0;
    my $tfh_mu  = 0;
    my $tnhe_mu = 0;
    my $tnhr_mu = 0;
    my $ttph_pr = 0;
    my $tth_pr  = 0;
    my $tfh_pr  = 0;
    my $tnhe_pr = 0;
    my $tnhr_pr = 0;
     
    my $ttime    = 0.;
    my $ttime_bl = 0.;
    my $ttime_mu = 0.;
    my $ttime_pr = 0.;

    my $nmsa       = 0;
    my $nmsa_bl    = 0;
    my $nmsa_mu    = 0;
    my $nmsa_pr    = 0;
    my $nmsatot    = 0;
    my $nmsatot_bl = 0;
    my $nmsatot_mu = 0;
    my $nmsatot_pr = 0;

    my $method;
    my $umethod = "";
    my $ignore;

    my $nmsa_detect = 0;

    my $idxf = 0;
    my $idxs = 0;
    open(OUTF,    ">$outfile")    || die;
    open(OUTF_PR, ">$outfile_pr") || die;
    open(OUTF_MU, ">$outfile_mu") || die;
    open(OUTF_BL, ">$outfile_bl") || die;
    open(BENCH,   "$benchfile")   || die;
    while (<BENCH>) {
	if    (/\#/) { next; }
	elsif (/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s*(.*)$/) {
	    my $name       = $1;
	    my $tph        = $2;
	    my $th         = $3;
	    my $fh         = $4;
	    my $nhe        = $5; # non homologous inferred
	    my $nhr        = $6; # non-homologous reference

	    my $rid        = $7; # the avgpid of the reference alignment 
	    my $eid        = $8;
	    my $rmatch     = $9;
	    my $ematch     = $10;
	    my $evodist    = $11;
	    my $time       = $12;
	    $method        = $13;

	    if ($tph == 0 && $fh == 0 && $nhe == 0 && $th > 0) { #this indicates the method did not detect this homology
	    } 
	    else { $nmsa_detect ++; }

	    #instead of using the pid of the ref aligment, use that of the actual homology domains
	    #given in the summary file
	    my $ID = $rid;
	    if ($opt_D) { $ID = pid_from_domain($idxf, $name, $NS, \@name, \@d1pid, \@d2pid, \@d1n, \@d2n, \@len); }
	    #printf("ID %f rid %f\n", $ID, $rid);

	    $idxf ++;
	    my $r_alen     = $14;
	    my $e_alen     = $15;

	    my $r_avgsqlen = $16;
	    my $r_stdsqlen = $17;
	    my $e_avgsqlen = $18;
	    my $e_stdsqlen = $19;

	    my $r_avginum  = $20;
	    my $e_avginum  = $21;
	    my $r_avgilen  = $22;
	    my $e_avgilen  = $23;
	    my $sc         = $24;
	    my $rest       = $25;

	    # for outputs that include also coverage measures
	    my $hr         = 0;
	    my $he         = 0;
	    my $hre        = 0;

	    if ($rest =~ /^\s*(\d+)\s+(\d+)\s+(\d+)\s*$/) {
		$hr         = $1;
		$he         = $2;
		$hre        = $3;
	    }

	    my $s;
	    my $p;
	    my $F;
	    FUNCS::calculateF($tph, $th, $fh,  \$s,  \$p,  \$F);
	    my $Cs;
	    my $Cp;
	    my $CF;
	    FUNCS::calculateF($hre, $hr, $he,  \$Cs,  \$Cp,  \$CF);
	    #printf("sen %f ppv %f F %f | Csen %f Cppv %f CF %f\n", $s, $p, $F, $Cs, $Cp, $CF);

	    my $SPE = ($nhr > 0)? 100.0*$nhe/$nhr : 0.0;

	    if ($method =~ /^1$muscle$/) {
		$nmsatot_mu ++;
		if ($ID < $minid || $ID > $maxid) { next; }
		$nmsa_mu  ++;
		$ttph_mu  += $tph;
		$tth_mu   += $th;
		$tfh_mu   += $fh;
		$tnhe_mu  += $nhe;
		$tnhr_mu  += $nhr;
		$ttime_mu += $time; 

		printf OUTF_MU "%f %f %f %d %d %d %d %d %f %f %f %f %d %d %d %f %f %f \n", $ID, $eid, $evodist, $tph, $th, $fh, $nhe, $nhr, $s, $p, $F, $SPE, $hr, $he, $hre, $Cs, $Cp, $CF;
	    }
	    elsif ($method =~ /^1$msaprobs$/) {
		$nmsatot_pr ++;
		if ($ID < $minid || $ID > $maxid) { next; }
		$nmsa_pr  ++;
		$ttph_pr  += $tph;
		$tth_pr   += $th;
		$tfh_pr   += $fh;
		$tnhe_pr  += $nhe;
		$tnhr_pr  += $nhr;
		$ttime_pr += $time;
		#printf "MS %d %d %d | %d %d %d | %f %f | %f %f \n", $tph, $th, $fh, $ttph_pr, $tth_pr, $tfh_pr, $tph/$th, $tph/$fh, $ttph_pr/$tth_pr, $ttph_pr/$tfh_pr;

		printf OUTF_PR "%f %f %f %d %d %d %d %d %f %f %f %f %d %d %d %f %f %f \n", $ID, $eid, $evodist, $tph, $th, $fh, $nhe, $nhr, $s, $p, $F, $SPE, $hr, $he, $hre, $Cs, $Cp, $CF;
	    }
	    else {
		$nmsatot ++;
		if ($ID < $minid || $ID > $maxid) { next; }
		$nmsa    ++;
		$ttph    += $tph;
		$tth     += $th;
		$tfh     += $fh;
		$tnhe    += $nhe;
		$tnhr    += $nhr;
		$ttime   += $time; 
		$umethod  = $method;

		printf OUTF "%f %f %f %d %d %d %d %d %f %f %f %f %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", 
		$ID, $eid, $evodist, $tph, $th, $fh, $nhe, $nhr, $s, $p, $F, $SPE, $hr, $he, $hre, $Cs, $Cp, $CF,
		$rmatch, $ematch, 
		$r_alen, $e_alen, 
		$r_avgsqlen, $r_stdsqlen, $e_avgsqlen, $e_stdsqlen, 
		$r_avginum, $e_avginum, $r_avgilen, $e_avgilen,
		$sc;
	    }
	}
    }
    close (BENCH);
    close (OUTF);
    close (OUTF_PR);
    close (OUTF_MU);
    close (OUTF_BL);
    
    my $sen    = 0.;
    my $ppv    = 0.;
    my $F      = 0.;
    my $SPE    = 0.;

    my $sen_bl = 0.;
    my $ppv_bl = 0.;
    my $F_bl   = 0.;
    my $SPE_bl = 0.;

    my $sen_mu = 0.;
    my $ppv_mu = 0.;
    my $F_mu   = 0.;
    my $SPE_mu = 0.;

    my $sen_pr = 0.;
    my $ppv_pr = 0.;
    my $F_pr   = 0.;
    my $SPE_pr = 0.;
    
    FUNCS::calculateF($ttph,    $tth,    $tfh,    \$sen,    \$ppv,    \$F);
    FUNCS::calculateF($ttph_bl, $tth_bl, $tfh_bl, \$sen_bl, \$ppv_bl, \$F_bl);
    FUNCS::calculateF($ttph_mu, $tth_mu, $tfh_mu, \$sen_mu, \$ppv_mu, \$F_mu);
    FUNCS::calculateF($ttph_pr, $tth_pr, $tfh_pr, \$sen_pr, \$ppv_pr, \$F_pr);
    
    $SPE    = ($tnhr    > 0.)? 100.*$tnhe/$tnhr       : 0.0;
    $SPE_bl = ($tnhr_bl > 0.)? 100.*$tnhe_bl/$tnhr_bl : 0.0;
    $SPE_mu = ($tnhr_mu > 0.)? 100.*$tnhe_mu/$tnhr_mu : 0.0;
    $SPE_pr = ($tnhr_pr > 0.)? 100.*$tnhe_pr/$tnhr_pr : 0.0;

    printf("\nNMSA %d/%d nmsa_detected %d ID[%f,%f] total_homologies %d \n", $nmsa, $nmsatot, $nmsa_detect, $minid, $maxid, $tth);
    if ($nmsatot_bl != $nmsatot) { printf("nmsa_BLAST %d/%d nmsa_detected %d ID[%f,%f]  total_homologies %d\n",       $nmsa_bl, $nmsatot_bl, $nmsa_detect, $minid, $maxid, $tth_bl); }
    if ($nmsatot_mu != $nmsatot) { printf("nmsa_MUSCLE %d/%d nmsa_detected %d ID[%f,%f]  total_homologies %d\n",      $nmsa_mu, $nmsatot_mu, $nmsa_detect, $minid, $maxid, $tth_mu); }
    if ($nmsatot_pr != $nmsatot) { printf("nmsa_MSAProbs %d/%d nmsa_detected %d ID[%f,%f]  total_homologies %d\n",    $nmsa_pr, $nmsatot_pr, $nmsa_detect, $minid, $maxid, $tth_pr); }
    printf("sen %.2f\tppv %.2f\tF %.2f\t\tSPE %.2f\t\t%s\tCPU Time (s) %.f\n",   $sen,    $ppv,    $F,    $SPE,    $umethod,  $ttime);
    printf("sen %.2f\tppv %.2f\tF %.2f\t\tSPE %.2f\t\t%s\t\tCPU Time (s) %.f\n", $sen_bl, $ppv_bl, $F_bl, $SPE_bl, $blast,    $ttime_bl);
    printf("sen %.2f\tppv %.2f\tF %.2f\t\tSPE %.2f\t\t%s\t\tCPU Time (s) %.f\n", $sen_mu, $ppv_mu, $F_mu, $SPE_mu, $muscle,   $ttime_mu);
    printf("sen %.2f\tppv %.2f\tF %.2f\t\tSPE %.2f\t\t%s\tCPU Time (s) %.f\n",   $sen_pr, $ppv_pr, $F_pr, $SPE_pr, $msaprobs, $ttime_pr);
    
}


sub parse_summary {

   my ($summaryfile, $ret_n, $name_ref, $d1pid_ref, $d2pid_ref, $d1n_ref, $d2n_ref, $len_ref) = @_;

   my $n = 0;
   my $ndom  = -1;
   my $d1pid = -1;
   my $d2pid = -1;
   my $d1n   = -1;
   my $d2n   = -1;
   my $len   = -1;

    open(SUMM,   "$summaryfile")   || die;
    while (<SUMM>) {
	if    (/\#/) { next; }
	elsif (/^(\S+)\s+(\d+)\s+\S+\s+\d+\s+(\d+)\s+(\S+)\s+\d+(.+)$/) {
	    $name_ref->[$n]  = $1;
	    $len_ref->[$n]   = $2;
	    $d1n_ref->[$n]   = $3;
	    $d1pid_ref->[$n] = $4;
	    my $rest         = $5;

	    $d2n_ref->[$n]   = -1.;
	    $d2pid_ref->[$n] = -1.;
	    $ndom = 1;

	    if ($rest =~ /^\s*(\d+)\s+(\S+)\s+$/) {
		$d2n_ref->[$n]   = $1;
		$d2pid_ref->[$n] = $2;
		$ndom = 2;
	    }
	    $n ++;
	}
    }
   close(SUMM);

   $$ret_n = $n;

   return $ndom;
}

sub pid_from_domain {

    my ($idx, $name, $NS, $name_ref, $d1pid_ref, $d2pid_ref, $d1n_ref, $d2n_ref, $len_ref) = @_;

    my $ID = -1;
    my $x = $idx;

    my $len = -1;
    my $d1n = -1;
    my $d2n = -1;
    my $d1pid = -1;
    my $d2pid = -1;

    #printf "name $name || %s\n", $name_ref->[$x];
    if ($name =~ /^$name_ref->[$x]$/) {
	$len   = $len_ref->[$x];
	$d1n   = $d1n_ref->[$x];
	$d2n   = $d2n_ref->[$x];
	$d1pid = $d1pid_ref->[$x];
	$d2pid = $d2pid_ref->[$x];
	#printf "   ^^x $x name $name d1pid $d1pid d2pid $d2pid d1n $d1n d2n $d2n len $len\n";
    }
    else {
	while ($x < $NS-1) {
	    $x ++;
	    #printf " -->name $name %s\n", $name_ref->[$x];
	    if ($name =~ /^$name_ref->[$x]$/) {
		$len   = $len_ref->[$x];
		$d1n   = $d1n_ref->[$x];
		$d2n   = $d2n_ref->[$x];
		$d1pid = $d1pid_ref->[$x];
		$d2pid = $d2pid_ref->[$x];
		#printf "      ++++x $x name $name d1pid $d1pid d2pid $d2pid\n";
		last;
	    }
	}
	
	if ($x == $NS) { find_dompid_in_summary($name, $NS, $name_ref, $d1pid_ref, $d2pid_ref, $d1n_ref, $d2n_ref, $len_ref, \$d1pid, \$d2pid, \$d1n, \$d2n, \$len); }
    }
    
    if ($d2pid >= 0) { 
	$ID  = ($d1pid + $d2pid) * 0.5; 
    }
    else { 
	$ID  = $d1pid;
    }
    
    return $ID;
}

sub find_dompid_in_summary {
    my ($name, $n, $name_ref, $d1pid_ref, $d2pid_ref, $d1n_ref, $d2n_ref, $len_ref, $ret_d1pid, $ret_d2pid, $ret_d1n, $ret_d2n, $ret_len) = @_;

    my $len   = -1;
    my $d1n   = -1;
    my $d2n   = -1;
    my $d1pid = -1;
    my $d2pid = -1;

    my $found = 0;
    for (my $i = 0; $i < $n; $i ++) {
	if ($name =~ /^$name_ref->[$i]$/) {
	    $len   = $len_ref->[$i];
	    $d1n   = $d1n_ref->[$i];
	    $d2n   = $d2n_ref->[$i];
	    $d1pid = $d1pid_ref->[$i];
	    $d2pid = $d2pid_ref->[$i];
	    $found ++;
	}
    }
    if ($found > 1) { print "oops $name found more than once in summary file\n"; die; }

    $$ret_len   = $len;
    $$ret_d1n   = $d1n;
    $$ret_d2n   = $d2n;
    $$ret_d1pid = $d1pid;
    $$ret_d2pid = $d2pid;
}

sub plot_id2divergence {

    my ($file, $field, $xlabel, $viewplots) = @_;

    my $psfile  = "$file.ps";
    my $pdffile = "$file.pdf";
    my $ylabel  = "\% ID of aligned region";
 
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    print GP "set output '$psfile'\n";
    print GP "set key right top\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    #print GP "set yrange  [0:100]\n";
    print GP "set xrange [100:0]\n";
    #print GP "set logscale x\n";
    #print GP "set nokey\n";
    my $cmd = "";
    $cmd  = "'$file' using $field:1 title 'reference msa' ls 1, '$file' using $field:2 title 'infered msa' ls 2";
    print GP "plot $cmd\n";
  
    close(GP);

    system("ps2pdf $psfile\n");
    system("rm $psfile\n");
    if ($viewplots) { system("open $pdffile&\n"); }
}


sub plot_id2bench_histo {

    my ($outfile, $binsize, $shift, $viewplots) = @_;

    my $statfile = "$outfile.stat";

    my $Sfileps   = "$statfile.SEN.ps";
    my $Pfileps   = "$statfile.PPV.ps";
    my $Ffileps   = "$statfile.F.ps";
    my $SPEfileps = "$statfile.SPE.ps";

    my $SCfileps = "$statfile.sc.ps";
    my $Mfileps  = "$statfile.M.ps";
    my $Ifileps  = "$statfile.I.ps";

    my $hidfileps    = "$statfile.histoID.ps";
    my $hmatchfileps = "$statfile.histoMatch.ps";

    my $xlabel = "\% ID of aligned region";
    my $ylabel = "";
    my $title  = "";
    my $key    = "";

    #histogram
    my $N     = 100;
    my $k     = 1.0/$binsize;

    my @hisID;
    my @hisMATCH;

    my @tph_histo;
    my @th_histo;
    my @fh_histo;
    my @nhe_histo;
    my @nhr_histo;
    my @hr_histo;
    my @he_histo;
    my @hre_histo;

    my @ave_senhisID;
    my @std_senhisID;
    my @ave_ppvhisID;
    my @std_ppvhisID;
    my @ave_FhisID;
    my @std_FhisID;
    my @ave_SPEhisID;
    my @std_SPEhisID;
    my @ave_CsenhisID;
    my @std_CsenhisID;
    my @ave_CppvhisID;
    my @std_CppvhisID;
    my @ave_CFhisID;
    my @std_CFhisID;

    FUNCS::init_histo_array($N, $k, \@hisID);
    FUNCS::init_histo_array($N, $k, \@hisMATCH);
    FUNCS::init_histo_array($N, $k, \@tph_histo);
    FUNCS::init_histo_array($N, $k, \@th_histo);
    FUNCS::init_histo_array($N, $k, \@fh_histo);
    FUNCS::init_histo_array($N, $k, \@nhe_histo);
    FUNCS::init_histo_array($N, $k, \@nhr_histo);
    FUNCS::init_histo_array($N, $k, \@hr_histo);
    FUNCS::init_histo_array($N, $k, \@he_histo);
    FUNCS::init_histo_array($N, $k, \@hre_histo);
    FUNCS::init_histo_array($N, $k, \@ave_senhisID);
    FUNCS::init_histo_array($N, $k, \@std_senhisID);
    FUNCS::init_histo_array($N, $k, \@ave_ppvhisID);
    FUNCS::init_histo_array($N, $k, \@std_ppvhisID);
    FUNCS::init_histo_array($N, $k, \@ave_FhisID);
    FUNCS::init_histo_array($N, $k, \@std_FhisID);
    FUNCS::init_histo_array($N, $k, \@ave_SPEhisID);
    FUNCS::init_histo_array($N, $k, \@std_SPEhisID);
    FUNCS::init_histo_array($N, $k, \@ave_CsenhisID);
    FUNCS::init_histo_array($N, $k, \@std_CsenhisID);
    FUNCS::init_histo_array($N, $k, \@ave_CppvhisID);
    FUNCS::init_histo_array($N, $k, \@std_CppvhisID);
    FUNCS::init_histo_array($N, $k, \@ave_CFhisID);
    FUNCS::init_histo_array($N, $k, \@std_CFhisID);

    my $longversion = 0;
    my @ave_sclhisID;
    my @std_sclhisID;
    my @ave_matchhisID;
    my @std_matchhisID;
    my @ave_openghisID;
    my @std_openghisID;
    FUNCS::init_histo_array($N, $k, \@ave_sclhisID);
    FUNCS::init_histo_array($N, $k, \@std_sclhisID);
    FUNCS::init_histo_array($N, $k, \@ave_matchhisID);
    FUNCS::init_histo_array($N, $k, \@std_matchhisID);
    FUNCS::init_histo_array($N, $k, \@ave_openghisID);
    FUNCS::init_histo_array($N, $k, \@std_openghisID);

    open(FILE, "$outfile")   || die;
    while (<FILE>) {
	if    (/\#/) { next; }
	elsif (/^(\S+)\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s*(\S.+)$/) {
	    my $rid  = $1;
	    my $tph  = $2;
	    my $th   = $3;
	    my $fh   = $4;
	    my $nhe  = $5;
	    my $nhr  = $6;
	    my $sen  = $7;
	    my $ppv  = $8;
	    my $F    = $9;
	    my $SPE  = $10;
	    my $hr   = $11;
	    my $he   = $12;
	    my $hre  = $13;
	    my $Csen = $14;
	    my $Cppv = $15;
	    my $CF   = $16;
	    #print("sen $sen ppv $ppv F $F | Csen $Csen Cppv $Cppv CF $CF\n");

	    my $rest = $17;

	    FUNCS::fill_histo_array(1,           $rid, $N, $k, $shift, \@hisID);
	    FUNCS::fill_histo_array($tph,        $rid, $N, $k, $shift, \@tph_histo);
	    FUNCS::fill_histo_array($th,         $rid, $N, $k, $shift, \@th_histo);
	    FUNCS::fill_histo_array($fh,         $rid, $N, $k, $shift, \@fh_histo);
	    FUNCS::fill_histo_array($nhe,        $rid, $N, $k, $shift, \@nhe_histo);
	    FUNCS::fill_histo_array($nhr,        $rid, $N, $k, $shift, \@nhr_histo);
	    FUNCS::fill_histo_array($hr,         $rid, $N, $k, $shift, \@hr_histo);
	    FUNCS::fill_histo_array($he,         $rid, $N, $k, $shift, \@he_histo);
	    FUNCS::fill_histo_array($hre,        $rid, $N, $k, $shift, \@hre_histo);

	    FUNCS::fill_histo_array($sen,        $rid, $N, $k, $shift, \@ave_senhisID);
	    FUNCS::fill_histo_array($sen*$sen,   $rid, $N, $k, $shift, \@std_senhisID);    	    
	    FUNCS::fill_histo_array($ppv,        $rid, $N, $k, $shift, \@ave_ppvhisID);
	    FUNCS::fill_histo_array($ppv*$ppv,   $rid, $N, $k, $shift, \@std_ppvhisID);    	    
	    FUNCS::fill_histo_array($F,          $rid, $N, $k, $shift, \@ave_FhisID);
	    FUNCS::fill_histo_array($F*$F,       $rid, $N, $k, $shift, \@std_FhisID);  
	    FUNCS::fill_histo_array($SPE,        $rid, $N, $k, $shift, \@ave_SPEhisID);
	    FUNCS::fill_histo_array($SPE*$SPE,   $rid, $N, $k, $shift, \@std_SPEhisID);  

	    FUNCS::fill_histo_array($Csen,       $rid, $N, $k, $shift, \@ave_CsenhisID);
	    FUNCS::fill_histo_array($Csen*$Csen, $rid, $N, $k, $shift, \@std_CsenhisID);    	    
	    FUNCS::fill_histo_array($Cppv,       $rid, $N, $k, $shift, \@ave_CppvhisID);
	    FUNCS::fill_histo_array($Cppv*$Cppv, $rid, $N, $k, $shift, \@std_CppvhisID);    	    
	    FUNCS::fill_histo_array($CF,         $rid, $N, $k, $shift, \@ave_CFhisID);
	    FUNCS::fill_histo_array($CF*$CF,     $rid, $N, $k, $shift, \@std_CFhisID);  

	    if ($rest =~ /^(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$/) {
		
		my $rmatch   = $1;
		my $ralen    = $2;
		my $ealen    = $3;
		my $rinum    = $4;
		my $einum    = $5;
		my $rilen    = $6;
		my $eilen    = $7;
		my $sc       = $8;
		$longversion = 0;
      
		my $scl   = ($ealen > 0)? $sc/$ealen : 0.0;
		my $match = $rmatch;
		my $openg = ($rmatch > 0. && $ralen > 0)? 10000. * $rinum / ($rmatch*$ralen) : 0.0;

		FUNCS::fill_histo_array(1,          $rmatch, $N, $k, $shift, \@hisMATCH);

		FUNCS::fill_histo_array($scl,          $rid, $N, $k, $shift, \@ave_sclhisID);
		FUNCS::fill_histo_array($scl*$scl,     $rid, $N, $k, $shift, \@std_sclhisID);
		FUNCS::fill_histo_array($match,        $rid, $N, $k, $shift, \@ave_matchhisID);
		FUNCS::fill_histo_array($match*$match, $rid, $N, $k, $shift, \@std_matchhisID);
		FUNCS::fill_histo_array($openg,        $rid, $N, $k, $shift, \@ave_openghisID);
		FUNCS::fill_histo_array($openg*$openg, $rid, $N, $k, $shift, \@std_openghisID);
	    }
 	}
    }
    close(FILE);
    my $medianID    = FUNCS::histogram_median($N, $k, \@hisID);
    my $medianMATCH = FUNCS::histogram_median($N, $k, \@hisMATCH);

    FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_senhisID,  \@std_senhisID,  , 0);
    FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_ppvhisID,  \@std_ppvhisID,  , 0);
    FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_FhisID,    \@std_FhisID,    , 0);
    FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_SPEhisID,  \@std_SPEhisID,  , 0);
    FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_CsenhisID, \@std_CsenhisID, , 0);
    FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_CppvhisID, \@std_CppvhisID, , 0);
    FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_CFhisID,   \@std_CFhisID,   , 0);
    if ($longversion) {
	FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_sclhisID,   \@std_sclhisID,   , 0);
	FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_matchhisID, \@std_matchhisID, , 0);
	FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_openghisID, \@std_openghisID, , 0);
    }
    
    final_stats_file($statfile, $N, $k, $shift, 
		     \@hisID, 
		     \@hisMATCH,
		     \@tph_histo, \@th_histo, \@fh_histo, 
		     \@nhe_histo, \@nhr_histo, 
		     \@hr_histo,  \@he_histo, \@hre_histo, 
		     \@ave_senhisID,    \@std_senhisID, 
		     \@ave_ppvhisID,    \@std_ppvhisID, 
		     \@ave_FhisID,      \@std_FhisID,  
		     \@ave_SPEhisID,    \@std_SPEhisID, 
		     \@ave_CsenhisID,   \@std_CsenhisID, 
		     \@ave_CppvhisID,   \@std_CppvhisID, 
		     \@ave_CFhisID,     \@std_CFhisID,  
		     \@ave_sclhisID,    \@std_sclhisID, 
		     \@ave_matchhisID,  \@std_matchhisID, 
		     \@ave_openghisID,  \@std_openghisID);
    
    printf("median_ID = $medianID\n");
    printf("median_MATCH = $medianMATCH\n");

    my $iscum = 0;
    $title  = "$outfile";
    $ylabel = "Number of alignments";
    my $label = "median_ID = $medianID";
    FUNCS::gnuplot_histo($statfile, 4, 1, $hidfileps, $title, $xlabel, $ylabel, $label, $iscum, 1);

    $label = "median_MATCH = $medianMATCH";
    $xlabel = "% MATCH";
    FUNCS::gnuplot_histo($statfile, 4, 2, $hmatchfileps, $title, $xlabel, $ylabel, $label, $iscum, $viewplots);

    $xlabel = "% ID of aligned region";
    $ylabel = "SEN (\%)";
    $key = "";
    FUNCS::gnuplot_ave_histo_with_dots("SEN", $outfile, $statfile, $Sfileps, $title, $key, $xlabel, $ylabel, 0.0, 100.0, 0.0, 100.0, $viewplots);
    $ylabel = "PPV (\%)";
    $key = "";
    FUNCS::gnuplot_ave_histo_with_dots("PPV", $outfile, $statfile, $Pfileps, $title, $key, $xlabel, $ylabel, 0.0, 100.0, 0.0, 100.0, $viewplots);
    $ylabel = "F (\%)";
    $key = "";
    FUNCS::gnuplot_ave_histo_with_dots("F",   $outfile, $statfile, $Ffileps, $title, $key, $xlabel, $ylabel, 0.0, 100.0, 0.0, 100.0, $viewplots);
    $ylabel = "HOMOLOGY SPE (\%)";
    $key = "";
    FUNCS::gnuplot_ave_histo_with_dots("SPE", $outfile, $statfile, $SPEfileps, $title, $key, $xlabel, $ylabel, 0.0, 100.0, 0.0, 100.0, $viewplots);
    
    $ylabel = "CSEN (\%)";
    $key = "";
    FUNCS::gnuplot_ave_histo_with_dots("CSEN", $outfile, $statfile, $Sfileps, $title, $key, $xlabel, $ylabel, 0.0, 100.0, 0.0, 100.0, $viewplots);
    $ylabel = "CPPV (\%)";
    $key = "";
    FUNCS::gnuplot_ave_histo_with_dots("CPPV", $outfile, $statfile, $Pfileps, $title, $key, $xlabel, $ylabel, 0.0, 100.0, 0.0, 100.0, $viewplots);
    $ylabel = "CF (\%)";
    $key = "";
    FUNCS::gnuplot_ave_histo_with_dots("CF",   $outfile, $statfile, $Ffileps, $title, $key, $xlabel, $ylabel, 0.0, 100.0, 0.0, 100.0, $viewplots);

    if ($longversion) {
	$ylabel = "score per position";
	$key = "";
	FUNCS::gnuplot_ave_histo_with_dots("SCL",   $outfile, $statfile, $SCfileps, $title, $key, $xlabel, $ylabel, 0.0, 100.0, -10.0, 0.0, $viewplots);
	$ylabel = "MATCHES (%)";
	$key = "";
	FUNCS::gnuplot_ave_histo_with_dots("MATCH", $outfile, $statfile, $Mfileps,  $title, $key, $xlabel, $ylabel, 0.0, 100.0, 80.0, 100.0, $viewplots);
	$ylabel = "GAP OPEN per MATCH (%)";
	$key = "";
	FUNCS::gnuplot_ave_histo_with_dots("OPENG", $outfile, $statfile, $Ifileps,  $title, $key, $xlabel, $ylabel, 0.0, 100.0, 0.0, 10.0, $viewplots);
    }
    
}



sub plot_id2bench_dots {

    my ($file, $field, $name, $ylabel, $viewplots) = @_;

    my $psfile  = "$file.$name.dot.ps";
    my $pdffile = "$file.$name.dot.pdf";
    my $xlabel  = "\% ID of aligned region";
    
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    print GP "set output '$psfile'\n";
    print GP "set key right top\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    print GP "set yrange  [0:100]\n";
    print GP "set xrange [100:0]\n";
    print GP "set title '$file'\n";

    #print GP "set nokey\n";
    my $cmd = "";
    #$cmd  = "'$file' using 1:$field title 'reference msa' ls 1, '$file' using 2:$field title 'infered msa' ls 2";
    $cmd  = "'$file' using 1:$field title 'reference msa' ls 1";
    print GP "plot $cmd\n";
   
    close(GP);

    system("ps2pdf $psfile\n");
    system("rm $psfile\n");
    if ($viewplots) { system("open $pdffile&\n"); }

}


