#!/usr/bin/perl -w
#ptransmarkbench.pl
 
use strict;
use Class::Struct;
use lib '/groups/eddy/home/rivase/projects/evohmm/scripts';
use lib '/Users/rivase/projects/evohmm/scripts';
use FUNCS;
#use constant GNUPLOT => '/usr/bin/gnuplot';
use constant GNUPLOT => '/opt/local/bin/gnuplot';
#use constant GNUPLOT => '/sw/bin/gnuplot';

use vars qw ($opt_b $opt_E $opt_i $opt_I $opt_p $opt_t);  # required if strict used
use Getopt::Std;
getopts ('b:Ei:I:pt');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  ptransmarkbench.pl [options] N <benchfile1>..<benchfileN>  \n\n";
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
    if ($filename =~ /\/([^\/]+)\/([^\/]+).bench/)       { $filename = "$1_$2"; }
    $name .= "$filename-";
    $filename[$n] = $filename;
    print "$filename\n";
}
$file[$N-1] = shift;
$filename = $file[$N-1];
if ($filename =~ /\/([^\/]+)\/([^\/]+).bench/)       { $filename = "$1_$2"; }
$name .= "$filename";
$filename[$N-1] = $filename;
print "$filename\n";

my $minid   = 0.0;
my $maxid   = 100.0;
my $binsizeid = 1.0;
if ($opt_b) { $binsizeid = $opt_b; if ($binsizeid <= 0.0) { print "ilegal binsizeid $binsizeid\n"; } }
if ($opt_i) { $minid     = $opt_i; if ($minid     <  0.0) { print "ilegal minid $minid\n"; } }
if ($opt_I) { $maxid     = $opt_I; if ($maxid     <  0.0) { print "ilegal maxid $maxid\n"; } }

my $mindiv = 0.0;
my $maxdiv = 7.0;
my $binsizediv = 1./10.0; 

my $efficiency = 0;
if ($opt_E) { $efficiency = 1; }

my $viewplots = 1;
if ($opt_p) { $viewplots = 0; }

my $together = 0;
if ($opt_t) { $together = 1; }

# calculate efficiency
# N >= 2
# we assume that the last file is the one with the "optimal" scores
#
if ($efficiency && $N > 1) {
    ptransmarkbench_efficiency($N, \@file, \@filename, $name, $binsizeid, $viewplots);
}

if (!$efficiency) {
    if (!$together) {
	for (my $n = 0; $n < $N; $n ++) {
	    ptransmarkbench_onefile($file[$n], $minid, $maxid, $binsizeid, $mindiv, $maxdiv, $binsizediv, $viewplots);
	}
    }
    
# compare
    $viewplots = 1;
    if ($N > 1) {
	my @statfile;
	for (my $n = 0; $n < $N; $n ++) {
	    $statfile[$n] = "$file[$n].out.stat";
	}
	ptransmarkbench_together($N, \@statfile, \@filename, $name, $viewplots);
    }
}

sub ptransmarkbench_efficiency{
    my ($Nf, $file_ref, $filename_ref, $name, $binsize, $viewplots) = @_;
    
    my @effifile;
    
    my $ntest;
    my @name;
    my @optsc;
    my @opttime;     
    my @optpid;     
    benchfile_getoptsc($file_ref->[$Nf-1], \$ntest, \@name, \@optsc, \@opttime, \@optpid);
    printf("ntest = %d\n", $ntest);

    for (my $n = 0; $n < $Nf-1; $n ++) {
	$effifile[$n]  = "$file_ref->[$n].E.n$n";
	
	efficiency($effifile[$n], $file_ref->[$n], $binsize, $ntest, \@name, \@optsc, \@opttime, \@optpid);
    }
    
    plot_efficiency($Nf, \@effifile, "E");
    plot_efficiency($Nf, \@effifile, "scE");

    
    for (my $n = 0; $n < $Nf-1; $n ++) { system("rm $effifile[$n]\n"); }
}

sub plot_efficiency
{

    my ($Nf, $effifile_ref, $which) = @_;

    my $psfile  = "$name.$which.ps";
    my $pdffile = "$name.$which.pdf";
    my $xlabel  = "\% ID";
    my $ylabel;
    
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
    
    print GP "set output '$psfile'\n";
    print GP "set key right top\n";
      
    my $ymin;
    my $field;
    if ($which =~ /^E$/) {
	$ylabel  = "Efficiency (\%)";
	$field = 2;
	$ymin = 0.;
    }
    elsif ($which =~ /^scE$/) {
	$ylabel  = "Efficiency per position(\%)";
	$field = 4;
	$ymin = 70.;
    }
    my $field_std = $field + 1;
    
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    #print GP "set yrange  [$ymin:*]\n";
    print GP "set xrange  [100:0]\n";
    print GP "set key bottom\n";

    my $cmd = "";
    my $x = 2;
    for (my $n = 0; $n < $Nf-2; $n ++) {
	#$cmd  .= "'$effifile_ref->[$n]'   using 1:$field with lines title '' ls $x, ";
	$cmd  .= "'$effifile_ref->[$n]'   using 1:$field with lines title '' ls $x, '$effifile_ref->[$n]'   using 1:$field:$field_std with yerrorbars title '' ls $x, ";
	$x ++; if ($x > 10) { $x = 1; }
    }
    #$cmd  .= "'$effifile_ref->[$Nf-2]'   using 1:$field with lines title '' ls $x";
    $cmd  .= "'$effifile_ref->[$Nf-2]'   using 1:$field  with lines title '' ls $x, '$effifile_ref->[$Nf-2]'   using 1:$field:$field_std with yerrorbars title '' ls $x";
    print GP "plot $cmd\n";
    close(GP);


    system("ps2pdf $psfile\n");
    system("rm $psfile\n");
    if ($viewplots) { system("open $pdffile&\n"); }

}


sub benchfile_getoptsc{
   my ($benchfile, $ntest_ref, $name_ref, $optsc_ref, $opttime_ref, $optpid_ref) = @_;
    
   my $hmm_name;
   my $hmm_time;
   my $hmm_avgpid;
   my $hmm_me;
   my $hmm_mre;
   my $hmm_alen;
   my $hmm_avgsqlen;

   my $n = 0;

   open(BENCH,   "$benchfile")   || die;
   while (<BENCH>) {
       if    (/\#/) { next; }
       elsif (/^ehmm\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$/) {
	   $hmm_name     = $1;
	   $hmm_time     = $2;
	   $hmm_alen     = $3;
	   $hmm_avgsqlen = $4;
	   $hmm_avgpid   = $5;
	   $hmm_me       = $6;
	   $hmm_mre      = $7;
       }
       elsif (/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s*$/) {
	   my $name       = $1;
	   my $tph        = $2;
	   my $th         = $3;
	   my $fh         = $4;
	   my $nhe        = $5; # non homologous inferred
	   my $nhr        = $6; # non-homologous reference
	   my $rid        = $7;
	   my $eid        = $8;
	   my $rmatch     = $9;
	   my $ematch     = $10;
	   my $evodist    = $11;
	   my $time       = $12;
	   my $method     = $13;
	   
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
	   
	   my $scn = $sc / $hmm_alen; # length-normalized score
	   if ($method =~ /^NA$/) {
	       $name_ref->[$n]    = $name;
	       $optsc_ref->[$n]   = $sc;
	       $opttime_ref->[$n] = $hmm_time;
	       $optpid_ref->[$n]  = $hmm_avgpid;
	       $n ++;

	   } 
	   else { printf "bad parsing of benchfile\n"; die; }
      }
   }
   close (BENCH);

   $$ntest_ref = $n;
 }

sub efficiency{
    my ($effifile, $benchfile, $binsize, $ntest, $name_ref, $optsc_ref, $opttime_ref, $optpid_ref) = @_;
    
    my $N = 100;
    my $k = 1.0/$binsize;
    my $shift = 0;
    
    my @hisID;
    my @ave_EffID;
    my @std_EffID;
    my @ave_EffscaledID;
    my @std_EffscaledID;
    
    FUNCS::init_histo_array($N, $k, \@hisID);
    FUNCS::init_histo_array($N, $k, \@ave_EffID);
    FUNCS::init_histo_array($N, $k, \@std_EffID);
    FUNCS::init_histo_array($N, $k, \@ave_EffscaledID);
    FUNCS::init_histo_array($N, $k, \@std_EffscaledID);
    
    my $hmm_name;
    my $hmm_time;
    my $hmm_avgpid;
    my $hmm_me;
    my $hmm_mre;
    my $hmm_alen;
    my $hmm_avgsqlen;
    
    my $nhmm = 0;
    my $nfound = 0;

    my $nt = 0;
    my $avgEff = 0;
    my $stdEff = 0;

    my $name_prv = "";

    open(BENCH,   "$benchfile")   || die;
    while (<BENCH>) {
	if    (/\#/) { next; }
	elsif (/^ehmm\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$/) {
	    $hmm_name     = $1;
	    $hmm_time     = $2;
	    $hmm_alen     = $3;
	    $hmm_avgsqlen = $4;
	    $hmm_avgpid   = $5;
	    $hmm_me       = $6;
	    $hmm_mre      = $7;
	}
	elsif (/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s*$/) {
	    my $name       = $1;
	    my $tph        = $2;
	    my $th         = $3;
	    my $fh         = $4;
	    my $nhe        = $5; # non homologous inferred
	    my $nhr        = $6; # non-homologous reference
	    my $rid        = $7;
	    my $eid        = $8;
	    my $rmatch     = $9;
	    my $ematch     = $10;
	    my $evodist    = $11;
	    my $time       = $12;
	    my $method     = $13;
	    
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
	    
	    my $scn = $sc / $hmm_alen; # length-normalized score

	    if ($name =~ /^$name_prv$/) { } else { $nfound ++; }
	    $name_prv = $name;

	    if ($method =~ /^NA$/) {
		my $eff;
		my $effsc;
		my $n;
		for ($n = 0; $n < $ntest; $n ++) {
		    if ($name =~ /^$name_ref->[$n]$/) {
			$effsc = 100 * exp($sc - $optsc_ref->[$n]);
			$eff   = 100 * exp($sc/$hmm_alen - $optsc_ref->[$n]/$hmm_alen);
			
			$nt ++;
			$avgEff += $effsc;
			$stdEff += $effsc*$effsc;
			
			FUNCS::fill_histo_array(1,             $hmm_avgpid, $N, $k, $shift, \@hisID);
			FUNCS::fill_histo_array($eff,          $hmm_avgpid, $N, $k, $shift, \@ave_EffID);
			FUNCS::fill_histo_array($eff*$eff,     $hmm_avgpid, $N, $k, $shift, \@std_EffID);    	    
			FUNCS::fill_histo_array($effsc,        $hmm_avgpid, $N, $k, $shift, \@ave_EffscaledID);
			FUNCS::fill_histo_array($effsc*$effsc, $hmm_avgpid, $N, $k, $shift, \@std_EffscaledID);    	    
			last;
		    }
		}
		if ($n == $ntest) { last; }
		
	    } 
	    else { printf "bad parsing of benchfile\n"; die; }
	}
    }
    close (BENCH);
    
    if ($nfound != $ntest) { print "wrong number of tests is $nfound should be $ntest\n"; die; }
    if ($nt     != $ntest) { print "wrong number of tests is $nt should be $ntest\n"; die; }

   FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_EffID,       \@std_EffID,       , 0);
   FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_EffscaledID, \@std_EffscaledID, , 0);

   FUNCS::calculate_averages(\$avgEff, \$stdEff, $nt);
   printf("Eff per positions %f +\- %f\n", $avgEff, $stdEff);

   open(EFILE, ">$effifile") || die;
   my $dim = $N * $k;
   
   for (my $i=0; $i<=$dim; $i++) { 
       my $idl = $i/$k     - $shift;
       my $idh = (($i+1)/$k < $N)? ($i+1)/$k - $shift : $N-$shift;
       
      printf EFILE "%f %f %f %f %f\n", $idl+0.5*($idh-$idl), $ave_EffID[$i], $std_EffID[$i], $ave_EffscaledID[$i], $std_EffscaledID[$i]; 
   }
   close(EFILE);
}

sub ptransmarkbench_together{
    my ($N, $statfile_ref, $filename_ref, $name, $viewplots) = @_;

    my $which     = "totalF";
    FUNCS::plot_id2totalF($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "totalSEN";
    FUNCS::plot_id2totalF($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "totalPPV";
    FUNCS::plot_id2totalF($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "totalSPE";
    FUNCS::plot_id2totalF($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);

    $which     = "F";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "SEN";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "PPV";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    $which     = "SPE";
    FUNCS::plot_id2F($N, $statfile_ref, $filename_ref, $which, $name, $viewplots);
    
}

sub ptransmarkbench_onefile{
    my ($benchfile, $minid, $maxid, $binsizeid, $mindiv, $maxdiv, $binsizediv, $viewplots) = @_;

    my $outf    = "$benchfile.out";
    parse_benchfile($benchfile, $outf, $minid, $maxid);
    plot_id2divergence($outf, $viewplots);
    plot_id2bench_histo($outf,  $binsizeid,  $viewplots);    
    #plot_div2bench_histo($outf, $binsizediv, $viewplots);    
}

sub final_stats_file {
    my ($outfile, $N, $k, $shift, 
	$his_ref,
	$tph_histo_ref, $th_histo_ref, $fh_histo_ref, 
	$nhe_histo_ref, $nhr_histo_ref, 
	$ave_senhis_ref,   $std_senhis_ref, 
	$ave_ppvhis_ref,   $std_ppvhis_ref, 
	$ave_Fhis_ref,     $std_Fhis_ref, 
	$ave_SPEhis_ref,   $std_SPEhis_ref) = @_;
    
    open(OUTF, ">$outfile") || die;
    printf "$outfile\n";

    my $dim = $N * $k;
    
    my $AUCsen = 0.0;
    my $AUCppv = 0.0;
    my $AUCF   = 0.0;
    my $AUCSPE = 0.0;
 
    my $AUCavgsen = 0.0;
    my $AUCavgppv = 0.0;
    my $AUCavgF   = 0.0;
    my $AUCavgSPE = 0.0;
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
	

	my $nmsa = $his_ref->[$i];
	
	if ($nmsa > 0) {
	    $num ++;
	    $AUCsen += $sen;
	    $AUCppv += $ppv;
	    $AUCF   += $F;
	    $AUCSPE += $SPE;
	    
	    $AUCavgsen += $ave_senhis_ref->[$i];
	    $AUCavgppv += $ave_ppvhis_ref->[$i];
	    $AUCavgF   += $ave_Fhis_ref->[$i];
	    $AUCavgSPE += $ave_SPEhis_ref->[$i];
	}

	printf OUTF "%d %d %f %f ", $nmsa, $his_ref->[$i], $idl, $idl + 0.5*($idh-$idl); 
	if ($nmsa > 0) {
	    printf OUTF "%f %f %f %f %f %f %f %f %f %f %f %f",
	    $sen, $ave_senhis_ref->[$i], $std_senhis_ref->[$i], 
	    $ppv, $ave_ppvhis_ref->[$i], $std_ppvhis_ref->[$i], 
	    $F,   $ave_Fhis_ref->[$i],   $std_Fhis_ref->[$i],
	    $SPE, $ave_SPEhis_ref->[$i], $std_SPEhis_ref->[$i]; 
	}
	printf OUTF "\n"; 
	
    }

    my $denom = $num;
    $denom = ($denom>0)? 1.0/$denom : 1.0;

    $AUCsen *= $denom;
    $AUCppv *= $denom;
    $AUCF   *= $denom;
    $AUCSPE *= $denom;
    
    $AUCavgsen *= $denom;
    $AUCavgppv *= $denom;
    $AUCavgF   *= $denom;
    $AUCavgSPE *= $denom;
    
    printf OUTF "#AUC %.1f %.1f %.1f %.1f | %.1f %.1f %.1f %.1f\n", $AUCsen, $AUCppv, $AUCF, $AUCSPE, $AUCavgsen, $AUCavgppv, $AUCavgF, $AUCavgSPE;  
    printf      "#AUC %.1f %.1f %.1f %.1f | %.1f %.1f %.1f %.1f\n", $AUCsen, $AUCppv, $AUCF, $AUCSPE, $AUCavgsen, $AUCavgppv, $AUCavgF, $AUCavgSPE; 

    close (OUTF);
    
}


sub parse_benchfile{
    my ($benchfile, $outfile, $minid, $maxid) = @_;
    
    my $ttph = 0;
    my $tth  = 0;
    my $tfh  = 0;
    my $tnhr = 0;
    my $tnhe = 0;
     
    my $ttime    = 0.;

    my $nmsa       = 0;
    my $nmsatot    = 0;

    my $ignore;

    my $hmm_name;
    my $hmm_time;
    my $hmm_avgpid;
    my $hmm_me;
    my $hmm_mre;
    my $hmm_alen;
    my $hmm_avgsqlen;

    my $hmm_name_prv = "";
    my $nhmm = 0;

    my $name_prv = "";
    my $ntest = 0;

    open(OUTF,    ">$outfile")    || die;
    open(BENCH,   "$benchfile")   || die;
    while (<BENCH>) {
	if    (/\#/) { next; }
	elsif (/^ehmm\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$/) {
	    $hmm_name     = $1;
	    $hmm_time     = $2;
	    $hmm_alen     = $3;
	    $hmm_avgsqlen = $4;
	    $hmm_avgpid   = $5;
	    $hmm_me       = $6;
	    $hmm_mre      = $7;
	    if ($hmm_name =~ /^$hmm_name_prv$/) { } else { $nhmm ++; }
	    $hmm_name_prv = $hmm_name;
	}
	elsif (/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s*$/) {
	    my $name       = $1;
	    my $tph        = $2;
	    my $th         = $3;
	    my $fh         = $4;
	    my $nhe        = $5; # non homologous inferred
	    my $nhr        = $6; # non-homologous reference
	    my $rid        = $7;
	    my $eid        = $8;
	    my $rmatch     = $9;
	    my $ematch     = $10;
	    my $evodist    = $11;
	    my $time       = $12;
	    my $method     = $13;

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

	    if ($name =~ /^$hmm_name\-\d+$/) { } else { printf "hmm %s?? test %s\n", $hmm_name, $name; die; }
	    if ($method =~ /^NA$/) { } else { printf "bad parsing of benchfile\n"; die; }

	    if ($name =~ /^$name_prv$/) { } else { $ntest ++; }
	    $name_prv = $name;

	    my $s;
	    my $p;
	    my $F;
	    FUNCS::calculateF($tph, $th, $fh,  \$s,  \$p,  \$F);

	    my $SPE = ($nhr > 0)? 100.0*$nhe/$nhr : 0.0;

	    $nmsatot ++;
	    if ($hmm_avgpid < $minid || $hmm_avgpid > $maxid) { next; }
	    $nmsa    ++;
	    $ttph    += $tph;
	    $tth     += $th;
	    $tfh     += $fh;
	    $tnhe    += $nhe;
	    $tnhr    += $nhr;
	    $ttime   += $time; 
	    
	    printf OUTF "%f %f %f %f %d %d %d %d %d %f %f %f %f %f\n", 
	    $hmm_time, $hmm_avgpid, $hmm_me, $hmm_mre, 
	    $tph, $th, $fh, $nhe, $nhr, 
	    $s, $p, $F, $SPE, $sc;
	}
    }
    close (BENCH);
    close (OUTF);
    
    printf("Nhmm     $nhmm\n");
    printf("Ntestmsa $ntest\n");

    my $sen    = 0.;
    my $ppv    = 0.;
    my $F      = 0.;
    my $SPE    = 0.;
    
    FUNCS::calculateF($ttph,    $tth,    $tfh,    \$sen,    \$ppv,    \$F);
    
    $SPE    = ($tnhr    > 0.)? 100.*$tnhe/$tnhr       : 0.0;
 
    printf("nmsa %d/%d id [%f,%f] total_homologies %d\n", $nmsa, $nmsatot, $minid, $maxid, $tth);
    printf("sen %.2f\tppv %.2f\tF %.2f\t\tSPE %.2f\t\tNA\tCPU Time (s) %.f\n",   $sen,    $ppv,    $F,    $SPE,    $ttime);
    
}


sub plot_id2divergence {

    my ($file, $viewplots) = @_;

    my $psfile  = "$file.ps";
    my $pdffile = "$file.pdf";
    my $xlabel;
    my $ylabel;
    my $cmd;
    my $field;

    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    print GP "set terminal postscript color 14\n";
    FUNCS::gnuplot_define_styles (*GP);
 
    print GP "set output '$psfile'\n";
    #print GP "set yrange  [0:100]\n";

    # x-axis pid
    $xlabel  = "HMM avg \% ID";
    print GP "set xlabel '$xlabel'\n";
    print GP "set xrange [100:0]\n";
    $ylabel = "HMM Mean entropy";
    print GP "set yrange [0:*]\n";
    print GP "set ylabel '$ylabel'\n";
    $field = 3;
    $cmd  = "'$file' using 2:$field title '' ls 2";
    print GP "plot $cmd\n";

    $ylabel = "HMM Mean Relative entropy";
    print GP "set ylabel '$ylabel'\n";
    $field = 4;
    $cmd  = "'$file' using 2:$field title '' ls 2";
    print GP "plot $cmd\n";

    $ylabel = "divergence time";
    print GP "set ylabel '$ylabel'\n";
    print GP "set yrange [0:10]\n";
    $field = 1;
    $cmd  = "'$file' using 2:$field title '' ls 2";
    print GP "plot $cmd\n";

    # x-axis divergence time
    $xlabel  = "HMM divergence time";
    print GP "set xrange [0:3]\n";
    print GP "set xlabel '$xlabel'\n";

    $ylabel = "HMM avg \%ID";
    print GP "set ylabel '$ylabel'\n";
    print GP "set yrange [100:0]\n";
    $field = 2;
    $cmd  = "'$file' using 1:$field title '' ls 2";
    print GP "plot $cmd\n";

    $ylabel = "HMM Mean entropy";
    print GP "set ylabel '$ylabel'\n";
    print GP "set yrange [0:4]\n";
    $field = 3;
    $cmd  = "'$file' using 1:$field title '' ls 2";
    print GP "plot $cmd\n";

    $ylabel = "HMM Mean Relative entropy";
    print GP "set ylabel '$ylabel'\n";
    $field = 4;
    $cmd  = "'$file' using 1:$field title '' ls 2";
    print GP "plot $cmd\n";
    close(GP);

    system("ps2pdf $psfile\n");
    system("rm $psfile\n");
    if ($viewplots) { system("open $pdffile&\n"); }
}


sub plot_id2bench_histo {

    my ($outfile, $binsize, $viewplots) = @_;

    my $statfile = "$outfile.IDstat";

    my $Sfileps   = "$statfile.SEN.ps";
    my $Pfileps   = "$statfile.PPV.ps";
    my $Ffileps   = "$statfile.F.ps";
    my $SPEfileps = "$statfile.SPE.ps";

    my $hidfileps    = "$statfile.histoID.ps";

    my $xlabel = "\% ID";
    my $ylabel = "";
    my $title  = "";
    my $key    = "";

    #histogram
    my $N     = 100;
    my $k     = 1.0/$binsize;
    my $shift = 0;

    my @hisID;

    my @tph_histo;
    my @th_histo;
    my @fh_histo;
    my @nhe_histo;
    my @nhr_histo;

    my @ave_senhisID;
    my @std_senhisID;
    my @ave_ppvhisID;
    my @std_ppvhisID;
    my @ave_FhisID;
    my @std_FhisID;
    my @ave_SPEhisID;
    my @std_SPEhisID;
    FUNCS::init_histo_array($N, $k, \@hisID);
    FUNCS::init_histo_array($N, $k, \@tph_histo);
    FUNCS::init_histo_array($N, $k, \@th_histo);
    FUNCS::init_histo_array($N, $k, \@fh_histo);
    FUNCS::init_histo_array($N, $k, \@nhe_histo);
    FUNCS::init_histo_array($N, $k, \@nhr_histo);
    FUNCS::init_histo_array($N, $k, \@ave_senhisID);
    FUNCS::init_histo_array($N, $k, \@std_senhisID);
    FUNCS::init_histo_array($N, $k, \@ave_ppvhisID);
    FUNCS::init_histo_array($N, $k, \@std_ppvhisID);
    FUNCS::init_histo_array($N, $k, \@ave_FhisID);
    FUNCS::init_histo_array($N, $k, \@std_FhisID);
    FUNCS::init_histo_array($N, $k, \@ave_SPEhisID);
    FUNCS::init_histo_array($N, $k, \@std_SPEhisID);

    open(FILE, "$outfile")   || die;
    while (<FILE>) {
	if    (/\#/) { next; }
	elsif (/^(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*(\S+)$/) {
	    my $time = $1;
	    my $id   = $2;
	    my $tph  = $3;
	    my $th   = $4;
	    my $fh   = $5;
	    my $nhe  = $6;
	    my $nhr  = $7;
	    my $sen  = $8;
	    my $ppv  = $9;
	    my $F    = $10;
	    my $SPE  = $11;
	    my $sc  = $12;

	    FUNCS::fill_histo_array(1,         $id, $N, $k, $shift, \@hisID);
	    FUNCS::fill_histo_array($tph,      $id, $N, $k, $shift, \@tph_histo);
	    FUNCS::fill_histo_array($th,       $id, $N, $k, $shift, \@th_histo);
	    FUNCS::fill_histo_array($fh,       $id, $N, $k, $shift, \@fh_histo);
	    FUNCS::fill_histo_array($nhe,      $id, $N, $k, $shift, \@nhe_histo);
	    FUNCS::fill_histo_array($nhr,      $id, $N, $k, $shift, \@nhr_histo);

	    FUNCS::fill_histo_array($sen,      $id, $N, $k, $shift, \@ave_senhisID);
	    FUNCS::fill_histo_array($sen*$sen, $id, $N, $k, $shift, \@std_senhisID);    	    
	    FUNCS::fill_histo_array($ppv,      $id, $N, $k, $shift, \@ave_ppvhisID);
	    FUNCS::fill_histo_array($ppv*$ppv, $id, $N, $k, $shift, \@std_ppvhisID);    	    
	    FUNCS::fill_histo_array($F,        $id, $N, $k, $shift, \@ave_FhisID);
	    FUNCS::fill_histo_array($F*$F,     $id, $N, $k, $shift, \@std_FhisID);  
	    FUNCS::fill_histo_array($SPE,      $id, $N, $k, $shift, \@ave_SPEhisID);
	    FUNCS::fill_histo_array($SPE*$SPE, $id, $N, $k, $shift, \@std_SPEhisID);  
 	}
    }
    close(FILE);
    my $medianID    = FUNCS::histogram_median($N, $k, \@hisID);

    FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_senhisID, \@std_senhisID, , 0);
    FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_ppvhisID, \@std_ppvhisID, , 0);
    FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_FhisID,   \@std_FhisID,   , 0);
    FUNCS::write_ave_histogram($N, $k, $shift, \@hisID, \@ave_SPEhisID, \@std_SPEhisID,   , 0);
    
    final_stats_file($statfile, $N, $k, $shift, 
		     \@hisID, 
		     \@tph_histo,       \@th_histo, \@fh_histo, 
		     \@nhe_histo,       \@nhr_histo, 
		     \@ave_senhisID,    \@std_senhisID, 
		     \@ave_ppvhisID,    \@std_ppvhisID, 
		     \@ave_FhisID,      \@std_FhisID,  
		     \@ave_SPEhisID,    \@std_SPEhisID);
    
    printf("median_ID = $medianID\n");

    my $iscum = 0;
    $title  = "$outfile";
    $ylabel = "Number of alignments";
    my $label = "median_ID = $medianID";
    FUNCS::gnuplot_histo($statfile, 4, 1, $hidfileps, $title, $xlabel, $ylabel, $label, $iscum, 1, 100.0, 0.0);

    $xlabel = "% ID";
    $ylabel = "SEN (\%)";

    my $xmin = 100.0;
    my $xmax = 0.;
    my $ymin = 0.0;
    my $ymax = 100.;
    $key = "";
    gnuplot_ave_histoID_with_dots("SEN", $outfile, $statfile, $Sfileps, $title, $key, $xlabel, $ylabel, $xmin, $xmax, $ymin, $ymax, $viewplots);
    $ylabel = "PPV (\%)";
    $key = "";
    gnuplot_ave_histoID_with_dots("PPV", $outfile, $statfile, $Pfileps, $title, $key, $xlabel, $ylabel, $xmin, $xmax, $ymin, $ymax, $viewplots);
    $ylabel = "F (\%)";
    $key = "";
    gnuplot_ave_histoID_with_dots("F",   $outfile, $statfile, $Ffileps, $title, $key, $xlabel, $ylabel, $xmin, $xmax, $ymin, $ymax, $viewplots);
    $ylabel = "HOMOLOGY SPE (\%)";
    $key = "";
    gnuplot_ave_histoID_with_dots("SPE", $outfile, $statfile, $SPEfileps, $title, $key, $xlabel, $ylabel, $xmin, $xmax, $ymin, $ymax, $viewplots);    
}

sub plot_div2bench_histo {

    my ($outfile, $binsize, $viewplots) = @_;

    my $statfile   = "$outfile.DIVstat";

    my $Sfileps    = "$statfile.SEN.ps";
    my $Pfileps    = "$statfile.PPV.ps";
    my $Ffileps    = "$statfile.F.ps";
    my $SPEfileps  = "$statfile.SPE.ps";

    my $hdivfileps = "$statfile.histoDIV.ps";

    my $xlabel = "divergence time";
    my $ylabel = "";
    my $title  = "";
    my $key    = "";

    #histogram
    my $N     = $maxdiv-$mindiv;
    my $k     = 1.0/$binsize;
    my $shift = 0;

    my @hisDIV;

    my @tph_histo;
    my @th_histo;
    my @fh_histo;
    my @nhe_histo;
    my @nhr_histo;

    my @ave_senhisDIV;
    my @std_senhisDIV;
    my @ave_ppvhisDIV;
    my @std_ppvhisDIV;
    my @ave_FhisDIV;
    my @std_FhisDIV;
    my @ave_SPEhisDIV;
    my @std_SPEhisDIV;
    FUNCS::init_histo_array($N, $k, \@hisDIV);
    FUNCS::init_histo_array($N, $k, \@tph_histo);
    FUNCS::init_histo_array($N, $k, \@th_histo);
    FUNCS::init_histo_array($N, $k, \@fh_histo);
    FUNCS::init_histo_array($N, $k, \@nhe_histo);
    FUNCS::init_histo_array($N, $k, \@nhr_histo);
    FUNCS::init_histo_array($N, $k, \@ave_senhisDIV);
    FUNCS::init_histo_array($N, $k, \@std_senhisDIV);
    FUNCS::init_histo_array($N, $k, \@ave_ppvhisDIV);
    FUNCS::init_histo_array($N, $k, \@std_ppvhisDIV);
    FUNCS::init_histo_array($N, $k, \@ave_FhisDIV);
    FUNCS::init_histo_array($N, $k, \@std_FhisDIV);
    FUNCS::init_histo_array($N, $k, \@ave_SPEhisDIV);
    FUNCS::init_histo_array($N, $k, \@std_SPEhisDIV);

    open(FILE, "$outfile")   || die;
    while (<FILE>) {
	if    (/\#/) { next; }
	elsif (/^(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*(\S+)$/) {
	    my $time = $1;
	    my $id   = $2;
	    my $tph  = $3;
	    my $th   = $4;
	    my $fh   = $5;
	    my $nhe  = $6;
	    my $nhr  = $7;
	    my $sen  = $8;
	    my $ppv  = $9;
	    my $F    = $10;
	    my $SPE  = $11;
	    my $sc   = $12;

	    my $div = $time;

	    FUNCS::fill_histo_array(1,         $div, $N, $k, $shift, \@hisDIV);
	    FUNCS::fill_histo_array($tph,      $div, $N, $k, $shift, \@tph_histo);
	    FUNCS::fill_histo_array($th,       $div, $N, $k, $shift, \@th_histo);
	    FUNCS::fill_histo_array($fh,       $div, $N, $k, $shift, \@fh_histo);
	    FUNCS::fill_histo_array($nhe,      $div, $N, $k, $shift, \@nhe_histo);
	    FUNCS::fill_histo_array($nhr,      $div, $N, $k, $shift, \@nhr_histo);

	    FUNCS::fill_histo_array($sen,      $div, $N, $k, $shift, \@ave_senhisDIV);
	    FUNCS::fill_histo_array($sen*$sen, $div, $N, $k, $shift, \@std_senhisDIV);    	    
	    FUNCS::fill_histo_array($ppv,      $div, $N, $k, $shift, \@ave_ppvhisDIV);
	    FUNCS::fill_histo_array($ppv*$ppv, $div, $N, $k, $shift, \@std_ppvhisDIV);    	    
	    FUNCS::fill_histo_array($F,        $div, $N, $k, $shift, \@ave_FhisDIV);
	    FUNCS::fill_histo_array($F*$F,     $div, $N, $k, $shift, \@std_FhisDIV);  
	    FUNCS::fill_histo_array($SPE,      $div, $N, $k, $shift, \@ave_SPEhisDIV);
	    FUNCS::fill_histo_array($SPE*$SPE, $div, $N, $k, $shift, \@std_SPEhisDIV);  
 	}
    }
    close(FILE);    my $medianDIV    = FUNCS::histogram_median($N, $k, \@hisDIV);

    FUNCS::write_ave_histogram($N, $k, $shift, \@hisDIV, \@ave_senhisDIV, \@std_senhisDIV, , 0);
    FUNCS::write_ave_histogram($N, $k, $shift, \@hisDIV, \@ave_ppvhisDIV, \@std_ppvhisDIV, , 0);
    FUNCS::write_ave_histogram($N, $k, $shift, \@hisDIV, \@ave_FhisDIV,   \@std_FhisDIV,   , 0);
    FUNCS::write_ave_histogram($N, $k, $shift, \@hisDIV, \@ave_SPEhisDIV, \@std_SPEhisDIV,   , 0);
    
    final_stats_file($statfile, $N, $k, $shift, 
		     \@hisDIV, 
		     \@tph_histo,       \@th_histo, \@fh_histo, 
		     \@nhe_histo,       \@nhr_histo, 
		     \@ave_senhisDIV,    \@std_senhisDIV, 
		     \@ave_ppvhisDIV,    \@std_ppvhisDIV, 
		     \@ave_FhisDIV,      \@std_FhisDIV,  
		     \@ave_SPEhisDIV,    \@std_SPEhisDIV);
    
    printf("median_DIV = $medianDIV\n");

    my $iscum = 0;
    $title  = "$outfile";
    $ylabel = "Number of alignments";
    my $label = "median_DIV = $medianDIV";
    FUNCS::gnuplot_histo($statfile, 4, 1, $hdivfileps, $title, $xlabel, $ylabel, $label, $iscum, 1, 0.0, 10.0);

    $xlabel = "time divergence";
    $ylabel = "SEN (\%)";

    my $xmin = $mindiv;
    my $xmax = $maxdiv;
    my $ymin = 0.0;
    my $ymax = 100.;
    $key = "";
    gnuplot_ave_histoDIV_with_dots("SEN", $outfile, $statfile, $Sfileps, $title, $key, $xlabel, $ylabel, $xmin, $xmax, $ymin, $ymax, $viewplots);
    $ylabel = "PPV (\%)";
    $key = "";
    gnuplot_ave_histoDIV_with_dots("PPV", $outfile, $statfile, $Pfileps, $title, $key, $xlabel, $ylabel, $xmin, $xmax, $ymin, $ymax, $viewplots);
    $ylabel = "F (\%)";
    $key = "";
    gnuplot_ave_histoDIV_with_dots("F",   $outfile, $statfile, $Ffileps, $title, $key, $xlabel, $ylabel, $xmin, $xmax, $ymin, $ymax, $viewplots);
    $ylabel = "HOMOLOGY SPE (\%)";
    $key = "";
    gnuplot_ave_histoDIV_with_dots("SPE", $outfile, $statfile, $SPEfileps, $title, $key, $xlabel, $ylabel, $xmin, $xmax, $ymin, $ymax, $viewplots);    
}


sub gnuplot_ave_histoID_with_dots {

    my ($which, $file, $hfile, $psfile, $title, $key, $xlabel, $ylabel, $xmin, $xmax, $ymin, $ymax, $viewplot) = @_;

    my $n;
    my $cmd;

    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    
    my $pdffile = $psfile;
    if ($pdffile =~ /^(\S+).ps$/) { $pdffile = "$1.pdf"; }

    print GP "set terminal postscript color solid 14\n";
    FUNCS::gnuplot_define_styles (*GP);
    
    print GP "set output '$psfile'\n";
    #print GP "set nokey\n";
    print GP "set xlabel '$xlabel'\n";
    #print GP "set title \"$name2\\n\\n$title\"\n";
    print GP "set title '$title'\n";
    print GP "set ylabel '$ylabel'\n";
    #print GP "set boxwidth 5.0\n";
  
    print GP "set xrange [$xmin:$xmax]\n";
    print GP "set yrange [$ymin:$ymax]\n";

    # plot 
    if ($which =~ /^SEN$/ || $which =~ /^PPV$/ ||$which =~ /^F$/ || $which =~ /^SPE$/) {
	my $field;
	my $field1;
	
	if    ($which =~/^SEN$/) { $field = 10; $field1 = 5;  }
	elsif ($which =~/^PPV$/) { $field = 11; $field1 = 8;  }
	elsif ($which =~/^F$/)   { $field = 12; $field1 = 11; }
	elsif ($which =~/^SPE$/) { $field = 13; $field1 = 14; }
	my $field_avg = $field1 + 1;
	my $field_std = $field1 + 2;
	$cmd = "";
	
	$cmd .= "'$file'  using 2:$field                title '' ls 1, ";
	$cmd .= "'$hfile' using 4:$field1               with boxes title '' ls 3, ";
	$cmd .= "'$hfile' using 4:$field1               with lines title '' ls 7, ";
	$cmd .= "'$hfile' using 4:$field_avg:$field_std with yerrorbars title '$key' ls 2";
	print GP "plot $cmd\n";
	
	close (GP);
    }

    system("ps2pdf $psfile\n");
    system("rm $psfile\n");
    if ($viewplot) { 
	system ("open $pdffile&\n"); 
    }
}

sub gnuplot_ave_histoDIV_with_dots {

    my ($which, $file, $hfile, $psfile, $title, $key, $xlabel, $ylabel, $xmin, $xmax, $ymin, $ymax, $viewplot) = @_;

    my $n;
    my $cmd;

    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    
    my $pdffile = $psfile;
    if ($pdffile =~ /^(\S+).ps$/) { $pdffile = "$1.pdf"; }

    print GP "set terminal postscript color solid 14\n";
    FUNCS::gnuplot_define_styles (*GP);
    
    print GP "set output '$psfile'\n";
    #print GP "set nokey\n";
    print GP "set xlabel '$xlabel'\n";
    #print GP "set title \"$name2\\n\\n$title\"\n";
    print GP "set title '$title'\n";
    print GP "set ylabel '$ylabel'\n";
    #print GP "set boxwidth 5.0\n";
  
    print GP "set xrange [$xmin:$xmax]\n";
    print GP "set yrange [$ymin:$ymax]\n";

    # plot 
    if ($which =~ /^SEN$/ || $which =~ /^PPV$/ ||$which =~ /^F$/ || $which =~ /^SPE$/) {
	my $field;
	my $field1;
	
	if    ($which =~/^SEN$/) { $field = 10; $field1 = 5;  }
	elsif ($which =~/^PPV$/) { $field = 11; $field1 = 8;  }
	elsif ($which =~/^F$/)   { $field = 12; $field1 = 11; }
	elsif ($which =~/^SPE$/) { $field = 13; $field1 = 14; }
	my $field_avg = $field1 + 1;
	my $field_std = $field1 + 2;
	$cmd = "";
	
	$cmd .= "'$file'  using 1:$field                title '' ls 1, ";
	$cmd .= "'$hfile' using 4:$field1               with boxes title '' ls 3, ";
	$cmd .= "'$hfile' using 4:$field1               with lines title '' ls 7, ";
	$cmd .= "'$hfile' using 4:$field_avg:$field_std with yerrorbars title '$key' ls 2";
	print GP "plot $cmd\n";
	
	close (GP);
    }

    system("ps2pdf $psfile\n");
    system("rm $psfile\n");
    if ($viewplot) { 
	system ("open $pdffile&\n"); 
    }
}
