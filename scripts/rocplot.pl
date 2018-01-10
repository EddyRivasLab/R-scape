#!/usr/bin/perl -w
# rocplot.pl

use strict;
use Class::Struct;
use POSIX;


# find directory where the script is installed
use FindBin;
use lib $FindBin::Bin;
use PDBFUNCS;
use FUNCS;

use vars qw ($opt_C $opt_D $opt_f $opt_L $opt_p $opt_P $opt_R $opt_r $opt_s $opt_T $opt_W $opt_v $opt_x);  # required if strict used
use Getopt::Std;
getopts ('C:D:f:L:p:PRrs:T:W:vx:');


# Print a helpful message if the user provides no input file.
if (!@ARGV) {
    print "usage:  rocplot.pl [options] <F> <file1> <stofile1> <strfile1> .. <fileF> <stofileF> <strfileF> <rscapebin> <gnuplot> <key>  \n\n";
    print "options:\n";
    exit;
}

my $F = shift;
my @prefile;
my @stofile;
my @strfile;
my @prename;
my @rocfile;
my @scat1file;
my @scat2file;
my $ID = 1.0;
my $outname;
 for (my $f = 0; $f < $F; $f ++){
    $prefile[$f]  = shift;
    $stofile[$f]  = shift;
    $strfile[$f]  = shift;
    $prename[$f]  = $prefile[$f];
    if ($prename[$f] =~ /^(\S+)\.sorted.out$/) {
	$outname = $1;
 
	$prename[$f] = "$outname.rscape";
	if ($outname     =~ /\/([^\/]+)\s*$/) { $outname = $1; }
	if ($prename[$f] =~ /\/(ID\S+)\//)    { $ID = $1; $outname .= ".$ID"; }
    }
}
my $stoname  = $prefile[0];
if ($stoname =~ /\/([^\/]+)$/) { $stoname = $1; }
if ($stoname =~ /([^\.]+)\./)  { $stoname = $1; }

my $rscapebin = shift;
my $gnuplot   = shift;
my $key       = shift; # tag to identify the benchmark

my $byali = 0;
#print "BYALI $byali\n";

my $currdir = $ENV{PWD};

my $maxD = 8;
if ($opt_D) { $maxD = $opt_D; }

my $minL = 1;
if ($opt_L) { $minL = $opt_L; }

my $target_factor = 1.5;
if ($opt_T) { $target_factor = $opt_T; }

my $dornaview = 0;
if ($opt_r) { $dornaview = 1; }

my $which = "MIN"; #options: CA CB MIN AVG NOH / C1' (for RNA suggested by Westhof)
if ($opt_W) { $which = "$opt_W"; }

# name of the output files
# .roc
for (my $f = 0; $f < $F; $f ++) {  
    $rocfile[$f]   = "$prename[$f].maxD$maxD.minL$minL.type$which.roc"; 
    $scat1file[$f] = "$prename[$f].maxD$maxD.minL$minL.type$which.scat1"; 
    $scat2file[$f] = "$prename[$f].maxD$maxD.minL$minL.type$which.scat2";
}

my $usechain = "";
if ($opt_C) { $usechain = "$opt_C"; }

my $seeplots = 0;

my $verbose = 0;
if ($opt_v) { $verbose = 1; }

# Map the  pdb2msa structure to the input alignment
my @pdb2msa;
for (my $f = 0; $f < $F; $f ++) {
    if ($strfile[$f]) { 
	my $found = 0;
	for (my $g = 0; $g < $f; $g ++) {
	    if ($stofile[$f] =~ /^$stofile[$g]$/ && $strfile[$f] =~ /^$strfile[$g]$/) {
		$pdb2msa[$f] = $pdb2msa[$g];
		$found = 1;
		last;
	    }
	}
	if ($found == 0) {
	    if ($strfile[$f] =~ /\.pdb/) {
		PDBFUNCS::pdb2msa($gnuplot, $rscapebin, $strfile[$f], $stofile[$f], \$pdb2msa[$f], $usechain, $maxD, $minL, $byali, $which, $dornaview, $seeplots);
	    }
	    elsif ($strfile[$f] =~ /\.EC\.interaction/) {
		structure_from_ecfile($strfile[$f], \$pdb2msa[$f], $maxD, $minL);
		$byali = -1; # minL cutoff already used in the annotation
	    }
	    elsif ($strfile[$f] =~ /\.cmap/) {
		structure_from_contactmapfile($strfile[$f], \$pdb2msa[$f], $maxD, $minL);
		$byali = -1; # minL cutoff already used in the annotation
	    }
	}
    }
    else {
	structure_from_msa($stofile[$f], \$pdb2msa[$f], $maxD, $minL);
    }
    
    printf "ncnt %d nbp %d nwc %d\n", $pdb2msa[$f]->{"PDB2MSA::ncnt"}, $pdb2msa[$f]->{"PDB2MSA::nbp"}, $pdb2msa[$f]->{"PDB2MSA::nwc"} ;	
}

# add a random file
my $dorandom = ($opt_P)? 1:0;
if ($dorandom) {
    $pdb2msa[$F]   = $pdb2msa[0];
    $prename[$F]   = "results/random/$stoname.random";
    $prefile[$F]   = "results/random/$stoname.random";
    $rocfile[$F]   = "results/random/$stoname.random.maxD$maxD.minL$minL.type$which.roc";
    $scat1file[$F] = "results/random/$stoname.random.maxD$maxD.minL$minL.type$which.scat1";
    $scat2file[$F] = "results/random/$stoname.random.maxD$maxD.minL$minL.type$which.scat2";
    $F ++;
}


my $alenDCA = -1;
my @mapDCA;
for (my $f = 0; $f < $F; $f ++) {

    my $pdbfile  = ($strfile[$f] =~ /\.pdb/)? 1 : 0;
    my $ecfile   = ($strfile[$f] =~ /\.EC\.interaction/)? 1 : 0;
    my $cmapfile = ($strfile[$f] =~ /\.cmap/)? 1 : 0;
    
    my $exist_rocfile = (-e $rocfile[$f])? 1 : 0;

    my $target_ncnt = 1e+10;
    if ($pdbfile =~ /\S+/) { $target_ncnt = floor($target_factor*$pdb2msa[$f]{"PDB2MSA::pdblen"}); }
    
    my $mapfile_pred = "$prename[$f].maxD$maxD.minL$minL.type$which.pred.map"; 
    my $mapfile_tp   = "$prename[$f].maxD$maxD.minL$minL.type$which.tp.map"; 
    
    my $method = "";
    if    ($prefile[$f] =~ /results\/(\S+)_filtered\//) { $method = $1; }
    elsif ($prefile[$f] =~ /results\/([^\/]+)\//)       { $method = $1; }
    elsif ($prefile[$f] =~ /test\/([^\/]+)\//)          { $method = $1; }
    
    my @his;
    my $N = 1000;
    my $k = 1;
    my $shift = 0;
    my $fmax = 100;
    FUNCS::init_histo_array($N, $k, \@his);

    print "\n$method: $prefile[$f]\n";
    if ($method =~ /^R-scape/ || $method =~ /^GTp/ || $method =~ /^PTFp/ || $method =~ /^neogremlin/) {
	## both functions below should produce exactly the same results (only if the minL used running R-scape is the same than here)
	
	if ($pdbfile || $ecfile || $cmapfile) {
	    if (!$exist_rocfile || !-e $mapfile_pred || !-e $mapfile_tp) {
		print "         $rocfile[$f]\n";
		create_rocfile_rscape_withpdb($rocfile[$f], $scat1file[$f], $scat2file[$f], $mapfile_pred, $mapfile_tp, $prefile[$f], $pdb2msa[$f], $target_ncnt, 
					      $N, $k, $shift, \@his, $fmax);
	    }
	}
	else {
	    if (!$exist_rocfile || !-e $mapfile_pred || !-e $mapfile_tp) { 
		print "         $rocfile[$f]\n";
		create_rocfile_rscape($rocfile[$f], $scat1file[$f], $scat2file[$f], $mapfile_pred, $mapfile_tp, $prefile[$f], $target_ncnt, 
				      $N, $k, $shift, \@his, $fmax);
	    }
	    $dorandom = 0;
	}
	if ($pdbfile) { predictions_plot($mapfile_pred, $mapfile_tp, $pdb2msa[$f], $target_ncnt); }
    }
    elsif ($method =~ /^mfDCA$/) {
	if (!$exist_rocfile || !-e $mapfile_pred || !-e $mapfile_tp) { 
	    print "         $rocfile[$f]\n";
	    create_rocfile_mfDCA($rocfile[$f], $scat1file[$f], $scat2file[$f], $mapfile_pred, $mapfile_tp, $prefile[$f], $stofile[$f], $pdb2msa[$f], 
				 \$alenDCA, \@mapDCA, $target_ncnt, $N, $k, $shift, \@his, $fmax);
	  }  
	if ($pdbfile) { predictions_plot($mapfile_pred, $mapfile_tp, $pdb2msa[$f],$target_ncnt); }
    }
    elsif ($method =~ /^plmDCA$/) {
	if (!$exist_rocfile || !-e $mapfile_pred || !-e $mapfile_tp) { 
	    print "         $rocfile[$f]\n";
	    create_rocfile_plmDCA($rocfile[$f], $scat1file[$f], $scat2file[$f], $mapfile_pred, $mapfile_tp, $prefile[$f], $stofile[$f], $pdb2msa[$f], 
				  \$alenDCA, \@mapDCA, $target_ncnt, $N, $k, $shift, \@his, $fmax);
	}
	if ($pdbfile) { predictions_plot($mapfile_pred, $mapfile_tp, $pdb2msa[$f], $target_ncnt); }
    }
    elsif ($method =~ /^gremlin$/) {
	if (!$exist_rocfile || !-e $mapfile_pred || !-e $mapfile_tp) { 
	    print "         $rocfile[$f]\n";
	    create_rocfile_gremlin($rocfile[$f], $scat1file[$f], $scat2file[$f], $mapfile_pred, $mapfile_tp, $prefile[$f], $pdb2msa[$f], $target_ncnt, 
				   $N, $k, $shift, \@his, $fmax);
	}
	if ($pdbfile) { predictions_plot($mapfile_pred, $mapfile_tp, $pdb2msa[$f],$target_ncnt); }
    }
    elsif ($method =~ /^plmc$/) {
	if (!$exist_rocfile || !-e $mapfile_pred || !-e $mapfile_tp) { 
	    print "         $rocfile[$f]\n";
	    create_rocfile_plmc($rocfile[$f], $scat1file[$f], $scat2file[$f], $mapfile_pred, $mapfile_tp, $prefile[$f], $pdb2msa[$f], $target_ncnt, 
				$N, $k, $shift, \@his, $fmax);
	}
	
	if ($pdbfile) { predictions_plot($mapfile_pred, $mapfile_tp, $pdb2msa[$f], $target_ncnt); }
    }
    elsif ($method =~ /random/) {
	if (!$exist_rocfile || !-e $mapfile_pred || !-e $mapfile_tp) { 
	    print "         $rocfile[$f]\n";
	    create_rocfile_random($rocfile[$f], $scat1file[$f], $scat2file[$f], $mapfile_pred, $mapfile_tp, $prefile[$f], $pdb2msa[$f], $target_ncnt, 
				  $N, $k, $shift, \@his, $fmax);
	}
	if ($pdbfile) { predictions_plot($mapfile_pred, $mapfile_tp, $pdb2msa[$f], $target_ncnt); }
   }
    else { print "method $method not implemented yet\n"; die; }
    
    my $hfile = "$rocfile[$f].his";
    FUNCS::write_histogram($N, $k, $shift, \@his, 1, $hfile, 0);
    
    my $pdbname = ($pdbfile)? $pdb2msa[$f]->pdbname : "";
    my $stoname = ($pdbfile)? $pdb2msa[$f]->stoname : $prename[$f];
    my $maxD    = ($pdbfile)? $pdb2msa[$f]->maxD    : $maxD;
    my $which   = ($pdbfile)? $pdb2msa[$f]->which   : $which;
    my $title   = "method: $method   PDB: $pdbname   MSA: $stoname   type $which maxD: $maxD   minL: 1  fmax: $fmax";
    my $xlabel  = "distance in PDB sequence";
    my $ylabel  = "number of contacts";
    my $key     = "";
    my $psfile  = "$hfile.ps";
    my $xleft   = 1;
    my $xright  = 200;
    my $ymax    = -1;
    my $xfield  = 1;
    my $yfield  = 2;
    $seeplots   = 0;
    FUNCS::gnuplot_histo($hfile, $xfield, $yfield, $psfile, $title, $xlabel, $ylabel, $key, 0, $seeplots, $xleft, $xright, $ymax, $gnuplot);
}

my $xmax = 100;
if ($opt_x) { $xmax = $opt_x; }

my $viewplots = 0;
my $isrna = 0;
if ($opt_R) { $isrna = 1; }

my $maxpp  = 5;
if ($opt_p) { $maxpp = $opt_p; }

my $maxppv = 102;
my $maxsen = 30;
if ($opt_s) { $maxsen = $opt_s; }

my $maxF   = 40;
if ($opt_f) { $maxF = $opt_f; }


$maxpp  = 0.6;
$maxsen = 102;
$maxF   = 102;
rocplot($key, $gnuplot, $stoname, $F, \@rocfile, \@prename, $maxD, $minL, $which, $xmax, $isrna, $maxpp, $maxsen, $maxppv, $maxF, $viewplots);




####################### routines

sub  create_rocfile_rscape {
    my ($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $file, $target_ncnt, $N, $k, $shift, $his_ref, $fmax) = @_;

    open my $fp,  '>', $rocfile   || die "Can't open $rocfile: $!";
    open my $sp1, '>', $scat1file || die "Can't open $scat1file: $!";
    open my $sp2, '>', $scat2file || die "Can't open $scat2file: $!";

    my $f   = 0;
    my $f_c = 0;
    my $f_b = 0;
    my $f_w = 0;
    my $t_c = 0;
    my $t_b = 0;
    my $t_w = 0;

    my $pdblen = 0;
    my $alen = 0;
    my $type = "";
    open(FILE, "$file") || die;
    while(<FILE>) {
	if (/^\# contacts\s+(\d+)\s+\((\d+)\s+bpairs\s+(\d+)\s+wc/) {
	    $t_c = $1;
	    $t_b = $2;
	    $t_w = $3;
	}
	elsif (/^\# pdblen\s+(\S+)\s*/) {
	    $pdblen = $1;
	}
	elsif (/^\#.+alen\s+(\d+)\s+/) {
	    $alen = $1;
	}
	elsif (/^\#/) {
	}
	elsif (/^(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s*$/) {
	    $type        = $1;
	    my $i        = $2;
	    my $j        = $3;
	    my $distance = $j-$i+1; # distance in the alignment
	
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
	    if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }
	    writeline($fp, $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, ($pdblen>0)?$pdblen:$alen);
	    
	}
 	elsif (/\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s*$/) { # Not annotated as a contact		
	    $f   ++;
	    my $i        = $1;
	    my $j        = $2;
	    my $distance = $j-$i+1; # distance in the alignment
	    
	    if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }
	    writeline($fp, $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, ($pdblen>0)?$pdblen:$alen);	    
	}
    }
    close(FILE);
    close($fp);
    close($sp1);
    close($sp2);
}

sub create_rocfile_rscape_withpdb {
    my ($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $file, $pdb2msa, $target_ncnt, $N, $k, $shift, $his_ref, $fmax) = @_;

    my @revmap = @{$pdb2msa->revmap};
    open my $fp,  '>', $rocfile   || die "Can't open $rocfile: $!";
    open my $sp1, '>', $scat1file || die "Can't open $scat1file: $!";
    open my $sp2, '>', $scat2file || die "Can't open $scat2file: $!";
   
    my $f   = 0;
    my $f_c = 0;
    my $f_b = 0;
    my $f_w = 0;
    my $t_c = $pdb2msa->ncnt;
    my $t_b = $pdb2msa->nbp;
    my $t_w = $pdb2msa->nwc;
    
    my $ncnt_rscape   = 0;
    my $ncnt_rscape_f = 0;
    my @cnt_rscape;
    my @cnt_rscape_f;

    my $type = "";
    open(FILE, "$file") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/^\S*\s+(\d+)\s+(\d+)\s+(\S+)\s+\S+\s*$/) {
	    my $i          = $1;
	    my $j          = $2;
	    my $score      = $3;
	    my $distance   = $j-$i+1; # distance in the alignment
	    my $pdbi       = (@revmap)? (($revmap[$i-1]>0)? $revmap[$i-1]+1 : 0) : 0;
	    my $pdbj       = (@revmap)? (($revmap[$j-1]>0)? $revmap[$j-1]+1 : 0) : 0;
	    my $chri       = "N";
	    my $chrj       = "N";
	    my $cdistance  = -1;
	    if (@revmap && ($pdbi == 0 || $pdbj == 0)) { next; }
	    
	    if (discard_pair_by_minL($i, $j, $pdbi, $pdbj, $minL, $byali)) { 
		#print "       DISCARD    i $i $pdbi j $j $pdbj \n"; 
		next; 
	    }
	    
	    
	    $f ++;
	    printf $sp1 "$f $pdbi $pdbj %d $score\n", $pdbj-$pdbi+1;
	    if ($ncnt_rscape <= $target_ncnt) {
		$cnt_rscape[$ncnt_rscape] = CNT->new();
		$cnt_rscape[$ncnt_rscape]->{"CNT::i"}        = $pdbi;
		$cnt_rscape[$ncnt_rscape]->{"CNT::j"}        = $pdbj;
		$cnt_rscape[$ncnt_rscape]->{"CNT::posi"}     = $i;
		$cnt_rscape[$ncnt_rscape]->{"CNT::posj"}     = $j;
		$cnt_rscape[$ncnt_rscape]->{"CNT::chri"}     = $chri;
		$cnt_rscape[$ncnt_rscape]->{"CNT::chrj"}     = $chrj;
		$cnt_rscape[$ncnt_rscape]->{"CNT::bptype"}   = "CONTACT";
		$cnt_rscape[$ncnt_rscape]->{"CNT::distance"} = $cdistance;
		$ncnt_rscape ++;
	    }

	    if (PDBFUNCS::found_alicoords_in_contactlist($i, $j, $pdb2msa->{"PDB2MSA::minL"}, $byali, $pdb2msa->{"PDB2MSA::ncnt"}, 
							 \@{$pdb2msa->{"PDB2MSA::cnt"}}, \$type, \$pdbi, \$pdbj, \$chri, \$chrj, \$cdistance)) 
	    {
		if ($type ==  0) { $f_w ++; }
		if ($type <  12) { $f_b ++; }
		$f_c ++;
		#print "^^HIT $f_c/$f i $i $pdbi j $j $pdbj \n";

		if (@revmap && ($pdbi <= 0 || $pdbj <= 0)) { print "bad contact found: pdbi $pdbi pdbj $pdbj\n"; die; }
		
		printf $sp2 "$f_c $pdbi $pdbj %d $score\n", $pdbj-$pdbi+1;
		
		if ($ncnt_rscape <= $target_ncnt) {
		    $cnt_rscape_f[$ncnt_rscape_f] = CNT->new();
		    $cnt_rscape_f[$ncnt_rscape_f]->{"CNT::i"}        = $pdbi;
		    $cnt_rscape_f[$ncnt_rscape_f]->{"CNT::j"}        = $pdbj;
		    $cnt_rscape_f[$ncnt_rscape_f]->{"CNT::posi"}     = $i;
		    $cnt_rscape_f[$ncnt_rscape_f]->{"CNT::posj"}     = $j;
		    $cnt_rscape_f[$ncnt_rscape_f]->{"CNT::chri"}     = $chri;
		    $cnt_rscape_f[$ncnt_rscape_f]->{"CNT::chrj"}     = $chrj;
		    $cnt_rscape_f[$ncnt_rscape_f]->{"CNT::bptype"}   = "CONTACT";
		    $cnt_rscape_f[$ncnt_rscape_f]->{"CNT::distance"} = $cdistance;
		    $ncnt_rscape_f ++;
		}
	    }
	    #else { print "    ^^NOHIT $f_c/$f i $i $pdbi j $j $pdbj \n"; }
	    
	    if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }
	    writeline($fp,      $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdb2msa->{"PDB2MSA::pdblen"});
	}
    }
    close(FILE);
    close($fp);
    close($sp1);
    close($sp2);

    if (@revmap) { predictions_create($mapfile_pred, $mapfile_tp, $pdb2msa, $file, $ncnt_rscape, \@cnt_rscape, $ncnt_rscape_f, \@cnt_rscape_f, $maxD, $minL, $which); }
    if (@revmap) { plot_scat("R-scape", $scat1file, $scat2file); }
}




sub  create_rocfile_mfDCA {
    my ($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $prefile, $stofile, $pdb2msa, $ret_alenDCA, $mapDCA_ref, $target_ncnt, $N, $k, $shift, $his_ref, $fmax) = @_;

    my $method = "mfDCA";
    my $which  = "DI";
    if ($rocfile =~ /mfDCA\.MI/) {
	$which = "MI";
	$method .= "-$which";
    }

    my $alenDCA = $$ret_alenDCA;
    my $alen;
    if ($alenDCA < 0) {
	mapDCA2MSA($stofile, $mapDCA_ref, \$alenDCA);
    }

    parse_mfDCA($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $prefile, $pdb2msa, $mapDCA_ref, $alenDCA, $target_ncnt, $N, $k, $shift, $his_ref, $fmax, $which);
    
    $$ret_alenDCA = $alenDCA; 
}

sub  create_rocfile_plmDCA {
    my ($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $prefile, $stofile, $pdb2msa, $ret_alenDCA, $mapDCA_ref, $target_ncnt, $N, $k, $shift, $his_ref, $fmax) = @_;

    my $method = "plmDCA";
  
    my $alenDCA = $$ret_alenDCA;
    if ($alenDCA < 0) {
	mapDCA2MSA($stofile, $mapDCA_ref, \$alenDCA);
    }

    parse_plmDCA($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $prefile, $pdb2msa, $mapDCA_ref, $alenDCA, $target_ncnt, $N, $k, $shift, $his_ref, $fmax);

    $$ret_alenDCA = $alenDCA;
}

sub  create_rocfile_gremlin {
    my ($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $prefile, $pdb2msa, $target_ncnt, $N, $k, $shift, $his_ref, $fmax) = @_;

    my $method = "gremlin";
  
    parse_gremlin($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $prefile, $pdb2msa, $target_ncnt, $N, $k, $shift, $his_ref, $fmax);

}

sub  create_rocfile_plmc {
    my ($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $prefile, $pdb2msa, $target_ncnt, $N, $k, $shift, $his_ref, $fmax) = @_;

    my $method = "plmc";
  
    parse_plmc($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $prefile, $pdb2msa, $target_ncnt, $N, $k, $shift, $his_ref, $fmax);

}

sub  create_rocfile_random {
    my ($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $prefile, $pdb2msa, $target_ncnt, $N, $k, $shift, $his_ref, $fmax) = @_;
    
    my @map   = @{$pdb2msa->map};
    open my $fp,  '>', $rocfile   || die "Can't open $rocfile: $!";
    open my $sp1, '>', $scat1file || die "Can't open $scat1file: $!";
    open my $sp2, '>', $scat2file || die "Can't open $scat2file: $!";

    my $f   = 0;
    my $f_c = 0;
    my $f_b = 0;
    my $f_w = 0;
    my $t_c    = $pdb2msa->{"PDB2MSA::ncnt"};
    my $t_b    = $pdb2msa->{"PDB2MSA::nbp"};
    my $t_w    = $pdb2msa->{"PDB2MSA::nwc"};
    my $pdblen = $pdb2msa->{"PDB2MSA::pdblen"};
    my $npre   = $pdblen*($pdblen-1)/2;
    
    my $ncnt_ran   = 0;
    my $ncnt_ran_f = 0;
    my @cnt_ran;
    my @cnt_ran_f;

    my $type = "";
    while ($f < $npre && $f_c < $t_c) {
	my $i          = int(rand($pdblen-1))+1;
	my $j          = int(rand($pdblen-1))+1;
	while ($j == $i) { $j = int(rand($pdblen-1)+1); }
       	my $distance   = $j-$i+1; # distance in the pdbfile
	my $posi       = $map[$i-1]+1;
	my $posj       = $map[$j-1]+1;
	my $chri       = "N";
	my $chrj       = "N";
	my $cdistance  = -1;
	
	if ($j-$i+1 < $minL) { next; }

	$f ++;
	printf $sp1 "$f $i $j %d 0\n", $j-$i+1;

	if ($ncnt_ran < $target_ncnt) {
	    $cnt_ran[$ncnt_ran] = CNT->new();
	    $cnt_ran[$ncnt_ran]->{"CNT::i"}        = $i;
	    $cnt_ran[$ncnt_ran]->{"CNT::j"}        = $j;
	    $cnt_ran[$ncnt_ran]->{"CNT::posi"}     = $posi;
	    $cnt_ran[$ncnt_ran]->{"CNT::posj"}     = $posj;
	    $cnt_ran[$ncnt_ran]->{"CNT::chri"}     = $chri;
	    $cnt_ran[$ncnt_ran]->{"CNT::chrj"}     = $chrj;
	    $cnt_ran[$ncnt_ran]->{"CNT::bptype"}   = "CONTACT";
	    $cnt_ran[$ncnt_ran]->{"CNT::distance"} = $cdistance;
	    $ncnt_ran ++;
	}

	if (PDBFUNCS::found_pdbcoords_in_contactlist(($i<$j)?$i:$j, ($j>$i)?$j:$i, $pdb2msa->{"PDB2MSA::minL"}, $byali, $pdb2msa->{"PDB2MSA::ncnt"}, 
						     \@{$pdb2msa->{"PDB2MSA::cnt"}}, \$type, \$posi, \$posj, \$chri, \$chrj, \$cdistance)) 
	{
	    if ($type ==  0) { $f_w ++; }
	    if ($type <  12) { $f_b ++; }
	    $f_c ++;
	    
	    printf $sp2 "$f_c $i $j %d 0\n", $j-$i+1;

	    if ($ncnt_ran < $target_ncnt) {
		$cnt_ran_f[$ncnt_ran_f] = CNT->new();
		$cnt_ran_f[$ncnt_ran_f]->{"CNT::i"}        = $i;
		$cnt_ran_f[$ncnt_ran_f]->{"CNT::j"}        = $j;
		$cnt_ran_f[$ncnt_ran_f]->{"CNT::posi"}     = $posi;
		$cnt_ran_f[$ncnt_ran_f]->{"CNT::posj"}     = $posj;
		$cnt_ran_f[$ncnt_ran_f]->{"CNT::chri"}     = $chri;
		$cnt_ran_f[$ncnt_ran_f]->{"CNT::chrj"}     = $chrj;
		$cnt_ran_f[$ncnt_ran_f]->{"CNT::bptype"}   = "CONTACT";
		$cnt_ran_f[$ncnt_ran_f]->{"CNT::distance"} = $cdistance;
		$ncnt_ran_f ++;
	    }

	}
	if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }

	writeline($fp, $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdblen);
    }

    close($fp);
    close($sp1);
    close($sp2);
    
    predictions_create($mapfile_pred, $mapfile_tp, $pdb2msa, $prefile, $ncnt_ran, \@cnt_ran, $ncnt_ran_f, \@cnt_ran_f, $maxD, $minL, $which);
    if (@map) { plot_scat("random", $scat1file, $scat2file); }
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
    my ($stofile, $mapDCA_ref, $ret_alenDCA) = @_;
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
}


sub parse_mfDCA {
    my ($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $file, $pdb2msa, $mapDCA_ref, $alenDCA, $target_ncnt, $N, $k, $shift, $his_ref, $fmax, $which) = @_;

    my $sortfile = sort_mfDCA($file, $which);
    my @revmap   = @{$pdb2msa->revmap};

    open my $fp,  '>', $rocfile   || die "Can't open $rocfile: $!";
    open my $sp1, '>', $scat1file || die "Can't open $scat1file: $!";
    open my $sp2, '>', $scat2file || die "Can't open $scat2file: $!";

    my $f   = 0;
    my $f_c = 0;
    my $f_b = 0;
    my $f_w = 0;
    my $t_c = $pdb2msa->ncnt;
    my $t_b = $pdb2msa->nbp;
    my $t_w = $pdb2msa->nwc;

    my $ncnt_mfDCA   = 0;
    my $ncnt_mfDCA_f = 0;
    my @cnt_mfDCA;
    my @cnt_mfDCA_f;

    my $type;
    open(FILE, "$sortfile") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+(\S+)\s*$/) {
	    my $idca       = $1;
	    my $jdca       = $2;
	    my $score      = $3;
	    my $i          = $mapDCA_ref->[$idca];
	    my $j          = $mapDCA_ref->[$jdca];
	    my $distance   = $j-$i+1; # distance in the alignment
	    my $pdbi       = (@revmap)?(($revmap[$i-1]>0)? $revmap[$i-1]+1 : 0):0;
	    my $pdbj       = (@revmap)?(($revmap[$j-1]>0)? $revmap[$j-1]+1 : 0):0;
	    my $chri       = "N";
	    my $chrj       = "N";
	    my $cdistance  = -1;
	    if (@revmap && ($pdbi == 0 || $pdbj == 0)) { next; }

	    if (discard_pair_by_minL($i, $j, $pdbi, $pdbj, $minL, $byali)) { next; }


	    $f ++;
	    printf $sp1 "$f $pdbi $pdbj %d $score\n", $pdbj-$pdbi+1;
	    if ($ncnt_mfDCA < $target_ncnt) {
		$cnt_mfDCA[$ncnt_mfDCA] = CNT->new();
		$cnt_mfDCA[$ncnt_mfDCA]->{"CNT::i"}        = $pdbi;
		$cnt_mfDCA[$ncnt_mfDCA]->{"CNT::j"}        = $pdbj;
		$cnt_mfDCA[$ncnt_mfDCA]->{"CNT::posi"}     = $i;
		$cnt_mfDCA[$ncnt_mfDCA]->{"CNT::posj"}     = $j;
		$cnt_mfDCA[$ncnt_mfDCA]->{"CNT::chri"}     = $chri;
		$cnt_mfDCA[$ncnt_mfDCA]->{"CNT::chrj"}     = $chrj;
		$cnt_mfDCA[$ncnt_mfDCA]->{"CNT::bptype"}   = "CONTACT";
		$cnt_mfDCA[$ncnt_mfDCA]->{"CNT::distance"} = $cdistance;
		$ncnt_mfDCA ++;
	    }

	    if (PDBFUNCS::found_alicoords_in_contactlist($i, $j, $pdb2msa->minL, $byali, $pdb2msa->ncnt, \@{$pdb2msa->{"PDB2MSA::cnt"}}, 
							 \$type, \$pdbi, \$pdbj, \$chri, \$chrj, \$cdistance)) {
		if (@revmap && ($pdbi <= 0 || $pdbj <= 0)) { print "bad contact found: pdbi $pdbi pdbj $pdbj\n"; die; }
		if ($type ==  0) { $f_w ++; }
		if ($type <  12) { $f_b ++; }
		$f_c ++;
		
		printf $sp2 "$f_c $pdbi $pdbj %d $score\n", $pdbj-$pdbi+1;

		if ($ncnt_mfDCA < $target_ncnt) {
		    $cnt_mfDCA_f[$ncnt_mfDCA_f] = CNT->new();
		    $cnt_mfDCA_f[$ncnt_mfDCA_f]->{"CNT::i"}        = $pdbi;
		    $cnt_mfDCA_f[$ncnt_mfDCA_f]->{"CNT::j"}        = $pdbj;
		    $cnt_mfDCA_f[$ncnt_mfDCA_f]->{"CNT::posi"}     = $i;
		    $cnt_mfDCA_f[$ncnt_mfDCA_f]->{"CNT::posj"}     = $j;
		    $cnt_mfDCA_f[$ncnt_mfDCA_f]->{"CNT::chri"}     = $chri;
		    $cnt_mfDCA_f[$ncnt_mfDCA_f]->{"CNT::chrj"}     = $chrj;
		    $cnt_mfDCA_f[$ncnt_mfDCA_f]->{"CNT::bptype"}   = "CONTACT";
		    $cnt_mfDCA_f[$ncnt_mfDCA_f]->{"CNT::distance"} = $cdistance;
		    $ncnt_mfDCA_f ++;
		}

	    }
	    if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }
	    writeline($fp,      $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdb2msa->pdblen);
	}
    }
    close(FILE);
    close($fp);
    close($sp1);
    close($sp2);

    predictions_create($mapfile_pred, $mapfile_tp, $pdb2msa, $file, $ncnt_mfDCA, \@cnt_mfDCA, $ncnt_mfDCA_f, \@cnt_mfDCA_f, $maxD, $minL, $which);
    if (@revmap) { plot_scat("mfDCA", $scat1file, $scat2file); }
    
    system("rm $sortfile\n");
}

sub parse_plmDCA {
    my ($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $file, $pdb2msa, $mapDCA_ref, $alenDCA, $target_ncnt, $N, $k, $shift, $his_ref, $fmax) = @_;

    my $sortfile = sort_plmDCA($file);
    my @revmap   = @{$pdb2msa->revmap};

    open my $fp,  '>', $rocfile   || die "Can't open $rocfile: $!";
    open my $sp1, '>', $scat1file || die "Can't open $scat1file: $!";
    open my $sp2, '>', $scat2file || die "Can't open $scat2file: $!";

    my $f   = 0;
    my $f_c = 0;
    my $f_b = 0;
    my $f_w = 0;
    my $t_c = $pdb2msa->ncnt;
    my $t_b = $pdb2msa->nbp;
    my $t_w = $pdb2msa->nwc;

    my $ncnt_plmDCA   = 0;
    my $ncnt_plmDCA_f = 0;
    my @cnt_plmDCA;
    my @cnt_plmDCA_f;

    my $type;
    open(FILE, "$sortfile") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+(\S+)\s*$/) {
	    my $idca       = $1;
	    my $jdca       = $2;
	    my $score      = $3;
	    my $i          = $mapDCA_ref->[$idca];
	    my $j          = $mapDCA_ref->[$jdca];
	    my $distance   = $j-$i+1; # distance in the alignment
	    my $pdbi       = (@revmap)?(($revmap[$i-1]>0)? $revmap[$i-1]+1 : 0):0;
	    my $pdbj       = (@revmap)?(($revmap[$j-1]>0)? $revmap[$j-1]+1 : 0):0;
	    my $chri       = "N";
	    my $chrj       = "N";
	    my $cdistance  = -1;
	    if (@revmap && ($pdbi == 0 || $pdbj == 0)) { next; }

	    if (discard_pair_by_minL($i, $j, $pdbi, $pdbj, $minL, $byali)) { next; }

	    $f ++;
	    printf $sp1 "$f $pdbi $pdbj %d $score\n", $pdbj-$pdbi+1;

	    if ($ncnt_plmDCA < $target_ncnt) {
		$cnt_plmDCA[$ncnt_plmDCA] = CNT->new();
		$cnt_plmDCA[$ncnt_plmDCA]->{"CNT::i"}        = $pdbi;
		$cnt_plmDCA[$ncnt_plmDCA]->{"CNT::j"}        = $pdbj;
		$cnt_plmDCA[$ncnt_plmDCA]->{"CNT::posi"}     = $i;
		$cnt_plmDCA[$ncnt_plmDCA]->{"CNT::posj"}     = $j;
		$cnt_plmDCA[$ncnt_plmDCA]->{"CNT::chri"}     = $chri;
		$cnt_plmDCA[$ncnt_plmDCA]->{"CNT::chrj"}     = $chrj;
		$cnt_plmDCA[$ncnt_plmDCA]->{"CNT::bptype"}   = "CONTACT";
		$cnt_plmDCA[$ncnt_plmDCA]->{"CNT::distance"} = $cdistance;
		$ncnt_plmDCA ++;
	    }

	    if (PDBFUNCS::found_alicoords_in_contactlist($i, $j, $pdb2msa->minL, $byali, $pdb2msa->ncnt, \@{$pdb2msa->{"PDB2MSA::cnt"}}, 
							 \$type, \$pdbi, \$pdbj, \$chri, \$chrj, \$cdistance)) {
		
		if (@revmap && ($pdbi <= 0 || $pdbj <= 0)) { print "bad contact found: pdbi $pdbi pdbj $pdbj\n"; die; }
		if ($type ==  0) { $f_w ++; }
		if ($type <  12) { $f_b ++; }
		$f_c ++;
		#print "HIT $f_c/$f i $i $pdbi j $j $pdbj \n";
		
		printf $sp2 "$f_c $pdbi $pdbj %d $score\n", $pdbj-$pdbi+1;
		
		if ($ncnt_plmDCA < $target_ncnt) {
		    $cnt_plmDCA_f[$ncnt_plmDCA_f] = CNT->new();
		    $cnt_plmDCA_f[$ncnt_plmDCA_f]->{"CNT::i"}        = $pdbi;
		    $cnt_plmDCA_f[$ncnt_plmDCA_f]->{"CNT::j"}        = $pdbj;
		    $cnt_plmDCA_f[$ncnt_plmDCA_f]->{"CNT::posi"}     = $i;
		    $cnt_plmDCA_f[$ncnt_plmDCA_f]->{"CNT::posj"}     = $j;
		    $cnt_plmDCA_f[$ncnt_plmDCA_f]->{"CNT::chri"}     = $chri;
		    $cnt_plmDCA_f[$ncnt_plmDCA_f]->{"CNT::chrj"}     = $chrj;
		    $cnt_plmDCA_f[$ncnt_plmDCA_f]->{"CNT::bptype"}   = "CONTACT";
		    $cnt_plmDCA_f[$ncnt_plmDCA_f]->{"CNT::distance"} = $cdistance;
		    $ncnt_plmDCA_f ++;
		}
	    }
	    #else { print "  NOHIT $f_c/$f i $i $pdbi j $j $pdbj \n"; }

	    if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }
	    #writeline(\*STDOUT, $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdb2msa->pdblen);
	    writeline($fp,      $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdb2msa->pdblen);
	}
    }
    close(FILE);
    close($fp);
    close($sp1);
    close($sp2);
    
    predictions_create($mapfile_pred, $mapfile_tp, $pdb2msa, $file, $ncnt_plmDCA, \@cnt_plmDCA, $ncnt_plmDCA_f, \@cnt_plmDCA_f, $maxD, $minL, $which);
    if (@revmap) { plot_scat("plmDCA", $scat1file, $scat2file); }
    system("rm $sortfile\n");
}

sub parse_gremlin {
    my ($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $file, $pdb2msa, $target_ncnt, $N, $k, $shift, $his_ref, $fmax) = @_;

    my $sortfile = sort_gremlin($file);
    my @revmap   = @{$pdb2msa->revmap};
    
    open my $fp,  '>', $rocfile   || die "Can't open $rocfile: $!";
    open my $sp1, '>', $scat1file || die "Can't open $scat1file: $!";
    open my $sp2, '>', $scat2file || die "Can't open $scat2file: $!";

    my $f   = 0;
    my $f_c = 0;
    my $f_b = 0;
    my $f_w = 0;
    my $t_c = $pdb2msa->ncnt;
    my $t_b = $pdb2msa->nbp;
    my $t_w = $pdb2msa->nwc;

    my $ncnt_grem   = 0;
    my $ncnt_grem_f = 0;
    my @cnt_grem;
    my @cnt_grem_f;
    
    my $type;
    open(FILE, "$sortfile") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+(\S+)\s*$/) {
	    my $i          = $1;
	    my $j          = $2;
	    my $score      = $3;
	    my $pdbi       = (@revmap)?(($revmap[$i-1]>0)? $revmap[$i-1]+1 : 0):0;
	    my $pdbj       = (@revmap)?(($revmap[$j-1]>0)? $revmap[$j-1]+1 : 0):0;
	    my $distance   = $j-$i+1; # distance in the alignment
	    my $chri       = "N";
	    my $chrj       = "N";
	    my $cdistance  = -1;
	    if (@revmap && ($pdbi == 0 || $pdbj == 0)) { next; }
	    
	    if (discard_pair_by_minL($i, $j, $pdbi, $pdbj, $minL, $byali)) { next; }

	    $f ++;
	    printf $sp1 "$f $pdbi $pdbj %d $score\n", $pdbj-$pdbi+1;

	    if ($ncnt_grem < $target_ncnt) {
		$cnt_grem[$ncnt_grem] = CNT->new();
		$cnt_grem[$ncnt_grem]->{"CNT::i"}        = $pdbi;
		$cnt_grem[$ncnt_grem]->{"CNT::j"}        = $pdbj;
		$cnt_grem[$ncnt_grem]->{"CNT::posi"}     = $i;
		$cnt_grem[$ncnt_grem]->{"CNT::posj"}     = $j;
		$cnt_grem[$ncnt_grem]->{"CNT::chri"}     = $chri;
		$cnt_grem[$ncnt_grem]->{"CNT::chrj"}     = $chrj;
		$cnt_grem[$ncnt_grem]->{"CNT::bptype"}   = "CONTACT";
		$cnt_grem[$ncnt_grem]->{"CNT::distance"} = $cdistance;
		$ncnt_grem ++;
	    }

	    if (PDBFUNCS::found_alicoords_in_contactlist($i, $j, $pdb2msa->minL, $byali, $pdb2msa->ncnt, \@{$pdb2msa->{"PDB2MSA::cnt"}}, 
							 \$type, \$pdbi, \$pdbj, \$chri, \$chrj, \$cdistance)) 
	    {
		if (@revmap && ($pdbi <= 0 || $pdbj <= 0)) { print "bad contact found: pdbi $pdbi pdbj $pdbj\n"; die; }
		if ($type ==  0) { $f_w ++; }
		if ($type <  12) { $f_b ++; }
		$f_c ++;

		printf $sp2 "$f_c $pdbi $pdbj %d $score\n", $pdbj-$pdbi+1;

		if ($ncnt_grem < $target_ncnt) {
		    $cnt_grem_f[$ncnt_grem_f] = CNT->new();
		    $cnt_grem_f[$ncnt_grem_f]->{"CNT::i"}        = $pdbi;
		    $cnt_grem_f[$ncnt_grem_f]->{"CNT::j"}        = $pdbj;
		    $cnt_grem_f[$ncnt_grem_f]->{"CNT::posi"}     = $i;
		    $cnt_grem_f[$ncnt_grem_f]->{"CNT::posj"}     = $j;
		    $cnt_grem_f[$ncnt_grem_f]->{"CNT::chri"}     = $chri;
		    $cnt_grem_f[$ncnt_grem_f]->{"CNT::chrj"}     = $chrj;
		    $cnt_grem_f[$ncnt_grem_f]->{"CNT::bptype"}   = "CONTACT";
		    $cnt_grem_f[$ncnt_grem_f]->{"CNT::distance"} = $cdistance;
		    $ncnt_grem_f ++;
		}
	    }
	    
	    if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }
	    writeline($fp,      $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdb2msa->pdblen);
	}
    }
    close(FILE);
    close($fp);
    close($sp1);
    close($sp2);

    predictions_create($mapfile_pred, $mapfile_tp, $pdb2msa, $file, $ncnt_grem, \@cnt_grem, $ncnt_grem_f, \@cnt_grem_f, $maxD, $minL, $which);
    if (@revmap) { plot_scat("gremlin", $scat1file, $scat2file); }
    system("rm $sortfile\n");
}

sub parse_plmc {
    my ($rocfile, $scat1file, $scat2file, $mapfile_pred, $mapfile_tp, $file, $pdb2msa, $target_ncnt, $N, $k, $shift, $his_ref, $fmax) = @_;

    my $sortfile = sort_plmc($file);
    my @revmap   = @{$pdb2msa->revmap};

    open my $fp,  '>', $rocfile   || die "Can't open $rocfile: $!";
    open my $sp1, '>', $scat1file || die "Can't open $scat1file: $!";
    open my $sp2, '>', $scat2file || die "Can't open $scat2file: $!";

    my $f   = 0;
    my $f_c = 0;
    my $f_b = 0;
    my $f_w = 0;
    my $t_c = $pdb2msa->ncnt;
    my $t_b = $pdb2msa->nbp;
    my $t_w = $pdb2msa->nwc;

    my $ncnt_plmc   = 0;
    my $ncnt_plmc_f = 0;
    my @cnt_plmc;
    my @cnt_plmc_f;
    
    my $type;
    open(FILE, "$sortfile") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+(\S+)\s*$/) {
	    my $i          = $1;
	    my $j          = $2;
	    my $score      = $3;
	    my $pdbi       = (@revmap)?(($revmap[$i-1]>0)? $revmap[$i-1]+1 : 0):0;
	    my $pdbj       = (@revmap)?(($revmap[$j-1]>0)? $revmap[$j-1]+1 : 0):0;
	    my $distance   = $j-$i+1; # distance in the alignment
	    my $chri       = "N";
	    my $chrj       = "N";
	    my $cdistance  = -1;
	    
	    if (@revmap && ($pdbi == 0 || $pdbj == 0)) { next; }
	    
	    if (discard_pair_by_minL($i, $j, $pdbi, $pdbj, $minL, $byali)) { next; }

	    $f ++;
	    printf $sp1 "$f $pdbi $pdbj %d $score\n", $pdbj-$pdbi+1;

	    if ($ncnt_plmc < $target_ncnt) {
		$cnt_plmc[$ncnt_plmc] = CNT->new();
		$cnt_plmc[$ncnt_plmc]->{"CNT::i"}        = $pdbi;
		$cnt_plmc[$ncnt_plmc]->{"CNT::j"}        = $pdbj;
		$cnt_plmc[$ncnt_plmc]->{"CNT::posi"}     = $i;
		$cnt_plmc[$ncnt_plmc]->{"CNT::posj"}     = $j;
		$cnt_plmc[$ncnt_plmc]->{"CNT::chri"}     = $chri;
		$cnt_plmc[$ncnt_plmc]->{"CNT::chrj"}     = $chrj;
		$cnt_plmc[$ncnt_plmc]->{"CNT::bptype"}   = "CONTACT";
		$cnt_plmc[$ncnt_plmc]->{"CNT::distance"} = $cdistance;
		$ncnt_plmc ++;
	    }

	    if (PDBFUNCS::found_alicoords_in_contactlist($i, $j, $pdb2msa->minL, $byali, $pdb2msa->ncnt, \@{$pdb2msa->{"PDB2MSA::cnt"}}, 
							 \$type, \$pdbi, \$pdbj, \$chri, \$chrj, \$cdistance)) 
	    {
		if (@revmap && ($pdbi <= 0 || $pdbj <= 0)) { print "bad contact found: pdbi $pdbi pdbj $pdbj\n"; die; }
		if ($type ==  0) { $f_w ++; }
		if ($type <  12) { $f_b ++; }
		$f_c ++;

		printf $sp2 "$f_c $pdbi $pdbj %d $score\n", $pdbj-$pdbi+1;

		if ($ncnt_plmc < $target_ncnt) {
		    $cnt_plmc_f[$ncnt_plmc_f] = CNT->new();
		    $cnt_plmc_f[$ncnt_plmc_f]->{"CNT::i"}        = $pdbi;
		    $cnt_plmc_f[$ncnt_plmc_f]->{"CNT::j"}        = $pdbj;
		    $cnt_plmc_f[$ncnt_plmc_f]->{"CNT::posi"}     = $i;
		    $cnt_plmc_f[$ncnt_plmc_f]->{"CNT::posj"}     = $j;
		    $cnt_plmc_f[$ncnt_plmc_f]->{"CNT::chri"}     = $chri;
		    $cnt_plmc_f[$ncnt_plmc_f]->{"CNT::chrj"}     = $chrj;
		    $cnt_plmc_f[$ncnt_plmc_f]->{"CNT::bptype"}   = "CONTACT";
		    $cnt_plmc_f[$ncnt_plmc_f]->{"CNT::distance"} = $cdistance;
		    $ncnt_plmc_f ++;
		}
	    }
	    
	    if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }
	    writeline($fp,      $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdb2msa->pdblen);
	}
    }
    close(FILE);
    close($fp);
    close($sp1);
    close($sp2);

    predictions_create($mapfile_pred, $mapfile_tp, $pdb2msa, $file, $ncnt_plmc, \@cnt_plmc, $ncnt_plmc_f, \@cnt_plmc_f, $maxD, $minL, $which);
    if (@revmap) { plot_scat("plmclin", $scat1file, $scat2file); }
    system("rm $sortfile\n");
}

sub predictions_create {
    my ($mapfile_pred, $mapfile_tp, $pdb2msa, $file, $ncnt_pred, $cnt_pred_ref, $ncnt_tp, $cnt_tp_ref, $maxD, $minL, $which) = @_;

    open(MAPPRED,  ">$mapfile_pred")  || die;
    PDBFUNCS::contactlist_print(\*MAPPRED, $ncnt_pred, $cnt_pred_ref, 1);
    close(MAPPRED);
    
    open(MAPTP,  ">$mapfile_tp")  || die;
    PDBFUNCS::contactlist_print(\*MAPTP, $ncnt_tp, $cnt_tp_ref, 1);
    close(MAPTP);    
}

sub predictions_plot {
    my ($mapfile_pred, $mapfile_tp, $pdb2msa, $target_ncnt) = @_;

    if (!-e $mapfile_pred || !-e $mapfile_tp) { return; }
    
    my $minpdbx = 1e+50;
    my $maxpdbx = 0;
    my $ncnt_pred = PDBFUNCS::contactlistfile_parse($mapfile_pred, \$minpdbx, \$maxpdbx);
    my $ncnt_tp   = PDBFUNCS::contactlistfile_parse($mapfile_tp,   \$minpdbx, \$maxpdbx);
    
    my $frac_ncnt = int($ncnt_tp/$target_ncnt*1000)/10;
    my $xfield  = 1;
    my $yfield  = 4;
    my $xylabel = "PDB position";
    my $title   = "First $target_factor*L($target_ncnt) predictions -- $ncnt_tp($frac_ncnt%) correct";
    print "         $mapfile_tp\n";
    print "         $title\n";
    
    my $nf = 3;
    my @mapfile;
    $mapfile[0] = $mapfile_tp;
    $mapfile[1] = $mapfile_pred;
    $mapfile[2] = $pdb2msa->mapfile;
    PDBFUNCS::plot_contact_map($nf, \@mapfile, $minpdbx, $maxpdbx, $xfield, $yfield, $title, $xylabel, $gnuplot, 0);
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

    if (0) { system("more $sortfile\n"); }
    
    return $sortfile;
}

sub sort_plmc {
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
	if (/^(\d+)\s\-\s(\d+)\s\-\s\S+\s+(\S+)\s*$/) {
	    $i{$n}  = $1;  # this plmc file uses notation: 1....
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

    if (0) { system("more $sortfile\n"); }
    
    return $sortfile;
}

sub rocplot {
    my ($key, $gnuplot, $stoname, $F, $file_ref, $prename_ref, $maxD, $minL, $which, $xmax, $isrna, $maxpp, $maxsen, $maxppv, $maxF, $seeplots) = @_;


   my $psfile = "results/$key-$stoname.N$F.maxD$maxD.minL$minL.type$which.ps";
    
    #if ($psfile =~ /\/([^\/]+)\s*$/) { $psfile = "$1"; }
    my $pdffile = $psfile;
    if ($pdffile =~ /^(\S+).ps$/) { $pdffile = "$1.pdf"; }
    print "\n rocFILE: $psfile\n";

    my $xlabel;
    my $ylabel;
    my $title  = "$stoname";
    my $x;
    my $y;
    my $x_max, my $x_min;
    my $y_max, my $y_min;
    
    open(my $gp, '|'."gnuplot") || die "Gnuplot: $!";
 
    print $gp "set terminal postscript color solid 14\n";
    print $gp "set output '$psfile'\n";

    print $gp "set style line 1   lt 1 lc rgb 'black'   pt 1 ps 0.5 lw 1\n";
    print $gp "set style line 2   lt 1 lc rgb 'brown'   pt 1 ps 0.5 lw 1\n";
    print $gp "set style line 3   lt 1 lc rgb 'grey'    pt 1 ps 0.5 lw 1\n";
    print $gp "set style line 4   lt 1 lc rgb 'purple'  pt 1 ps 0.5 lw 1\n";
    print $gp "set style line 5   lt 1 lc rgb 'orange'  pt 1 ps 0.5 lw 1\n";
    print $gp "set style line 6   lt 1 lc rgb 'blue'    pt 1 ps 0.5 lw 1\n";
    print $gp "set style line 7   lt 1 lc rgb 'cyan'    pt 1 ps 0.5 lw 1\n";
    #print $gp "set style line 7   lt 1 ls rgb '#1F78B4' pt 1 ps 0.5 lw 1\n"; # dark blue
    print $gp "set style line 8   lt 1 lc rgb '#005A32' pt 1 ps 0.5 lw 1\n"; # dark green
    print $gp "set style line 9   lt 1 lc rgb '#74C476' pt 1 ps 0.5 lw 1\n"; # light green
    print $gp "set style line 10  lt 1 lc rgb 'red'     pt 1 ps 0.5 lw 1\n";
 
    #print $gp "set style line 1 lc rgb '#084594' pt 65 ps 0.5 lw 3\n"; # very blue
    #print $gp "set style line 2 lc rgb '#1F78B4' pt 65 ps 0.5 lw 3\n"; # dark blue
    #print $gp "set style line 3 lc rgb '#F16913' pt 65 ps 0.5 lw 3\n"; # orange
    #print $gp "set style line 4 lc rgb '#005A32' pt 65 ps 0.5 lw 3\n"; # dark green
    #print $gp "set style line 5 lc rgb '#74C476' pt 65 ps 0.5 lw 3\n"; # light green
    #print $gp "set style line 6 lc rgb '#4A1486' pt 65 ps 0.5 lw 3\n"; # dark purple
    #print $gp "set style line 7 lc rgb '#BCBDDC' pt 65 ps 0.5 lw 3\n"; # light purple
    #print $gp "set style line 8 lc rgb 'red'     pt 65 ps 0.5 lw 3\n"; # red

    my $logscale = 0;
    $xlabel = "number of predictions per position";
    $ylabel = "PPV contacts (%)";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = $maxppv;
    $x = 9;
    $y = 17;
    #roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "SEN contacts";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = $maxsen;
    $x = 9;
    $y = 16;
    #roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "F contacts (%)";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = $maxF;
    $x = 9;
    $y = 18;
    #roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);

    # basepairs
    if ($isrna) {
	$xlabel = "number of predictions per position";
	$ylabel = "PPV bpairs (%)";
	$x_min = 0.001;
	$x_max = $maxpp;
	$y_min = 0;
	$y_max = $maxppv;
	$x = 9;
	$y = 20;
	#roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
	$xlabel = "number of predictions per position";
	$ylabel = "SEN bpairs (%)";
	$x_min = 0.001;
	$x_max = $maxpp;
	$y_min = 0;
	$y_max = $maxsen;
	$x = 9;
	$y = 19;
	#roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
	$xlabel = "number of predictions per position";
	$ylabel = "F bpairs";
	$x_min = 0.001;
	$x_max = $maxpp;
	$y_min = 0;
	$y_max = $maxF;
	$x = 9;
	$y = 21;
	#roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
	
	$xlabel = "number of predictions per position";
	$ylabel = "PPV WC (%)";
	$x_min = 0.001;
	$x_max = $maxpp;
	$y_min = 0;
	$y_max = $maxppv;
	$x = 9;
	$y = 23;
	#roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
	$xlabel = "number of predictions per position";
	$ylabel = "SEN WC (%)";
	$x_min = 0.001;
	$x_max = $maxpp;
	$y_min = 0;
	$y_max = $maxsen;
	$x = 9;
	$y = 22;
	#roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
	$xlabel = "number of predictions per position";
	$ylabel = "F WC";
	$x_min = 0.001;
	$x_max = $maxpp;
	$y_min = 0;
	$y_max = $maxF;
	$x = 9;
	$y = 24;
	#roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    }
    
    $logscale = 0;
    $xlabel = "number of predictions";
    $ylabel = "PPV contacts (%)";
    $x_min = 1;
    $x_max = $xmax;
    $y_min = 0;
    $y_max = $maxppv;
    $x = 1;
    $y = 17;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions";
    $ylabel = "SEN contacts (%)";
    $x_min = 1;
    $x_max = $xmax;
    $y_min = 0;
    $y_max = $maxsen;
    #$y_max = 450;
    $x = 1;
    $y = 16;
    #$y = 2;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions";
    $ylabel = "F contacts (%)";
    $x_min = 1;
    $x_max = $xmax;
    $y_min = 0;
    $y_max = $maxF;
    $x = 1;
    $y = 18;
    #roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    
    # basepairs
    if ($isrna) {
	$logscale = 0;
	$xlabel = "number of predictions";
	$ylabel = "PPV bpairs (%)";
	$x_min = 1;
	$x_max = $xmax;
	$y_min = 0;
	$y_max = $maxppv;
	$x = 1;
	$y = 20;
	roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
	$xlabel = "number of predictions";
	$ylabel = "SEN bpairs (%)";
	$x_min = 1;
	$x_max = $xmax;
	$y_min = 0;
	$y_max = $maxsen;
	#$y_max = 450;
	$x = 1;
	$y = 19;
	#$y = 2;
	roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
	$xlabel = "number of predictions";
	$ylabel = "F contacts (%)";
	$x_min = 1;
	$x_max = $xmax;
	$y_min = 0;
	$y_max = $maxF;
	$x = 1;
	$y = 21;
	#roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);

	$logscale = 0;
	$xlabel = "number of predictions";
	$ylabel = "PPV WC (%)";
	$x_min = 1;
	$x_max = $xmax;
	$y_min = 0;
	$y_max = $maxppv;
	$x = 1;
	$y = 23;
	roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
	$xlabel = "number of predictions";
	$ylabel = "SEN WC (%)";
	$x_min = 1;
	$x_max = $xmax;
	$y_min = 0;
	$y_max = $maxsen;
	#$y_max = 450;
	$x = 1;
	$y = 22;
	#$y = 2;
	roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
	$xlabel = "number of predictions";
	$ylabel = "F WC (%)";
	$x_min = 1;
	$x_max = $xmax;
	$y_min = 0;
	$y_max = $maxF;
	$x = 1;
	$y = 24;
	#roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    }
    
    $logscale = 0;
    $xlabel = "SEN contacts (%)";
    $ylabel = "PPV contacts (%)";
    $x_min = 0;
    $x_max = $maxsen;
    $y_min = 0;
    $y_max = $maxppv;
    $x = 16;
    $y = 17;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    if ($isrna) {
	$xlabel = "SEN bpairs (%)";
	$ylabel = "PPV bpairs (%)";
	$x_min = 0;
	$x_max = $maxsen;
	$y_min = 0;
	$y_max = $maxppv;
	$x = 19;
	$y = 20;
	roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
	$xlabel = "SEN WC (%)";
	$ylabel = "PPV WC (%)";
	$x_min = 0;
	$x_max = $maxsen;
	$y_min = 0;
	$y_max = $maxppv;
	$x = 22;
	$y = 23;
	roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    }
    
    close($gp);

    if ($seeplots) { system ("open $psfile&\n"); }

    
}


sub roc_oneplot {
    my ($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $xmin, $xmax, $ymin, $ymax, $logscale, $nolines) = @_;
   
    my $cmd = "";
    my $m = 1;
    
    print $gp "set title  '$title'\n";
    print $gp "set xlabel '$xlabel'\n";
    print $gp "set ylabel '$ylabel'\n";
    print $gp "set xrange [$xmin:$xmax]\n";
    print $gp "set yrange [$ymin:$ymax]\n";
    if ($logscale) { print $gp "set logscale x\n"; }
    for (my $f = 0; $f < $F; $f++) {
	my $key = $prename_ref->[$f];
	if ($nolines) {
	    $cmd .= ($f == $F-1)? "'$file_ref->[$f]' using $x:$y  title '$key'            ls $m"   : "'$file_ref->[$f]' using $x:$y  title '$key'            ls $m, ";
	}
	else {
	    $cmd .= ($f == $F-1)? "'$file_ref->[$f]' using $x:$y  title ''                ls $m, " : "'$file_ref->[$f]' using $x:$y  title ''                ls $m, ";
	    $cmd .= ($f == $F-1)? "'$file_ref->[$f]' using $x:$y  title '$key' with lines ls $m"   : "'$file_ref->[$f]' using $x:$y  title '$key' with lines ls $m, ";
	}
	$m ++;
	if ($m == 11) { $m = 1; }

    }
    print $gp "plot $cmd\n";
    if ($logscale) { print $gp "unset logscale\n"; }
}

sub writeline {
    my ($fp, $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdblen) = @_;

    my $sen_c, my $sen_b, my $sen_w;
    my $ppv_c, my $ppv_b, my $ppv_w;
    my $F_c,   my $F_b,   my $F_w;
    
    FUNCS::calculateF($f_c, $t_c, $f, \$sen_c, \$ppv_c, \$F_c);
    FUNCS::calculateF($f_b, $t_b, $f, \$sen_b, \$ppv_b, \$F_b);
    FUNCS::calculateF($f_w, $t_w, $f, \$sen_w, \$ppv_w, \$F_w);

    # tab separated fields
    # ---------------------
    #
    # f         fc         fb         fw         tc         tb         tw         pdblen
    # 1         2          3          4          5          6          7          8
    #
    # f/pdblen  fc/pdblen  fb/pdblen  fw/pdblen  tc/pdblen  tb/pdblen  tw/pdblen 
    # 9         10         11         12         13         14         15  
    #
    # sen_c  ppv_c  F_c
    # 16     17     18
    #
    # sen_b  ppv_b  F_b
    # 19     20     21
    #
    # sen_w  ppv_w  F_w
    # 22     23     24
    #
    if ($pdblen > 0) {
	printf $fp "$f\t$f_c\t$f_b\t$f_w\t$t_c\t$t_b\t$t_w\t$pdblen\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
	$f/$pdblen, $f_c/$pdblen, $f_b/$pdblen, $f_w/$pdblen, $t_c/$pdblen, $t_b/$pdblen, $t_w/$pdblen, $sen_c, $ppv_c, $F_c, $sen_b, $ppv_b, $F_b, $sen_w, $ppv_w, $F_w;
    }
    else {
	printf $fp "$f\t$f_c\t$f_b\t$f_w\t$t_c\t$t_b\t$t_w\t$pdblen\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
	0, 0, 0, 0, 0, 0, 0, $sen_c, $ppv_c, $F_c, $sen_b, $ppv_b, $F_b, $sen_w, $ppv_w, $F_w;
    }
}



sub discard_pair_by_minL {
    my ($i, $j, $pdbi, $pdbj, $minL, $byaly) = @_;

    if ($byaly) {  if ($j-$i+1       < $minL) { return 1; } }
    else        {  if ($pdbj-$pdbi+1 < $minL) { return 1; } }
    return 0;	
}

sub plot_scat {
    my ($method, $scat1file, $scat2file) = @_;

    my $psfile = "$scat1file.ps";
    
    print "          $psfile\n";
   
    my $xlabel = "backbone distance";
    my $ylabel = "score";
    my $title  = "$method";
    my $x;
    my $y;
    open(my $gp, '|'."gnuplot") || die "Gnuplot: $!";
 
    print $gp "set terminal postscript color solid 14\n";
    print $gp "set output '$psfile'\n";
    FUNCS::gnuplot_define_styles ($gp);
    
    my $cmd = "";
    
    print $gp "set title  '$title'\n";
    print $gp "set xlabel '$xlabel'\n";
    print $gp "set ylabel '$ylabel'\n";
    #print $gp "set xrange [$xmin:$xmax]\n";
    #print $gp "set yrange [$ymin:$ymax]\n";
    $cmd = "'$scat1file' using 4:5  title 'F'            ls 5, ";
    $cmd .= "'$scat2file' using 4:5  title 'TP'          ls 4";
    print $gp "plot $cmd\n";
    close($gp);   
}


sub structure_from_msa {
    my ($stofile, $pdb2msa_ref, $maxD, $minL) = @_;

    $$pdb2msa_ref = PDB2MSA->new();
    my $stoname = $stofile;
    if ($stoname =~ /(PF[^\.]+)\./) { $stoname = $1; }
    if ($stoname =~ /(RF[^\.]+)\./) { $stoname = $1; }
    if ($stoname =~ /([^\.]+)\./)   { $stoname = $1; }
    
    my $nsq;
    my @name;
    my @sq;
    my @ct;
    my $ss;
    my $rfasq;
    my $msalen = FUNCS::parse_stofile($stofile, \$nsq, \@name, \@sq, \$ss, \@ct, \$rfasq);
    my $ncnt = 0;
    my @cnt;

    for (my $x = 0; $x < $msalen; $x ++) {
	if ($ct[$x] >= 0 && $ct[$x] > $x) {
	    
	    my $posi = $x+1;
	    my $posj = $ct[$x]+1;
	    if ($posj - $posi + 1 >= $minL) {
		$cnt[$ncnt] = CNT->new();
		$cnt[$ncnt]->{"CNT::i"}        = -1;
		$cnt[$ncnt]->{"CNT::j"}        = -1;
		$cnt[$ncnt]->{"CNT::posi"}     = $posi;
		$cnt[$ncnt]->{"CNT::posj"}     = $posj;
		$cnt[$ncnt]->{"CNT::chri"}     = substr($sq[0], $x,      1);
		$cnt[$ncnt]->{"CNT::chrj"}     = substr($sq[0], $ct[$x], 1);
		$cnt[$ncnt]->{"CNT::bptype"}   = "WWc";
		$cnt[$ncnt]->{"CNT::distance"} = -1;
		$ncnt ++;
	    }
	}
    }
    
    $$pdb2msa_ref->{"PDB2MSA::pdbname"}   = "";
    $$pdb2msa_ref->{"PDB2MSA::stoname"}   = $stoname;
    $$pdb2msa_ref->{"PDB2MSA::pdblen"}    = -1;
    $$pdb2msa_ref->{"PDB2MSA::msalen"}    = $msalen;
    
   
    $$pdb2msa_ref->{"PDB2MSA::ncnt"}      = $ncnt;
    $$pdb2msa_ref->{"PDB2MSA::nbp"}       = $ncnt;
    $$pdb2msa_ref->{"PDB2MSA::nwc"}       = $ncnt;
    $$pdb2msa_ref->{"PDB2MSA::maxD"}      = $maxD;
    $$pdb2msa_ref->{"PDB2MSA::minL"}      = $minL;
    $$pdb2msa_ref->{"PDB2MSA::byali"}     = 0;
    $$pdb2msa_ref->{"PDB2MSA::which"}     = 0;
    @{$$pdb2msa_ref->{"PDB2MSA::cnt"}}    = @cnt;

}


sub structure_from_ecfile {
    my ($ecfile, $pdb2msa_ref, $maxD, $minL) = @_;

    $$pdb2msa_ref = PDB2MSA->new();
    my $ecname = $ecfile;
    if ($ecname =~ /([^\.]+)\./) { $ecname = $1; }

    my $ncnt = 0;
    my $nbp  = 0;
    my $nwc  = 0;
    my @cnt;
    open(FILE, "$ecfile") || die;
    while(<FILE>) {
	if (/^(\S+)\s+(\d+)\s+(\d+)\s+\d+\s+(\S)\s+\d+\s+(\S)\s+\d+\s+(\S+)\s+(\S+)\s*$/) {
	    my $sc       = $1;
	    my $i        = $2;
	    my $j        = $3;
	    my $chri     = $4;
	    my $chrj     = $5;
	    my $distance = $6;
	    my $type     = $7;

	    if ($j-$i+1 >= $minL && $distance <= $maxD) {
		if    ($type =~ /^None$/) { next; }
		elsif ($type =~ /^cWW$/)  { $type = "WWc"; $nwc ++; $nbp ++ }
		else                                              { $nbp ++ }
		
		$cnt[$ncnt] = CNT->new();
		$cnt[$ncnt]->{"CNT::i"}        = -1;
		$cnt[$ncnt]->{"CNT::j"}        = -1;
		$cnt[$ncnt]->{"CNT::posi"}     = $i+1; # this file uses notation: 0... 
		$cnt[$ncnt]->{"CNT::posj"}     = $j+1;
		$cnt[$ncnt]->{"CNT::chri"}     = $chri;
		$cnt[$ncnt]->{"CNT::chrj"}     = $chrj;
		$cnt[$ncnt]->{"CNT::bptype"}   = $type;
		$cnt[$ncnt]->{"CNT::distance"} = $distance;
		$ncnt ++;
	    }
	}
    }
    close(FILE);
	
    $$pdb2msa_ref->{"PDB2MSA::pdbname"}   = "";
    $$pdb2msa_ref->{"PDB2MSA::stoname"}   = $ecname;
    $$pdb2msa_ref->{"PDB2MSA::pdblen"}    = -1;
    $$pdb2msa_ref->{"PDB2MSA::msalen"}    = -1;
    
   
    $$pdb2msa_ref->{"PDB2MSA::ncnt"}      = $ncnt;
    $$pdb2msa_ref->{"PDB2MSA::nbp"}       = $nbp;
    $$pdb2msa_ref->{"PDB2MSA::nwc"}       = $nwc;
    $$pdb2msa_ref->{"PDB2MSA::maxD"}      = $maxD;
    $$pdb2msa_ref->{"PDB2MSA::minL"}      = $minL;
    $$pdb2msa_ref->{"PDB2MSA::byali"}     = -1;
    $$pdb2msa_ref->{"PDB2MSA::which"}     = 0;
    @{$$pdb2msa_ref->{"PDB2MSA::cnt"}}    = @cnt;

}

sub structure_from_contactmapfile {
    my ($cmapfile, $pdb2msa_ref, $target_maxD, $target_minL) = @_;

    $$pdb2msa_ref = PDB2MSA->new();
    my $cmapname = $cmapfile;
    if ($cmapname =~ /([^\.]+)\./) { $cmapname = $1; }

    my $ncnt = 0;
    my $nbp  = 0;
    my $nwc  = 0;
    my @cnt;
    my $alen;
    my $pdblen;
    open(FILE, "$cmapfile") || die;
    while(<FILE>) {
	if (/^\#\s+(\d+)\s+(\d+)\s+\|\s+bptype\s+(\S+)\s*$/) {
	    my $i        = $1;
	    my $j        = $2;
	    my $type     = $3;

	    if    ($type =~ /^CONTACT$/) { }
	    elsif ($type =~ /^WWc$/)     { $nbp ++; $nwc ++ }
	    else                         { $nbp ++; }
	    
	    $cnt[$ncnt] = CNT->new();
	    $cnt[$ncnt]->{"CNT::i"}        = -1;
	    $cnt[$ncnt]->{"CNT::j"}        = -1;
	    $cnt[$ncnt]->{"CNT::posi"}     = $i; # this file uses notation: 1..len 
	    $cnt[$ncnt]->{"CNT::posj"}     = $j;
	    $cnt[$ncnt]->{"CNT::chri"}     = '';
	    $cnt[$ncnt]->{"CNT::chrj"}     = '';
	    $cnt[$ncnt]->{"CNT::bptype"}   = $type;
	    $cnt[$ncnt]->{"CNT::distance"} = -1;
	    $ncnt ++;	    
	}
	elsif (/^\# maxD\s+(\S+)/) {
	    my $maxD = $1;
	    if ($target_maxD != $maxD) { print "change maxD from $maxD to $target_maxD in cmapfile $cmapfile\n"; die; }
	}
	elsif (/^\# mind\s+(\S+)/) {
	    my $minL = $1;
	    if ($target_minL != $minL) { print "change minL from $minL to $target_minL in cmapfile $cmapfile\n"; die; }
	}
	elsif (/^\# alen\s+(\S+)/) {
	    $alen = $1;
	}
	elsif (/^\# pdblen\s+(\S+)/) {
	    $pdblen = $1;
	}
    }
    close(FILE);
	
    $$pdb2msa_ref->{"PDB2MSA::pdbname"}   = "";
    $$pdb2msa_ref->{"PDB2MSA::stoname"}   = $cmapname;
    $$pdb2msa_ref->{"PDB2MSA::pdblen"}    = $pdblen;
    $$pdb2msa_ref->{"PDB2MSA::msalen"}    = $alen;
       
    $$pdb2msa_ref->{"PDB2MSA::ncnt"}      = $ncnt;
    $$pdb2msa_ref->{"PDB2MSA::nbp"}       = $nbp;
    $$pdb2msa_ref->{"PDB2MSA::nwc"}       = $nwc;
    $$pdb2msa_ref->{"PDB2MSA::maxD"}      = $maxD;
    $$pdb2msa_ref->{"PDB2MSA::minL"}      = $minL;
    $$pdb2msa_ref->{"PDB2MSA::byali"}     = -1;
    $$pdb2msa_ref->{"PDB2MSA::which"}     = 0;
    @{$$pdb2msa_ref->{"PDB2MSA::cnt"}}    = @cnt;

}
