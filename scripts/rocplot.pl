#!/usr/bin/perl -w
# rocplot.pl

use strict;
use Class::Struct;

# find directory where the script is installed
use FindBin;
use lib $FindBin::Bin;
use PDBFUNCS;
use FUNCS;

use vars qw ($opt_D $opt_G $opt_L $opt_P $opt_r $opt_R $opt_W $opt_v);  # required if strict used
use Getopt::Std;
getopts ('D:G:L:P:rR:W:v');


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
my $stoname = $prefile[0];
if ($stoname =~ /(PF[^\.]+)\./) { $stoname = $1; }
if ($stoname =~ /(RF[^\.]+)\./) { $stoname = $1; }

my $rscapebin = shift;
my $gnuplot   = shift;
if ($gnuplot) { use constant GNUPLOT => '$gnuplot'; }


my $currdir = $ENV{PWD};

my $maxD = 8;
if ($opt_D) { $maxD = $opt_D; }

my $minL = 1;
if ($opt_L) { $minL = $opt_L; }
for (my $f = 0; $f < $F; $f ++) {  $rocfile[$f] = ($minL>0)? "$prename[$f].minL$minL.roc":"$prename[$f].roc"; }

my $dornaview = 0;
if ($opt_r) { $dornaview = 1; }

my $which = "MIN"; #options: CA C MIN AVG NOH / C1' (for RNA suggested by Westhof)
if ($opt_W) { $which = "$opt_W"; }
my $seeplots = 0;

my $verbose = 0;
if ($opt_v) { $verbose = 1; }

# Map the  pdb2msa structure to the input alignment
my $pdbfile = "";
my $pdb2msa;
if ($opt_P) { 
    $pdbfile = "$opt_P"; 
    PDBFUNCS::pdb2msa($gnuplot, $rscapebin, $pdbfile, $stofile, \$pdb2msa, $maxD, $minL, $which, $dornaview, $seeplots);
}

# Map the pdb2msa structure to the stofile used by gremlin
my $stofile_gremlin = "";
my $pdb2msa_gremlin;
if ($opt_G) { 
    $stofile_gremlin = "$opt_G"; 
    PDBFUNCS::pdb2msa($gnuplot, $rscapebin, $pdbfile, $stofile_gremlin, \$pdb2msa_gremlin, $maxD, $minL, $which, $dornaview, $seeplots);
}  

# Map the pdb2msa structure to the stofile used by rscape (if different from the main one)
my $stofile_rscape = "";
my $pdb2msa_rscape;
if ($opt_R) { 
    $stofile_rscape = "$opt_R"; 
    PDBFUNCS::pdb2msa($gnuplot, $rscapebin, $pdbfile, $stofile_rscape, \$pdb2msa_rscape, $maxD, $minL, $which, $dornaview, $seeplots);
}  

# add a random file
my $dorandom = 1;
if ($dorandom) {
    $prename[$F] = "$stoname.random";
    $prefile[$F] = "results/random/$stoname.random";
    $rocfile[$F] = "results/random/$stoname.random.minL$minL.roc";
    $F ++;
}

my $alenDCA = -1;
my @mapDCA;
for (my $f = 0; $f < $F; $f ++) {
    my$method = "";
    
    if    ($prefile[$f] =~ /results\/(\S+)_filtered\//) { $method = $1; }
    elsif ($prefile[$f] =~ /results\/([^\/]+)\//)       { $method = $1; }

    my @his;
    my $N = 1000;
    my $k = 1;
    my $shift = 0;
    my $fmax = 100;
    FUNCS::init_histo_array($N, $k, \@his);

    print "\n$method: $prefile[$f]\n";
    if ($method =~ /^R-scape$/) {
	## both functions below should produce exactly the same results (only if the minL used running R-scape is the same than here)
	if ($pdbfile =~ //) {
	    create_rocfile_rscape($rocfile[$f], $prefile[$f], $N, $k, $shift, \@his, $fmax);
	    $dorandom = 0;
	}
	else { create_rocfile_rscape_withpdb($rocfile[$f], $prefile[$f], ($stofile_rscape)?$pdb2msa_rscape:$pdb2msa,     $N, $k, $shift, \@his, $fmax); }
    }
    elsif ($method =~ /^mfDCA$/) {
	if ($pdbfile) { create_rocfile_mfDCA  ($rocfile[$f], $prefile[$f], $stofile, $pdb2msa, \$alenDCA, \@mapDCA,      $N, $k, $shift, \@his, $fmax); }
    }
    elsif ($method =~ /^plmDCA$/) {
	if ($pdbfile) { create_rocfile_plmDCA ($rocfile[$f], $prefile[$f], $stofile, $pdb2msa, \$alenDCA, \@mapDCA,      $N, $k, $shift, \@his, $fmax); }
    }
    elsif ($method =~ /^gremlin$/) {
	if ($pdbfile) { create_rocfile_gremlin($rocfile[$f], $prefile[$f], ($stofile_gremlin)?$pdb2msa_gremlin:$pdb2msa, $N, $k, $shift, \@his, $fmax); }
    }
    elsif ($method =~ /random/) {
	create_rocfile_random($rocfile[$f], $pdb2msa, $N, $k, $shift, \@his, $fmax);
    }
    else { print "method $method not implemented yet\n"; die; }
    
    my $hfile = "$rocfile[$f].his";
    FUNCS::write_histogram($N, $k, $shift, \@his, 1, $hfile, 0);
    
    my $pdbname = $pdb2msa->pdbname;
    my $stoname = $pdb2msa->stoname;
    my $maxD    = $pdb2msa->maxD;
    my $title  = "method: $method   PDB: $pdbname   MSA: $stoname   maxD: $maxD   minL: 1  fmax: $fmax";
    my $xlabel = "distance in PDB sequence";
    my $ylabel = "number of contacts";
    my $key    = "";
    my $psfile = "$hfile.ps";
    my $xleft  = 1;
    my $xright = 200;
    my $ymax   = -1;
    my $xfield = 1;
    my $yfield = 2;
    $seeplots  = 0;
    FUNCS::gnuplot_histo($hfile, $xfield, $yfield, $psfile, $title, $xlabel, $ylabel, $key, 0, $seeplots, $xleft, $xright, $ymax);
}

my $xmax = 1000;
my $viewplots = 0;
rocplot($stoname, $F, \@rocfile, \@prename, $xmax, $viewplots);






####################### routines

sub  create_rocfile_rscape {
    my ($rocfile, $file, $N, $k, $shift, $his_ref, $fmax) = @_;

    open my $fp, '>', $rocfile || die "Can't open $rocfile: $!";

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
	    writeline($fp, $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdblen);
	    
	}
 	elsif (/\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s*$/) { # Not annotated as a contact		
	    $f   ++;
	    my $i        = $1;
	    my $j        = $2;
	    my $distance = $j-$i+1; # distance in the alignment
	    
	    if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }
	    writeline($fp, $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdblen);
	    
	}
    }
    close(FILE);
    close($fp);
}

sub create_rocfile_rscape_withpdb {
    my ($rocfile, $file, $pdb2msa, $N, $k, $shift, $his_ref, $fmax) = @_;
    
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
	elsif (/^\S*\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s*$/) {
	    my $i        = $1;
	    my $j        = $2;
	    my $distance = $j-$i+1; # distance in the alignment
	    my $pdbi        = -1;
	    my $pdbj        = -1;
	    if (PDBFUNCS::found_alicoords_in_contactlist($i, $j, $pdb2msa->{"PDB2MSA::minL"}, $pdb2msa->{"PDB2MSA::ncnt"}, 
					       \@{$pdb2msa->{"PDB2MSA::cnt"}}, \$type, \$pdbi, \$pdbj)) 
	    {
		if    ($type ==  0) { $f_w ++; }
		elsif ($type <  12) { $f_b ++; }
		$f_c ++;
		if ($pdbi <= 0 || $pdbj <= 0) { print "bad contact found: pdbi $pdbi pdbj $pdbj\n"; die; }
	    }
	    $f ++;
	    if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }
	    writeline($fp,      $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdb2msa->{"PDB2MSA::pdblen"});
	}
    }
    close(FILE);
    close($fp);
}




sub  create_rocfile_mfDCA {
    my ($rocfile, $prefile, $stofile, $pdb2msa, $ret_alenDCA, $mapDCA_ref, $N, $k, $shift, $his_ref, $fmax) = @_;

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

    parse_mfDCA($rocfile, $prefile, $pdb2msa, $mapDCA_ref, $alenDCA, $N, $k, $shift, $his_ref, $fmax, $which);
    
    $$ret_alenDCA = $alenDCA; 
}

sub  create_rocfile_plmDCA {
    my ($rocfile, $prefile, $stofile, $pdb2msa, $ret_alenDCA, $mapDCA_ref, $N, $k, $shift, $his_ref, $fmax) = @_;

    my $method = "plmDCA";
  
    my $alenDCA = $$ret_alenDCA;
    if ($alenDCA < 0) {
	mapDCA2MSA($stofile, $mapDCA_ref, \$alenDCA);
    }

    parse_plmDCA($rocfile, $prefile, $pdb2msa, $mapDCA_ref, $alenDCA, $N, $k, $shift, $his_ref, $fmax);

    $$ret_alenDCA = $alenDCA;
}

sub  create_rocfile_gremlin {
    my ($rocfile, $prefile, $pdb2msa, $N, $k, $shift, $his_ref, $fmax) = @_;

    my $method = "gremlin";
  
    parse_gremlin($rocfile, $prefile, $pdb2msa, $N, $k, $shift, $his_ref, $fmax);

}

sub  create_rocfile_random {
    my ($rocfile, $pdb2msa, $N, $k, $shift, $his_ref, $fmax) = @_;

    open my $fp, '>', $rocfile || die "Can't open $rocfile: $!";

    my $f   = 0;
    my $f_c = 0;
    my $f_b = 0;
    my $f_w = 0;
    my $t_c    = $pdb2msa->{"PDB2MSA::ncnt"};
    my $t_b    = $pdb2msa->{"PDB2MSA::nbp"};
    my $t_w    = $pdb2msa->{"PDB2MSA::nwc"};
    my $pdblen = $pdb2msa->{"PDB2MSA::pdblen"};
    my $npre   = $pdblen*($pdblen-1)/2;
    
    my $type = "";
    while ($f < $npre) {
	my $i = int(rand($pdblen));
	my $j = int(rand($pdblen));
	while ($j == $i) { $j = int(rand($pdblen)); }
       
	my $distance = $j-$i+1; # distance in the pdbfile
	my $posi;
	my $posj;
	if (PDBFUNCS::found_pdbcoords_in_contactlist(($i<$j)?$i:$j, ($j>$i)?$j:$i, $pdb2msa->{"PDB2MSA::minL"}, $pdb2msa->{"PDB2MSA::ncnt"}, 
						     \@{$pdb2msa->{"PDB2MSA::cnt"}}, \$type, \$posi, \$posj)) 
	{
	    if    ($type ==  0) { $f_w ++; }
	    elsif ($type <  12) { $f_b ++; }
	    $f_c ++;
	}
	$f ++;
	if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }

	writeline($fp,      $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdblen);
    }

    close($fp);
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
    my ($rocfile, $file, $pdb2msa, $mapDCA_ref, $alenDCA, $N, $k, $shift, $his_ref, $fmax, $which) = @_;

    my $sortfile = sort_mfDCA($file, $which);
    open my $fp, '>', $rocfile || die "Can't open $rocfile: $!";

    my $f   = 0;
    my $f_c = 0;
    my $f_b = 0;
    my $f_w = 0;
    my $t_c = $pdb2msa->ncnt;
    my $t_b = $pdb2msa->nbp;
    my $t_w = $pdb2msa->nwc;

    my $type;
    open(FILE, "$sortfile") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+\S+\s*$/) {
	    my $idca = $1;
	    my $jdca = $2;
	    my $i        = $mapDCA_ref->[$idca];
	    my $j        = $mapDCA_ref->[$jdca];
	    my $distance = $j-$i+1; # distance in the alignment
	    my $pdbi        = -1;
	    my $pdbj        = -1;
	    if (PDBFUNCS::found_alicoords_in_contactlist($i, $j, $pdb2msa->minL, $pdb2msa->ncnt, \@{$pdb2msa->{"PDB2MSA::cnt"}}, \$type, \$pdbi, \$pdbj)) {
		if ($pdbi <= 0 || $pdbj <= 0) { print "bad contact found: pdbi $pdbi pdbj $pdbj\n"; die; }
		if    ($type ==  0) { $f_w ++; }
		elsif ($type <  12) { $f_b ++; }
		$f_c ++;
	    }
	    $f ++;
	    if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }
	    writeline($fp,      $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdb2msa->pdblen);
	}
    }
    close(FILE);
    close($fp);

    system("rm $sortfile\n");
}

sub parse_plmDCA {
    my ($rocfile, $file, $pdb2msa, $mapDCA_ref, $alenDCA, $N, $k, $shift, $his_ref, $fmax) = @_;

    my $sortfile = sort_plmDCA($file);
    open my $fp, '>', $rocfile || die "Can't open $rocfile: $!";

    my $f   = 0;
    my $f_c = 0;
    my $f_b = 0;
    my $f_w = 0;
    my $t_c = $pdb2msa->ncnt;
    my $t_b = $pdb2msa->nbp;
    my $t_w = $pdb2msa->nwc;

    my $type;
    open(FILE, "$sortfile") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+\S+\s*$/) {
	    my $idca     = $1;
	    my $jdca     = $2;
	    my $i        = $mapDCA_ref->[$idca];
	    my $j        = $mapDCA_ref->[$jdca];
	    my $distance = $j-$i+1; # distance in the alignment
	    my $pdbi     = -1;
	    my $pdbj     = -1;
	    if (PDBFUNCS::found_alicoords_in_contactlist($i, $j, $pdb2msa->minL, $pdb2msa->ncnt, \@{$pdb2msa->{"PDB2MSA::cnt"}}, \$type, \$pdbi, \$pdbj)) {
		if ($pdbi <= 0 || $pdbj <= 0) { print "bad contact found: pdbi $pdbi pdbj $pdbj\n"; die; }
		if    ($type ==  0) { $f_w ++; }
		elsif ($type <  12) { $f_b ++; }
		$f_c ++;
	    }
	    $f ++;
	    if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }
	    #writeline(\*STDOUT, $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdb2msa->pdblen);
	    writeline($fp,      $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdb2msa->pdblen);
	}
    }
    close(FILE);
    close($fp);

    system("rm $sortfile\n");
}

sub parse_gremlin {
    my ($rocfile, $file, $pdb2msa, $N, $k, $shift, $his_ref, $fmax) = @_;

    my $sortfile = sort_gremlin($file);
    
    open my $fp, '>', $rocfile || die "Can't open $rocfile: $!";
    my $f   = 0;
    my $f_c = 0;
    my $f_b = 0;
    my $f_w = 0;
    my $t_c = $pdb2msa->ncnt;
    my $t_b = $pdb2msa->nbp;
    my $t_w = $pdb2msa->nwc;

    my $type;
    open(FILE, "$sortfile") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+\S+\s*$/) {
	    my $i = $1;
	    my $j = $2;
	    my $pdbi        = -1;
	    my $pdbj        = -1;
	    my $distance = $j-$i+1; # distance in the alignment
	    if (PDBFUNCS::found_alicoords_in_contactlist($i, $j, $pdb2msa->minL, $pdb2msa->ncnt, \@{$pdb2msa->{"PDB2MSA::cnt"}}, \$type, \$pdbi, \$pdbj)) {
		if ($pdbi <= 0 || $pdbj <= 0) { print "bad contact found: pdbi $pdbi pdbj $pdbj\n"; die; }
		if    ($type ==  0) { $f_w ++; }
		elsif ($type <  12) { $f_b ++; }
		$f_c ++;
	    }
	    $f ++;
	    if ($f <= $fmax) { FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, $his_ref); }
	    writeline($fp,      $f, $f_c, $f_b, $f_w, $t_c, $t_b, $t_w, $pdb2msa->pdblen);
	}
    }
    close(FILE);
    close($fp);

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

sub rocplot {
    my ($stoname, $F, $file_ref, $prename_ref, $xmax, $seeplots) = @_;


   my $psfile = "results/$stoname.N$F.ps";
    
    #if ($psfile =~ /\/([^\/]+)\s*$/) { $psfile = "$1"; }
    my $pdffile = $psfile;
    if ($pdffile =~ /^(\S+).ps$/) { $pdffile = "$1.pdf"; }
    print "FILE: $psfile\n";

    my $maxpp = 2;

    my $xlabel;
    my $ylabel;
    my $title  = "$stoname";
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
    $xlabel = "SEN contacts (%)";
    $ylabel = "PPV contacts (%)";
    $x_min = 0;
    $x_max = 50;
    $y_min = 0;
    $y_max = 100;
    $x = 16;
    $y = 17;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "SEN bpairs (%)";
    $ylabel = "PPV bpairs (%)";
    $x_min = 0;
    $x_max = 50;
    $y_min = 0;
    $y_max = 100;
    $x = 19;
    $y = 20;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);

    $logscale = 0;
    $xlabel = "number of predictions per position";
    $ylabel = "PPV contacts (%)";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 9;
    $y = 17;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "SEN contacts";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 9;
    $y = 16;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "F contacts (%)";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 9;
    $y = 18;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);

    # basepairs
    $xlabel = "number of predictions per position";
    $ylabel = "PPV bpairs (%)";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 9;
    $y = 20;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "SEN bpairs (%)";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 9;
    $y = 19;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions per position";
    $ylabel = "F bpairs";
    $x_min = 0.001;
    $x_max = $maxpp;
    $y_min = 0;
    $y_max = 100;
    $x = 9;
    $y = 21;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
 
    $logscale = 1;
    $xlabel = "number of predictions";
    $ylabel = "PPV contacts (%)";
    $x_min = 1;
    $x_max = $xmax;
    $y_min = 0;
    $y_max = 100;
    $x = 1;
    $y = 17;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions";
    $ylabel = "SEN contacts (%)";
    $x_min = 1;
    $x_max = $xmax;
    $y_min = 0;
    $y_max = 100;
    $x = 1;
    $y = 16;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);
    $xlabel = "number of predictions";
    $ylabel = "F contacts (%)";
    $x_min = 1;
    $x_max = $xmax;
    $y_min = 0;
    $y_max = 100;
    $x = 1;
    $y = 18;
    roc_oneplot($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $x_min, $x_max, $y_min, $y_max, $logscale);


    close($gp);

    #system ("/usr/local/bin/ps2pdf $psfile $pdffile\n"); 
    #system("rm $psfile\n");
    if ($seeplots) { system ("open $psfile&\n"); }

    
}


sub roc_oneplot {
    my ($gp, $F, $file_ref, $prename_ref, $x, $y, $xlabel, $ylabel, $title, $xmin, $xmax, $ymin, $ymax, $logscale, $nolines) = @_;
   
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
	if ($nolines) {
	    $cmd .= ($f == $F-1)? "'$file_ref->[$f]' using $x:$y  title '$key'            ls $m"   : "'$file_ref->[$f]' using $x:$y  title '$key'            ls $m, ";
	}
	else {
	    $cmd .= ($f == $F-1)? "'$file_ref->[$f]' using $x:$y  title ''                ls $m, " : "'$file_ref->[$f]' using $x:$y  title ''                ls $m, ";
	    $cmd .= ($f == $F-1)? "'$file_ref->[$f]' using $x:$y  title '$key' with lines ls $m"   : "'$file_ref->[$f]' using $x:$y  title '$key' with lines ls $m, ";
	}
	$m ++; if ($m == 7) { $m = 8; }
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
    printf $fp "$f\t$f_c\t$f_b\t$f_w\t$t_c\t$t_b\t$t_w\t$pdblen\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
    $f/$pdblen, $f_c/$pdblen, $f_b/$pdblen, $f_w/$pdblen, $t_c/$pdblen, $t_b/$pdblen, $t_w/$pdblen, $sen_c, $ppv_c, $F_c, $sen_b, $ppv_b, $F_b, $sen_w, $ppv_w, $F_w;
}
