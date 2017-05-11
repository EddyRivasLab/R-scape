#!/usr/bin/perl -w
#compare_to_contactmap.pl

use strict;
use Class::Struct;

use vars qw ($opt_D $opt_L $opt_N $opt_R $opt_v);  # required if strict used
use Getopt::Std;
getopts ('D:L:N:Rv');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  compare_to_contactmap.pl [options]  N <file1>..<fileN> <afafile> <filtermapfile>\n\n";
        print "options:\n";
 	exit;
}

my $F = shift;

my @prefile;
my @prename;
my @plotfile;
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
my $afafile       = shift;
my $filtermapfile = shift;
my $pfamname = $prefile[0];
if ($pfamname =~ /(PF[^\.]+)\./) { $pfamname = $1; }
my $oafafile = "data/$pfamname.afa";
print "original file: $oafafile\n";

my $gnuplot = shift;
use constant GNUPLOT => '$gnuplot';


my $NCH = 0;
my @cmpfile;
my @corfile;
my $maxD = 8;
if ($opt_D) { $maxD = $opt_D; }

my $cmpfileA = "data/$pfamname.chainA.maxD$maxD.map";
my $corfileA = "data/$pfamname.chainA.maxD$maxD.cor";
my $cmpfileB = "data/$pfamname.chainB.maxD$maxD.map";
my $corfileB = "data/$pfamname.chainB.maxD$maxD.cor";
if (-f $cmpfileA && -f $corfileA) {
    $cmpfile[$NCH] = "$cmpfileA";
    $corfile[$NCH] = "$corfileA";
    $NCH ++;
}
else {
    $cmpfileA = "data/$pfamname.chainX.maxD8.map";
    $corfileA = "data/$pfamname.chainX.maxD8.cor";
    if (-f $cmpfileA && -f $corfileA) {
	$cmpfile[$NCH] = "$cmpfileA";
	$corfile[$NCH] = "$corfileA";
	$NCH ++;
    }
}
if (-f $cmpfileB && -f $corfileB) {
    $cmpfile[$NCH] = "$cmpfileB";
    $corfile[$NCH] = "$corfileB";
    $NCH ++;
}
print "NCH $NCH\n";
for (my $ch = 0; $ch < $NCH; $ch ++) {
    print "$cmpfile[$ch]\n";
}

my $dorandom = 1;

my $seeplots = 0;
my $verbose = 0;

my $min_pre = 1;
my $max_pre = 10000;
my $inc_pre = 40;
my $plot_pre;
my $maxplot_pre = $min_pre;
if ($opt_N) { $max_pre = $opt_N; }
my $minL = 1;
if ($opt_L) { $minL = $opt_L; }

# fraction of positions with >= minC contacts
my $minC = 3;

my $replot = 0;
if ($opt_R) { $replot = 1; } 
for (my $f = 0; $f < $F; $f ++) { $plotfile[$f] = "$prename[$f].cmp.minL$minL.plot"; }

for (my $f = 0; $f < $F; $f ++) {
    my $method = "";
    if ($prefile[$f] =~ /results\/(\S+)\//) { $method = $1; }
    if ($method =~ /^mfDCA/) {
	for (my $ff = $F; $ff > $f; $ff --) {
	    $plotfile[$ff+1] = $plotfile[$ff];
	    $prefile[$ff+1]  = $prefile[$ff];
	}
	$plotfile[$f+1] = "$prename[$f].MI.cmp.minL$minL.plot";
	$prefile[$f+1]  = "$prename[$f]";
	$F ++;
	$f ++;	
    }
}
# add a random file
if ($dorandom) {
    if ($prename[0] =~ /results\/[^\/]+_([^\/]+)\//) { 
	$plotfile[$F] = "results/random_$1/$pfamname.random.cmp.minL$minL.plot";
    }
    else { 
	$plotfile[$F] = "results/random/$pfamname.random.cmp.minL$minL.plot";
    }
}

if (!$replot) { for (my $f = 0; $f < $F; $f ++) { system("rm $plotfile[$f]\n"); }  if ($dorandom) { system("rm $plotfile[$F]\n"); } }

my $ncmp = 0;
my @cmp_i;
my @cmp_j;

if (!$replot) {
    for (my $ch = 0; $ch < $NCH; $ch ++) { parse_cmp($cmpfile[$ch], $minL, \$ncmp, \@cmp_i, \@cmp_j); }

    # map the coordenates in DCA output files to those of the input alignment
    #
    # mlDCA and plmDCA trim internally the alignmen by removing coluns with . or lower case
    # we need mapDCA() because mlDCA and plmDCA report coordinates relative
    # to this filtered alignment, not relative to the input alignment
    #
    my $alen     = 0;    # length of the input alignment
    my $alen_dca = 0;    # length of the columns used by DCA (no . no lower case)
    my @mapDCA;      # mapDCA[1..L] valued in 1..alen
    parse_afa($afafile, \@mapDCA, \$alen, \$alen_dca);
    $plot_pre = 2*$alen;

    # if data/filter_gremlin need to use the .map file
    # the .map file maps the input alignment to the original alignment
    #
    my $oalen;   # length of the original alignment
    my @foo;
    parse_afa($oafafile, \@foo, \$oalen,);
    my @map2original;
    my @original2map;
    parse_filtermap($filtermapfile, \@map2original, \@original2map, $alen, $oalen);
    
    # map the coordenates in the pdb file to those of the input alignment
    #
    # we need this to create the random predictions with |i-j| >= minL
    # because the pdf "cor" file is relative to the original alignment
    # we need to use @original2map()
    my $lenPDB;
    my @mapPDB;  # mapPDB[0..lenPDB-1] valued in 0..alen-1
    parse_cor($corfile[0], $minL, \@original2map, \$lenPDB, \@mapPDB);
  
    my $method;
    my $L = $alen;
    for (my $f = 0; $f < $F; $f ++) {
	$method = "";
	if    ($prefile[$f] =~ /results\/(\S+)_filtered\//) { $method = $1; }
	elsif ($prefile[$f] =~ /results\/(\S+)\//)          { $method = $1; }

	my $which = "DI";
	if ($method =~ /^mfDCA/ && $plotfile[$f] =~ /mfDCA\.MI/) {
	    $which = "MI";
	    $method .= "-$which";
	}
	if ($method =~ /DCA/)     { $L = $alen_dca; }
	if ($method =~ /R-scape/) { parse_rscape_alen($prefile[$f], \$L); }
	
	print "\n$method: $prefile[$f]\n";
	create_plotfile($min_pre, $max_pre, $inc_pre, $plot_pre, $plotfile[$f], $prefile[$f], $method, $L, $minL, $ncmp, \@cmp_i, \@cmp_j, 
			\@map2original, \@mapDCA, $which, $oalen);
    }
    if ($dorandom) {
	$method = "random";
	print "\n$method:\n";
	create_random_plotfile($min_pre, $max_pre, $inc_pre, $plot_pre, $plotfile[$F], $method, $L, $minL, $ncmp, \@cmp_i, \@cmp_j, 
			       \@map2original, $lenPDB, \@mapPDB, $oalen);
    }
}

if ($replot) { $maxplot_pre = $max_pre; }
$maxplot_pre = 2000;
$seeplots = 0;
plot_roc("$outname", $minL, ($dorandom)?$F+1:$F, \@plotfile, $maxplot_pre, $seeplots);

sub create_random_plotfile {
    my ($min_pre, $max_pre, $inc_pre, $plot_pre, $plotfile, $method, $L, $minL, $ncmp, $cmp_i_ref, $cmp_j_ref, $map2o_ref, $lenPDB, $mapPDB_ref, $oalen) = @_;
 
    my 	$pre_cutoff = $min_pre;
    my  $maxran = 1000;
    my $nmax = $lenPDB*($lenPDB-1)/2;
    if ($nmax > $maxran) { $nmax = $maxran; }
    print "NMAX $nmax\n";

    my $npre = 0;
    my @i;
    my @j;
    my @pre_i;
    my @pre_j;
    my $pr_npre;
    while ($pre_cutoff <= $nmax) {
	if    ($pre_cutoff <= 50)   { $inc_pre = 10;  }
	elsif ($pre_cutoff <= 1000) { $inc_pre = 100; }
	elsif ($pre_cutoff >  1000) { $inc_pre = 700; }
	
	my $npre = random_onecutoff($plotfile, $method, $L, $minL, $ncmp, $cmp_i_ref, $cmp_j_ref, $map2o_ref, $lenPDB, $mapPDB_ref, $pre_cutoff, $pr_npre, 
				    $plot_pre, \$npre, \@i, \@j, \@pre_i, \@pre_j, $oalen);

	if ($pre_cutoff == $nmax) { last; }
	if ($pre_cutoff < $nmax && $pre_cutoff+$inc_pre > $nmax) { $pre_cutoff  = $nmax;    }
	else                                                     { $pre_cutoff += $inc_pre; }
	$pr_npre = $npre;
    }
}

sub create_plotfile {
    my ($min_pre, $max_pre, $inc_pre, $plot_pre, $plotfile, $prefile, $method, $L, $minL, $ncmp, $cmp_i_ref, $cmp_j_ref, 
	$map2o_ref, $mapDCA_ref, $whichDCA, $oalen) = @_;

    my @pre_i;
    my @pre_j;
    my $npre = onefile_onecutoff($prefile, $method, $minL, $max_pre, \@pre_i, \@pre_j, $whichDCA, $oalen);
    if ($npre < $max_pre) { $max_pre = $npre; }
    
    my $pre_cutoff = $min_pre;
    my $pr_npre = 0;
    while ($pre_cutoff <= $max_pre) {
	if    ($pre_cutoff <= 50)   { $inc_pre = 10;   }
	elsif ($pre_cutoff <= 1000) { $inc_pre = 100;  }
	elsif ($pre_cutoff >  1000) { $inc_pre = 700; }
	$npre = write_stats_to_file($plotfile, $method, $L, $ncmp, $cmp_i_ref, $cmp_j_ref, $map2o_ref, $mapDCA_ref, \@pre_i, \@pre_j, 
				    $pre_cutoff, $pr_npre, $plot_pre, $oalen);

	if ($npre < $pre_cutoff) { if ($npre > $maxplot_pre) { $maxplot_pre = $npre; } last; }
	$pre_cutoff += $inc_pre;
	$pr_npre = $npre;
    }
}

 
sub random_onecutoff {
    my ($plotfile, $method, $L, $minL, $ncmp, $cmp_i_ref, $cmp_j_ref, $map2o_ref, $lenPDB, $mapPDB_ref, $pre_cutoff, $pr_npre, 
	$plot_pre, $ret_npre, $i_ref, $j_ref, $pre_i_ref, $pre_j_ref, $oalen) = @_;

    my $npre = $$ret_npre;
    my $maxit = 500000;
    my $nit = 0;

    # generate $pre_cutoff random pairs at a minL = $minL
    while ($nit < $maxit && $npre < $pre_cutoff) {
	#pick i 0...$lenPDB-1
	my $i = int(rand()*$lenPDB);
	while ($mapPDB_ref->[$i] < 0) {
	    $i = int(rand()*$lenPDB);
	}
	#pick j, such that |j-i| >= minL
	my $j = int(rand()*$lenPDB);
	while ($nit < $maxit && ($mapPDB_ref->[$j] < 0 || abs($j-$i) < $minL || !pair_is_new($i, $j, $npre, $i_ref, $j_ref)) ) {
	    $j = int(rand()*$lenPDB);
	    $nit ++;
	}
	$i_ref->[$npre] = ($i<$j)? $i:$j;
	$j_ref->[$npre] = ($i<$j)? $j:$i;
	$pre_i_ref->[$npre] = $mapPDB_ref->[$i_ref->[$npre]]+1;
	$pre_j_ref->[$npre] = $mapPDB_ref->[$j_ref->[$npre]]+1;

	$npre ++;
	$nit ++;
    }
    print "random: $npre\n";
    write_stats_to_file($plotfile, $method, $L, $ncmp, $cmp_i_ref, $cmp_j_ref, $map2o_ref, "", $pre_i_ref, $pre_j_ref, $npre, $pr_npre, $plot_pre, $oalen);

    $$ret_npre = $npre;
    return $npre;
}

sub pair_is_new {
    my ($i, $j, $npre, $i_ref, $j_ref) = @_;

    for (my $n = 0; $n < $npre; $n ++) {
	if ($i == $i_ref->[$n] && $j == $j_ref->[$n]) { return 0; }
	if ($j == $i_ref->[$n] && $i == $j_ref->[$n]) { return 0; }
    }
    return 1;
}

sub onefile_onecutoff{
    my ($prefile, $method, $minL, $pre_cutoff, $pre_i_ref, $pre_j_ref, $whichDCA) = @_;
    
    my $npre = 0;

    if    ($method =~ /^R-scape/) { parse_rscape  ($prefile, $minL, \$npre,  $pre_i_ref,  $pre_j_ref, $pre_cutoff); }
    elsif ($method =~ /^mfDCA/)   { parse_mfdca   ($prefile, $minL, \$npre,  $pre_i_ref,  $pre_j_ref, $pre_cutoff, $whichDCA);}
    elsif ($method =~ /^plmDCA/)  { parse_plmdca  ($prefile, $minL, \$npre,  $pre_i_ref,  $pre_j_ref, $pre_cutoff); }
    elsif ($method =~ /^gremlin/) { parse_gremplin($prefile, $minL, \$npre,  $pre_i_ref,  $pre_j_ref, $pre_cutoff); }
    else                          { print "method? $method\n"; die; }
 
    return $npre;
}

sub plot_map {
    my ($key, $mapfile_tf, $mapfile_fp, $cmpfile) = @_;
    
    my $psfile = "$mapfile_tf.ps";
    #if ($psfile =~ /\/([^\/]+)\s*$/) { $psfile = "$1"; }
    my $pdffile = $psfile;
    if ($pdffile =~ /^(\S+).ps$/) { $pdffile = "$1.pdf"; }

    my $xlabel = "alignment position";
    my $title  = "$mapfile_tf";
    
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    
    print GP "set terminal postscript color solid 14\n";
    FUNCS::gnuplot_define_styles (*GP);

    print GP "set output '$psfile'\n";
    print GP "unset key\n";
    print GP "set size ratio -1\n";
    print GP "set size 1,1\n";

    print GP "set nokey\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$xlabel'\n";
    #print GP "set xrange [0:$len]\n";
    #print GP "set yrange [0:$len]\n";

    print GP "set title \"$title\\n\\n$key\"\n";
    #print GP "set title '$title'\n";

    my $m = 1;
    my $cmd  = "'$cmpfile' using 1:3  title '' ls $m, '$cmpfile' using 3:1  title '' ls $m, "; $m ++;
    $cmd    .= "'$mapfile_tf' using 1:2  title '' ls $m, "; $m +=2;
    $cmd    .= "'$mapfile_fp' using 1:2  title '' ls $m"; 	
    
    print    "plot $cmd\n";
    print GP "plot $cmd\n";
    close (GP);

    system ("ps2pdf $psfile $pdffile\n"); 
    system("rm $psfile\n");
    if ($seeplots) { system ("open $pdffile&\n"); }
}

sub  write_stats_to_file {
    my ($plotfile, $method, $L, $ncmp, $cmp_i_ref, $cmp_j_ref, $map2o_ref, $mapDCA_ref, $pre_i_ref, $pre_j_ref, $pre_cutoff, $pr_npre, $plot_pre, $oalen) = @_;

    my $mapfile_tf = "$plotfile.npre$pre_cutoff.tf.map";
    my $mapfile_fp = "$plotfile.npre$pre_cutoff.fp.map";
    open (TF, ">$mapfile_tf") || die;
    open (FP, ">$mapfile_fp") || die;

    my @contact;
    for (my $l = 0; $l <= $oalen; $l ++) {
	$contact[$l] = 0;
    }
    
    my $tf = 0;
    my $fp = 0;
    my $c;
    my $npre = 0;
    for (my $p = 0; $p < $pre_cutoff; $p ++) {
	my $ip = $map2o_ref->[$pre_i_ref->[$p]];
	my $jp = $map2o_ref->[$pre_j_ref->[$p]];
	if ($method =~ /^mfDCA/ || $method =~ /^plmDCA/) {
	    $ip = $map2o_ref->[$mapDCA_ref->[$pre_i_ref->[$p]]];
	    $jp = $map2o_ref->[$mapDCA_ref->[$pre_j_ref->[$p]]];
	}
	
	$npre ++;
	
	for ($c = 0; $c < $ncmp; $c ++) {
	    my $ic = $cmp_i_ref->[$c];
	    my $jc = $cmp_j_ref->[$c];
	    
	    if ($ip == $ic && $jp == $jc) {
		$tf ++;

		$contact[$ip] ++;
		$contact[$jp] ++;

		#print    "$ip $jp | $contact[$ip] $contact[$jp] \n";
		print TF "$ip $jp\n";
		last;
	    }
	}
	if ($c == $ncmp) { $fp ++; print FP "$ip $jp\n"; }
    }
    close(TF);
    close(FP);
    
    my $sen;
    my $ppv;
    my $F;
    FUNCS::calculateF($tf, $ncmp, $npre, \$sen, \$ppv, \$F);

    my $nc = 0;
    for (my $l = 0; $l <= $oalen; $l ++) {
	if ($contact[$l] >= $minC) { $nc ++; }
    }
    
    open(PLOT, ">>$plotfile") || die;
    if ($npre<=100) { printf "true_found=$tf found=$npre true=$ncmp | sen=%.2f%% ppv=%.2f%% F=%.2f%% | >= $minC contacts %.2f\n", $sen, $ppv, $F, ($L>0)?$nc/$L:0; }
    printf PLOT "$npre $tf $npre $ncmp %.2f %.2f %.2f %f %f %d %d\n", $sen, $ppv, $F, ($L>0)?$npre/$L:0, ($L>0)?$nc/$L:0, $nc, $L;
    close(PLOT);

    $sen = int($sen*1000)/1000;
    $ppv = int($ppv*1000)/1000;
    $F   = int($F*1000)  /1000;
    if ($pr_npre < $plot_pre && $npre >= $plot_pre) {
	my $title = "TF=$tf FP=$fp F=$npre T=$ncmp sen=$sen ppv=$ppv F=$F";
	print "$title\n";
	plot_map($title, $mapfile_tf, $mapfile_fp, $cmpfile[0]); 
    }
    return $npre;
}
   

sub parse_rscape {
    my ($file, $minL, $ret_npre, $pre_i_ref, $pre_j_ref, $max_pre) = @_;
    my $npre = 0;

    open(FILE, "$file") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s*$/) {
	    if ($max_pre < 0 || $max_pre > $npre) {
		my $i = $1;
		my $j = $2;
		if (abs($i-$j) >= $minL) {
		    $pre_i_ref->[$npre] = $i;
		    $pre_j_ref->[$npre] = $j;
		    $npre ++;
		}
	    }
	}
    }
    close(FILE);

    if ($verbose) {
	print "n predictions: $npre\n";
	for (my $p = 0; $p < $npre; $p ++) {
	    print "$pre_i_ref->[$p] $pre_j_ref->[$p]\n";
	}
    }
    
    $$ret_npre = $npre;
}

sub parse_rscape_alen {
    my ($file, $ret_alen) = @_;
    my $alen = 0;
    
    open(FILE, "$file") || die;
    while(<FILE>) {
	if (/^\# MSA .+ alen (\d+)\s/) {
	    $alen = $1;
	}
    }
    close(FILE);
    
    $$ret_alen = $alen;
}

sub parse_mfdca {
    my ($file, $minL, $ret_npre, $pre_i_ref, $pre_j_ref, $max_pre, $which) = @_;
    my $npre = 0;

    my $sortfile = sort_mfDCA($file, $which);
    if (0) { system("more $sortfile\n"); }
    
    open(FILE, "$sortfile") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+\S+\s*$/) {
	    my $i = $1;
	    my $j = $2;
	    if ($max_pre < 0 || $max_pre > $npre) {
		if (abs($i-$j) >= $minL) {
		    $pre_i_ref->[$npre] = $i;
		    $pre_j_ref->[$npre] = $j;
		    $npre ++;
		}
	    }
	}
    }
    close(FILE);

    if ($verbose) {
	print "n predictions: $npre\n";
	for (my $p = 0; $p < $npre; $p ++) {
	    print "$pre_i_ref->[$p] $pre_j_ref->[$p]\n";
	}
    }
    
    $$ret_npre = $npre;
    system("rm $sortfile\n");
}

sub parse_plmdca {
    my ($file, $minL, $ret_npre, $pre_i_ref, $pre_j_ref, $max_pre) = @_;
    my $npre = 0;

    my $sortfile = sort_plmDCA($file);

    open(FILE, "$sortfile") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+\S+\s*$/) {
	    my $i = $1;
	    my $j = $2;
	    if ($max_pre < 0 || $max_pre > $npre) {
		if (abs($i-$j) >= $minL) {
		    $pre_i_ref->[$npre] = $i;
		    $pre_j_ref->[$npre] = $j;
		    $npre ++;
		}
	    }
	}
    }
    close(FILE);

    if ($verbose) {
	print "n predictions: $npre\n";
	for (my $p = 0; $p < $npre; $p ++) {
	    print "$pre_i_ref->[$p] $pre_j_ref->[$p]\n";
	}
    }
    
    $$ret_npre = $npre;
    system("rm $sortfile\n");
}

sub parse_gremplin {
    my ($file, $minL, $ret_npre, $pre_i_ref, $pre_j_ref, $max_pre) = @_;
    my $npre = 0;

    my $sortfile = sort_gremlin($file);

    open(FILE, "$sortfile") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/(\d+)\s+(\d+)\s+\S+\s*$/) {
	    my $i = $1;
	    my $j = $2;
	    if ($max_pre < 0 || $max_pre > $npre) {
		if (abs($i-$j) >= $minL) {
		    $pre_i_ref->[$npre] = $i;
		    $pre_j_ref->[$npre] = $j;
		    $npre ++;
		}
	    }
	}
    }
    close(FILE);

    if ($verbose) {
	print "n predictions: $npre\n";
	for (my $p = 0; $p < $npre; $p ++) {
	    print "$pre_i_ref->[$p] $pre_j_ref->[$p]\n";
	}
    }
    
    $$ret_npre = $npre;
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

sub parse_cmp {
    my ($file, $minL, $ret_ncmp, $cmp_i_ref, $cmp_j_ref) = @_;
    my $ncmp = $$ret_ncmp;

    my $n;
    open(FILE, "$file") || die;
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/^(\d+)\s+\S+\s+(\d+)\s+\S+\s+\S+\s*$/) {
	    my $i = $1;
	    my $j = $2;
	    if (abs($i-$j) >= $minL) {
		for ($n = 0; $n < $ncmp; $n ++) {
		    if ($i == $cmp_i_ref->[$n] && $j == $cmp_j_ref->[$n]) { last; }
		    if ($j == $cmp_i_ref->[$n] && $i == $cmp_j_ref->[$n]) { last; }
		}
		if ($n == $ncmp) {
		    $cmp_i_ref->[$ncmp] = $1;
		    $cmp_j_ref->[$ncmp] = $2;
		    $ncmp ++;
		}
	    }
	}
    }
    close(FILE);
    
    print "n contacts: $ncmp\n";
    if ($verbose) {
	print "n contacts: $ncmp\n";
	for (my $c = 0; $c < $ncmp; $c ++) {
	    print "$cmp_i_ref->[$c] $cmp_j_ref->[$c]\n";
	}
    }
    
    $$ret_ncmp = $ncmp;
}

sub parse_cor {
    my ($file, $minL, $o2map_ref, $ret_lenPDB, $mapPDB_ref) = @_;

    my $lenPDB = 0;
    my $i;
    my $j;
    open(FILE, "$file") || die; 
    while(<FILE>) {
	if (/^\#/) {
	}
	elsif (/^(\d+)\s+(\d+)\s*$/) {
	    $i = $1-1;
	    $j = $o2map_ref->[$2]-1;
	    
	    $mapPDB_ref->[$i] = $j;
	}
    }
    close(FILE);
    $lenPDB = $i + 1;

    if ($verbose) {
	print "PDB seq len=$lenPDB\n";
	for (my $l = 0; $l < $lenPDB; $l ++) {
	    print "$l $mapPDB_ref->[$l]\n";
	}
    }
    
    $$ret_lenPDB = $lenPDB;
}

sub plot_roc {
    my ($pfamname, $minL, $F, $file_ref, $xmax, $seeplots) = @_;

    my $psfile = "$pfamname.N$F.ps";
    
    #if ($psfile =~ /\/([^\/]+)\s*$/) { $psfile = "$1"; }
    my $pdffile = $psfile;
    if ($pdffile =~ /^(\S+).ps$/) { $pdffile = "$1.pdf"; }

    my $xlabel = "number of predictions";
    my $ylabel;
    my $title  = "$pfamname [minL $minL]";
    
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    
    print GP "set terminal postscript color solid 14\n";
    FUNCS::gnuplot_define_styles (*GP);

    print GP "set output '$psfile'\n";    
    #print GP "set nokey\n";
    #print GP "set logscale x\n";
    #print GP "set title \"$title\\n\\n$key\"\n";

    my $cmd;
    my $m;
    my $f;

    my $maxpp = 10;
    $cmd = "";
    $m = 1;
    $ylabel = "PPV";
    $xlabel = "SEN";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    print GP "set title '$title'\n";
    print GP "set ylabel '$ylabel'\n";
    print GP "set xrange [0:100]\n";
    print GP "set yrange [0:100]\n";
    for ($f = 0; $f < $F; $f++) {
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 5:6                      ls 7, " : "'$file_ref->[$f]' using 5:6                     ls $m, ";
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 5:6  title '' with lines ls 7"   : "'$file_ref->[$f]' using 5:6 title '' with lines ls $m, ";
	$m ++; if ($m == 7) { $m = 8; }
    }
    print GP "plot $cmd\n";


    $xlabel = "number of predictions per position";
    print GP "set xrange [0:$maxpp]\n";
    print GP "set yrange [0:100]\n";
    $cmd = "";
    $m = 1;
    $ylabel = "PPV";
    print GP "set title '$title'\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    for ($f = 0; $f < $F; $f++) {
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 8:6                      ls 7, " : "'$file_ref->[$f]' using 8:6                     ls $m, ";
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 8:6  title '' with lines ls 7"   : "'$file_ref->[$f]' using 8:6 title '' with lines ls $m, ";
	$m ++; if ($m == 7) { $m = 8; }
    }
    print GP "plot $cmd\n";
    $xlabel = "number of predictions";
    print GP "set xrange [1:$xmax]\n";
    print GP "set yrange [0:100]\n";
    $cmd = "";
    $m = 1;
    $ylabel = "PPV";
    print GP "set title '$title'\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    for ($f = 0; $f < $F; $f++) {
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:6                      ls 7, " : "'$file_ref->[$f]' using 1:6                     ls $m, ";
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:6  title '' with lines ls 7"   : "'$file_ref->[$f]' using 1:6 title '' with lines ls $m, ";
	$m ++; if ($m == 7) { $m = 8; }
    }
    print GP "plot $cmd\n";
    $cmd = "";
    $m = 1;
    $ylabel = "SEN";
    print GP "set title '$title'\n";
    print GP "set ylabel '$ylabel'\n";
    for ($f = 0; $f < $F; $f++) {
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:5                      ls 7, " : "'$file_ref->[$f]' using 1:5                      ls $m, ";
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:5  title '' with lines ls 7"   : "'$file_ref->[$f]' using 1:5  title '' with lines ls $m, ";
	$m ++; if ($m == 7) { $m = 8; }
    }
    print GP "plot $cmd\n";
    $cmd = "";
    $m = 1;
    $ylabel = "F";
    print GP "unset yrange\n";
    print GP "set title '$title'\n";
    print GP "set ylabel '$ylabel'\n";
    for ($f = 0; $f < $F; $f++) {
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:7                      ls 7, " : "'$file_ref->[$f]' using 1:7                      ls $m, ";
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:7  title '' with lines ls 7"   : "'$file_ref->[$f]' using 1:7  title '' with lines ls $m, ";
	$m ++; if ($m == 7) { $m = 8; }
    }
    print GP "plot $cmd\n";
    
    $xlabel = "number of predictions per position";
    $cmd = "";
    $m = 1;
    $ylabel = "fraction of positions with >= $minC contacts";
    print GP "set title '$title'\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    print GP "set yrange [0:*]\n";
    print GP "set xrange [0:$maxpp]\n";
    for ($f = 0; $f < $F; $f++) {
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 8:9                      ls 7, " : "'$file_ref->[$f]' using 8:9                     ls $m, ";
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 8:9  title '' with lines ls 7"   : "'$file_ref->[$f]' using 8:9 title '' with lines ls $m, ";
	$m ++; if ($m == 7) { $m = 8; }
    }
    print GP "plot $cmd\n";
    $xlabel = "number of predictions";
    print GP "set xrange [1:$xmax]\n";
    print GP "set yrange [0:*]\n";
    $cmd = "";
    $m = 1;
    $ylabel = "fraction of positions with >= $minC contacts";
    print GP "set title '$title'\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$ylabel'\n";
    for ($f = 0; $f < $F; $f++) {
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:9                      ls 7, " : "'$file_ref->[$f]' using 1:9                     ls $m, ";
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:9  title '' with lines ls 7"   : "'$file_ref->[$f]' using 1:9 title '' with lines ls $m, ";
	$m ++; if ($m == 7) { $m = 8; }
    }
    print GP "plot $cmd\n";

    print GP "set logscale x\n";
    $cmd = "";
    $m = 1;
    $ylabel = "PPV";
    print GP "set title '$title'\n";
    print GP "set ylabel '$ylabel'\n";
    for ($f = 0; $f < $F; $f++) {
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:6                     ls 7, " : "'$file_ref->[$f]' using 1:6                     ls $m, ";
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:6 title '' with lines ls 7"   : "'$file_ref->[$f]' using 1:6 title '' with lines ls $m, ";
	$m ++; if ($m == 7) { $m = 8; }
    }
    print GP "plot $cmd\n";
    $cmd = "";
    $m = 1;
    $ylabel = "SEN";
    print GP "set xrange [1:$xmax]\n";
    print GP "set yrange [0:100]\n";
    print GP "set title '$title'\n";
    print GP "set ylabel '$ylabel'\n";
    for ($f = 0; $f < $F; $f++) {
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:5                      ls 7, " : "'$file_ref->[$f]' using 1:5                      ls $m, ";
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:5  title '' with lines ls 7"   : "'$file_ref->[$f]' using 1:5  title '' with lines ls $m, ";
	$m ++; if ($m == 7) { $m = 8; }
    }
    print GP "plot $cmd\n";
    $cmd = "";
    $m = 1;
    $ylabel = "F";
    print GP "set title '$title'\n";
    print GP "set ylabel '$ylabel'\n";
    print GP "unset yrange\n";
    for ($f = 0; $f < $F; $f++) {
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:7                      ls 7, " : "'$file_ref->[$f]' using 1:7                      ls $m, ";
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:7  title '' with lines ls 7"   : "'$file_ref->[$f]' using 1:7  title '' with lines ls $m, ";
	$m ++; if ($m == 7) { $m = 8; }
    }
    print GP "plot $cmd\n";

    print GP "set logscale x\n";
    $cmd = "";
    $m = 1;
    $ylabel = "fraction of positions with >= $minC contacts";
    print GP "set title '$title'\n";
    print GP "set ylabel '$ylabel'\n";
    print GP "set yrange [0:*]\n";
 
    for ($f = 0; $f < $F; $f++) {
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:9                     ls 7, " : "'$file_ref->[$f]' using 1:9                     ls $m, ";
	$cmd .= ($f == $F-1)? "'$file_ref->[$f]' using 1:9 title '' with lines ls 7"   : "'$file_ref->[$f]' using 1:9 title '' with lines ls $m, ";
	$m ++; if ($m == 7) { $m = 8; }
    }
    print GP "plot $cmd\n";

    close (GP);

    system ("/bin/ps2pdf $psfile $pdffile\n"); 
    system("rm $psfile\n");
    if ($seeplots) { system ("open $pdffile&\n"); }
}

sub parse_afa {
    my ($afafile, $map_ref, $ret_alen, $ret_alen_dca) = @_;

    my $alen     = 0;
    my $asq      = "";
    my $alen_dca = 0;
    my $asq_dca  = "";

    # grab the first sequence
    my $n = 0;
    print "afafile:$afafile\n";
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
	    $asq_dca .= "$char"; 
	    $pos ++;     
	    $map_ref->[$pos] = $apos;
	}
    }
    $alen_dca = length($asq_dca);

    print "alen     $alen\n$asq\n";
    print "alen_dca $alen_dca\n$asq_dca\n";
   
    $$ret_alen_dca = $alen_dca;
    $$ret_alen     = $alen;
}

sub parse_filtermap {
    my ($filtermapfile, $map2o_ref, $o2map_ref, $alen, $oalen) = @_;

    if ($filtermapfile) {
	my $newalen = 0;
	for (my $i = 0; $i <= $alen;  $i ++) { $map2o_ref->[$i] = 0; }
	for (my $i = 0; $i <= $oalen; $i ++) { $o2map_ref->[$i] = 0; }
	
	my $n = 0;
	open (MAP, "$filtermapfile") || die;
	while (<MAP>) {
	    if ($n == 0 && /^\#.+alen=(\d+)/) {
		$newalen = $1;
		if  ($n==0 && $alen != $newalen) { print "bad map file. $newalen should be $alen\n"; die; }
	    }
	    elsif (/^(\d+)\s+(\d+)\s*$/) {
		my $i    = $1+1;
		my $mapi = $2+1;
		$map2o_ref->[$i]    = $mapi;
		$o2map_ref->[$mapi] = $i;
	    }  
	    $n ++;
	}
	close (MAP);
    }
    else {
	for (my $i = 0; $i <= $alen;  $i ++) { $map2o_ref->[$i] = $i; }
	for (my $i = 0; $i <= $oalen; $i ++) { $o2map_ref->[$i] = $i; }
    }
    if (0) {
	for (my $i = 0; $i <= $alen;  $i ++) { print "^^$i $map2o_ref->[$i]\n"; }
	for (my $i = 0; $i <= $oalen; $i ++) { print "++$i $o2map_ref->[$i]\n"; }
     }
}

