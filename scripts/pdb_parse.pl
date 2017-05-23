#!/usr/bin/perl -w
#pdb_parse.pl

use strict;
use Class::Struct;

use vars qw ($opt_C $opt_M $opt_W $opt_D $opt_L $opt_v);  # required if strict used
use Getopt::Std;
getopts ('C:M:W:D:L:v');

struct RES => {
char     => '$', # number of atoms
coor     => '$', # coordenate in seqref
nat      => '$', # number of atoms
type     => '@', # type for eacth atom
x        => '@', # x position for eacht atom
y        => '@', # y position for eacht atom
z        => '@', # z position for eacht atom
};

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  pdb_parse.pl [options] <pdbfile> <stofile> <rscapedir> <gnuplotdir> \n\n";
        print "options:\n";
 	exit;
}

my $pdbfile   = shift;
my $stofile   = shift;
my $rscapebin = shift;
my $gnuplot   = shift;
use constant GNUPLOT => '$gnuplot';

my $currdir = $ENV{PWD};

my $stoname = $stofile;
if ($stoname =~ /\/([^\/]+)\.[^\/]+$/) { $stoname = $1; }
my $pfamname = $stofile;
if ($pfamname =~ /\/([^\/]+)\.[^\/]+$/) { $pfamname = $1; }

print "STO:  $stoname\n";
print "PFAM: $pfamname\n";
my $hmmer       = "$rscapebin/../lib/hmmer";
my $hmmbuild    = "$hmmer/src/hmmbuild";
my $hmmersearch = "$hmmer/src/hmmsearch";

my $easel    = "$hmmer/easel";
my $sfetch   = "$easel/miniapps/esl-sfetch";
my $reformat = "$easel/miniapps/esl-reformat";

my $coorfile = "";
if ($opt_C) { 
    $coorfile = "$opt_C";
    open(COORF,  ">$coorfile")  || die;
}
my $mapfile = "";
if ($opt_M) { 
    $mapfile = "$opt_M";
    open(MAPF,  ">$mapfile")  || die;
}

my $maxD = 8;
if ($opt_D) { $maxD = $opt_D; }
$maxD = int($maxD*100)/100;

my $minL = 1;
if ($opt_L) { $minL = $opt_L; }

#options: CA C ALL NOH
my $which = "CA";
if ($opt_W) { $which = "$opt_W"; }

my $nch = 0;
my @chname = ();
my @sq = ();
my $len;
my $alen;

my $seeplots = 0;

my $pdbname = parse_pdb($pdbfile, \$nch, \@chname, \@sq);
print     "# PFAM: $stofile\n";
print     "# PDB:  $pdbname\n# $nch chains\n";

for (my $n = 0; $n < $nch; $n ++) {
    my $map0file = "$pdbfile.chain$chname[$n].maxD$maxD.map";
    my $map1file = "$pdbfile.chain$chname[$n].maxD$maxD.$pfamname.map";
    my $mapfile  = "$currdir/$stoname.chain$chname[$n].maxD$maxD.map";
    my $corfile  = "$currdir/$stoname.chain$chname[$n].maxD$maxD.cor";

    print "$mapfile\n";
    if ($coorfile) {
	print COORF "$corfile\n";
    }
    
    open(COR,  ">$corfile")  || die;
    open(MAP,  ">$mapfile")  || die;
    open(MAP0, ">$map0file") || die;
    open(MAP1, ">$map1file") || die;
    
    #print      "# PDB: $pdbname\n";
    print COR  "# PDB: $pdbname\n";
    print MAP  "# PDB: $pdbname\n"; if ($mapfile)  { print MAPF  "# PDB: $pdbname\n"; }
    print MAP0 "# PDB: $pdbname\n";
    print MAP1 "# PDB: $pdbname\n";
    #print      "# chain $chname[$n]\n# $sq[$n]\n";
    print COR  "# chain $chname[$n]\n";
    print MAP  "# chain $chname[$n]\n"; if ($mapfile)  { print MAPF  "# chain $chname[$n]\n"; }
    print MAP0 "# chain $chname[$n]\n";
    print MAP1 "# chain $chname[$n]\n";
    print MAP1 "# maps to $pfamname\n";
    
    $len = length($sq[$n]);
    $alen = parse_pdb_contact_map($pdbfile, $stofile, $chname[$n], $sq[$n], $which, $maxD, $minL);
    close(COR);
    close(MAP);
    close(MAP0);
    close(MAP1);
    if ($alen == 0) { next; }
    
    plot_contact_map($map0file, $len,  $seeplots);
    plot_contact_map($map1file, -1,    $seeplots);
    plot_contact_map($mapfile,  $alen, $seeplots);
}

if ($coorfile) { close(COORF); }
if ($mapfile)  { close(MAPF);  }

sub parse_pdb {
    my ($pdbfile, $ret_nch, $chname_ref, $sq_ref) = @_;

    my $pdfname;
    my $nch = 0;

    my $cur_chain;
    my $prv_chain = "";
    my $sq;
    my @sqlen;

    open(FILE, "$pdbfile") || die;
    while (<FILE>) {
	if (/^HEADER\s+.+\s+(\S+)\s*$/) {
	    $pdbname = $1;
	}
	elsif (/^SEQRES\s+\d+\s+(\S+)\s+(\d+)\s+(\S+.+)$/) {
	    $cur_chain   = $1;
	    my $sqlen    = $2;
	    $sq          = $3;
	    $sq =~ s/\n//g;
	    
	    if ($prv_chain =~ /^$/) {
		$sq_ref->[$nch]     = $sq;
		$chname_ref->[$nch] = $cur_chain;
		$sqlen[$nch]        = $sqlen;
	    }
	    elsif ($cur_chain =~ /^$prv_chain$/) {
		$sq_ref->[$nch] .= $sq;
	    }
	    else { 
		$nch ++; 
		$sq_ref->[$nch]     = $sq;
		$chname_ref->[$nch] = $cur_chain;
		$sqlen[$nch]        = $sqlen;
	    }
	    $prv_chain = $cur_chain;
	}
    }
    close(FILE);
    $nch ++;

    for (my $n = 0; $n < $nch; $n ++) {
	my $len = sq_conversion(\$sq_ref->[$n]);
	if ($len != $sqlen[$n]) { print "parse_pdb(): seq len is $len should be $sqlen[$n]\n"; die; }
    }
    
    $$ret_nch = $nch;
    return $pdbname;
}

sub map_pdbsq {
    my ($stofile, $pdbname, $chname, $pdbsq, $map_ref) = @_;

    my $len = length($pdbsq);
    for (my $l = 0; $l < $len; $l ++) {
	$map_ref->[$l] = -1;
    }

    my $pfam_name = "";
    my $pfam_asq = "";
    my $i;
    my $j;
    open(STO, "$stofile") || die;
    while (<STO>) {
	if (/^\#\=GS\s+(\S+)\s+DR PDB\;\s+\S*$pdbname\S*\s+$chname\;\s+(\d+)\-(\d+)\;/) {
	    $pfam_name = $1;
	    $i         = $2;
	    $j         = $3;
	    #print "# $pfam_name $i $j\n";
	}
	elsif (/^\#\=GS\s+(\S+)\s+AC\s+\S*$pdbname\S*\s*$/) {
	    $pfam_name = $1;
	    $i         = 1;
	    $j         = -1;
	    #print "# $pfam_name\n";
	}
	elsif (/^$pfam_name\s+(\S+)\s*$/) {
	    $pfam_asq .= uc $1;
	}
    }
    close(STO);

    my $from_pdb;
    my $from_pfam;
    my $ali_pdb;
    my $ali_pfam;
    my $alen = find_pdbsq_in_pfam($stofile, $chname, $pdbname, $pdbsq, \$pfam_name, \$from_pdb, \$ali_pdb, \$from_pfam, \$ali_pfam, \$pfam_asq);
    if ($alen == 0 || $pfamname =~ //) {
	print "could not find $pdbname chain $chname in sto file\n";
	return 0;
    }

    my @ali_pdb  = split(//,$ali_pdb);
    my @ali_pfam = split(//,$ali_pfam);
    my @pfam_asq = split(//,$pfam_asq);
    my $x = $from_pdb-1;
    my $n = 0;
    
    my $y = 0;
    while ($n < $from_pfam-1) {
	if ($pfam_asq[$y] =~ /^[\.\-]$/) { 
	    printf "# skip pfam gap $y %s \n", $pfam_asq[$y]; 
	    $y ++;
	}
	else {
	    $n ++; $y ++;
 	}
    }
    print "^^ x $x y $y\n";
    my $pos = 0;
    while ($pos < $alen) {
	my $pos_pdb  = uc($ali_pdb[$pos]);
	my $pos_pfam = uc($ali_pfam[$pos]);
	#printf("\n^^^^ pos $pos | $pos_pdb $pos_pfam | pfam $y %s\n", $pfam_asq[$y]);
	if ($pos_pdb =~ /^[\.\-]$/  && $pos_pfam =~ /^[\.\-]$/  ) { 
	    printf "# double gap at pos %d\n", $pos; 
	}
	elsif ($pos_pdb =~ /^[\.\-]$/)  { 
	    printf "# pdb gap | move pfam %d $pos_pfam\n", $y;
	    while($pfam_asq[$y] =~ /^[\.\-]$/) { 
		printf "# skip pfam gap $y %s \n", $pfam_asq[$y]; 
		$y ++; 
	    }
	    $y ++;	
	}
	elsif ($pos_pfam =~ /^[\.\-]$/)  { 
	    printf "# pfam gap | move pdb $x $pos_pdb \n"; 
	    $x ++; 
	}
	elsif ($pos_pfam =~ /^$pos_pdb$/)  { 
	    while($pfam_asq[$y] =~ /^[\.\-]$/) { 
		printf "# skip pfam gap $y %s \n", $pfam_asq[$y]; 
		$y ++; 
	    }
	    $map_ref->[$x] = $y;  
	    printf "# match $x $pos_pdb | %d $pos_pfam\n", $y; 
	    $x ++; $y ++; 
	}
	else {
	    while($pfam_asq[$y] =~ /^[\.\-]$/) { 
		printf "# skip pfam gap $y %s \n", $pfam_asq[$y]; 
		$y ++; 
	    }
	    $map_ref->[$x] = $y;  
	    printf "# mismach $x $pos_pdb | %d $pos_pfam\n", $y; 
	    $x ++; $y ++;
	}
	$pos ++;
    }
    
    return length($pfam_asq);
}

sub alipos_isgap {
    my ($char) = @_;

    if ($char =~ /^[\.\-]$/) {
	return 1;
    }
    return 0;
}

sub find_pdbsq_in_pfam {
    my ($stofile, $chain, $pdbname, $pdbsq, $ret_pfamname, $ret_from_pdb, $ret_ali_pdb, $ret_from_pfam, $ret_ali_pfam, $ret_pfam_asq) = @_;
 
    my $ali_pdb  = "";
    my $ali_pfam = "";
    my $from_pdb;
    my $from_pfam;
    my $pfam_asq = "";
    my $eval = 1000;

    my $pfamname =  $$ret_pfamname;
     if ($pfamname eq "") {    
	print "could not find $pdbname chain $chain in sto file... trying by homology\n";
    }

    my $pdbsqfile = "$currdir/$pdbname";
    
    open(F, ">$pdbsqfile") || die;
    print F ">$pdbname\n$pdbsq\n";
    close(F);
     
    my $hmm    = "$currdir/hmmfile";
    my $hmmout = "$currdir/hmmout";
 
    system("$hmmbuild             $hmm $pdbsqfile >  $hmmout\n");
    system("$hmmersearch -E $eval $hmm $stofile   >  $hmmout\n");
    system("/bin/echo $hmmersearch -E $eval --F1 10 --F2 10 --F3 10 $hmm $stofile   \n");
    
    #system("/usr/bin/more $hmmout\n");

    # take hit with best evalue
    my $pdbasq;
    my $asq;
    parse_hmmout_for_besthit($hmmout, $pdbname, \$pfamname, \$from_pdb, \$ali_pdb, \$from_pfam, \$ali_pfam);
    if ($pfamname eq "") { print "could not find best hit for chain $chain\n"; }
    
    $pfam_asq = get_asq_from_sto($stofile, $pfamname, 0);

    if ($pfamname ne "") {
	print "^^>$pdbname\n$pdbsq\n";
	print ">^^$pfamname\n$pfam_asq\n";
	print "^^$ali_pdb\n";
	print "^^$ali_pfam\n";
    }
    
    $$ret_pfamname  = $pfamname;
    $$ret_ali_pdb   = $ali_pdb;
    $$ret_ali_pfam  = $ali_pfam;
    $$ret_from_pdb  = $from_pdb;
    $$ret_from_pfam = $from_pfam;
    $$ret_pfam_asq  = $pfam_asq;

    system("/bin/rm $pdbsqfile\n");
    #system("/bin/rm $hmm\n");
    system("/bin/rm $hmmout\n");

    return length($ali_pfam);
}

    
sub parse_hmmout_for_besthit {
    my ($hmmout, $pdbname, $ret_pfamname, $ret_from_pdb, $ret_ali_pdb, $ret_from_pfam, $ret_ali_pfam) = @_;

    my $from_pdb  = -1;
    my $to_pdb    = -1;
    my $from_pfam = -1;
    my $to_pfam   = -1;
    my $ali_pdb   = "";
    my $ali_pfam  = "";
    my $pfamname  = "";

    my $n = 0;
    open(HF, "$hmmout") || die;
    while(<HF>) {
	if (/#/) {
	}
	elsif (/\>\>\s+(\S+)/ && $n == 0) {
	    $pfamname = $1;
	    $ali_pfam = "";
	    $ali_pdb  = "";
	    $from_pfam = 123456789;
	    $from_pdb  = 123456789;
	    $n ++;
	}
	elsif (/\>\>\s+(\S+)/ && $n > 0) {
	    last;
	}
	elsif ($n==1 && /^\s*$pfamname\s+(\d+)\s+(\S+)\s+(\d+)\s*$/) {
	    my $i  = $1;
	    $ali_pfam  .= $2;
	    $to_pfam    = $3;
	    $from_pfam  = ($i < $from_pfam)? $i:$from_pfam;
	}
	elsif ($n==1 && /^\s*$pdbname\s+(\d+)\s+(\S+)\s+(\d+)\s*$/) {
	    my $i     = $1;
	    $ali_pdb .= $2;
	    $to_pdb   = $3;
	    $from_pdb = ($i < $from_pdb)? $i:$from_pdb;
	}
    }
    close(HF);
    
    $$ret_from_pdb  = $from_pdb;
    $$ret_from_pfam = $from_pfam;
    $$ret_pfamname  = ($from_pfam<0)? "":$pfamname;
    $$ret_ali_pdb   = $ali_pdb;
    $$ret_ali_pfam  = $ali_pfam;
}

sub coordtrans_seq2asq{
    my ($i, $sq, $asq, $ret_ai) = @_;

    if ($i > length($sq)) { printf "sq2asq bad coordinate i=$i len=%d\n", length($sq); die; }

    my @sq  = split(//, $sq);
    my @asq = split(//, $asq);
    
    my $char;
    my $ai = 0;
    my $x  = 0;
    while ($asq) {
	$asq =~ s/^(\S+)//; $char = 1;
	if    ($char =~ /\./ || $char =~ /\-/ || $char =~ /\_/) { $ai ++; }
	elsif ($char =~ /^(\S+)$/)                              { $ai ++; $x ++; }
	if ($x == $i) { last; }
    }
    if ($sq[$i] =~ /^$asq[$ai]$/) { } else { printf "sq2asq seqchar %s should be asqchar %s", $sq[$i], $asq[$ai]; die; }
    
    $$ret_ai = $ai;
}

sub coordtrans_aseq2sq{
    my ($ai, $asq, $sq, $ret_i) = @_;

    if ($ai > length($asq)) { printf "asq2sq bad coordinate ai=$ai alen=%d\n", length($asq); die; }
    
    my @sq  = split(//, $sq);
    my @asq = split(//, $asq);

    my $char;
    my $i = 0;
    my $z  = 0;
    while ($asq) {
	$asq =~ s/^(\S+)//; $char = 1;
	if    ($char =~ /\./ || $char =~ /\-/ || $char =~ /\_/) { $z ++; }
	elsif ($char =~ /^(\S+)$/)                              { $z ++; $i ++; }
	if ($z == $ai) { last; }
    }
    if ($sq[$i] =~ /^$asq[$ai]$/) { } else { printf "asq2sq seqchar %s should be asqchar %s", $sq[$i], $asq[$ai]; die; }

    $$ret_i = $i;
}

sub get_asq_from_sto {
    my ($stofile, $name, $from) = @_;
    my $asq = "";

    if ($name =~ //) { return $asq; }
     
    my $afafile = "$stofile";
    if    ($afafile =~ /^(\S+).txt$/) { $afafile = "$1.afa"; }
    elsif ($afafile =~ /^(\S+).sto$/) { $afafile = "$1.afa"; }
    if (-f $afafile) { } else { system("$reformat afa $stofile > $afafile\n"); }

    my $found = 0;
    open(AFA, "$afafile") || die;
    while(<AFA>) {
	if (/^\>$name/) {
	    $found = 1;
	}
	elsif ($found == 1 && /^\>/) {
	    $found = 0;
	    last;
	}       
	elsif ($found == 1 && /^(\S+)$/) {
	    $asq .= $1;
	}
     }
    close(AFA);
    

    return substr($asq, $from);
}


sub parse_pdb_contact_map {
    my ($pdbfile, $stofile, $chain, $sq, $which, $maxD, $minL) = @_;

    my $len = length($sq);
    my @sq = split(//,$sq);
    
    #printf      "# maxD  $maxD\n";
    #printf      "# minL  $minL\n";
    printf COR  "# maxD  $maxD\n";
    printf COR  "# minL  $minL\n";
    printf MAP  "# maxD  $maxD\n";
    printf MAP  "# minL  $minL\n";
    if ($mapfile)  { 
	print MAPF  "# maxD  $maxD\n";
	print MAPF  "# minL  $minL\n";
    }
    printf MAP0 "# maxD  $maxD\n";
    printf MAP0 "# minL  $minL\n";
    printf MAP1 "# maxD  $maxD\n";
    printf MAP1 "# minL  $minL\n";

    my @map = (); # map sq to the sequence in the alignment  
    my $alen = map_pdbsq($stofile, $pdbname, $chain, $sq, \@map);
    if ($alen == 0) { return $alen; }
    for (my $x = 0; $x < $len; $x ++) {
	printf COR "%d %d\n", $x+1, $map[$x]+1;
    }
 
    my $nct = 0;
    my $nc = 0;
    my $nat1;
    my $nat2;
    my $coor1;
    my $coor2;
    my $char1;
    my $char2;
    my @type1;
    my @type2;
    my @x1;
    my @x2;
    my @y1;
    my @y2;
    my @z1;
    my @z2;
    my $l1;
    my $l2;
    my $distance;
    my $atom_offset = atom_offset($pdbfile, $chain);
    if ($atom_offset < 0) { print "could not find atom offset\n"; die; }
    print "#atom offset $atom_offset chain $chain\n";

    my @res;
    get_atoms_coord($pdbfile, $atom_offset, \@sq, $len, $chain, \@res);
 
    for ($l1 = 0; $l1 < $len; $l1 ++) {
	$nat1  = $res[$l1]->{"RES::nat"};
	$coor1 = $res[$l1]->{"RES::coor"};
	$char1 = $res[$l1]->{"RES::char"};
	@type1 = @{$res[$l1]->{"RES::type"}};
	@x1    = @{$res[$l1]->{"RES::x"}};
	@y1    = @{$res[$l1]->{"RES::y"}};
	@z1    = @{$res[$l1]->{"RES::z"}};

	if (aa_conversion($char1) ne $sq[$l1]) {
	    printf "#chain$chain;position %d/%d %d/%d: mismatched character %s should be $sq[$coor1]\n", 
	    ($l1+1), $len, $coor1+1, $len+$atom_offset, aa_conversion($char1); die; }
       
	for ($l2 = $l1+1; $l2 < $len; $l2 ++) {
	    $nat2  = $res[$l2]->{"RES::nat"};
	    $coor2 = $res[$l2]->{"RES::coor"};
	    $char2 = $res[$l2]->{"RES::char"};
	    @type2 = @{$res[$l2]->{"RES::type"}};
	    @x2    = @{$res[$l2]->{"RES::x"}};
	    @y2    = @{$res[$l2]->{"RES::y"}};
	    @z2    = @{$res[$l2]->{"RES::z"}};

	    if (aa_conversion($char2) ne $sq[$l2]) {
		printf "#chain$chain;position %d/%d %d/%d: mismatched character %s should be $sq[$coor2]\n", 
		($l2+1), $len, $coor2+1, $len+$atom_offset, aa_conversion($char2); die; }

	    $distance = distance($which, $nat1, \@type1, \@x1, \@y1, \@z1, $nat2, \@type2, \@x2, \@y2, \@z2);

	    if ($distance > 0) {
		if (abs($l1-$l2) >= $minL && $distance <= $maxD) {
		    printf MAP0 "%d %s %d %s %.2f\n", $l1+1, $sq[$l1], $l2+1, $sq[$l2], $distance;
		    $nct ++;
		    
		    if ($map[$l1] >= 0 && $map[$l2] >= 0) {
			$nc ++;
			printf      "%d %d %s | %d %d %s | %.2f\n", $l1+1, $map[$l1]+1, $sq[$l1], $l2+1, $map[$l2]+1, $sq[$l2], $distance;
			if ($mapfile) {
			    printf MAPF  "%d %s %d %s %.2f\n", $map[$l1]+1, $sq[$l1], $map[$l2]+1, $sq[$l2], $distance; 
			}
			printf MAP  "%d %s %d %s %.2f\n", $map[$l1]+1, $sq[$l1], $map[$l2]+1, $sq[$l2], $distance; 
			printf MAP1 "%d %s %d %s %.2f\n", $l1+1, $sq[$l1], $l2+1, $sq[$l2], $distance; 
		    }
		}
	    }
	}
    }
    #print      "\# contacts: $nc\n";
    print MAP  "\# contacts: $nc\n";
    if ($mapfile) {
	print MAP  "\# contacts: $nc\n";
    }
    print MAP0 "\# contacts: $nct\n";
    print MAP1 "\# contacts: $nc\n";

    return $alen;
}

sub atom_offset {
    my ($pdbfile, $chain) = @_;

    my $atom_offset = 1;
    open(FILE, "$pdbfile") || die;
    while (<FILE>) {
	if (/^DBREF1\s+\S+\s+$chain\s+(\S+)\s+\S+\s+\S+\s*/) {
	    $atom_offset = $1;
	    last;
	}
	elsif (/^DBREF\s+\S+\s+$chain\s+(\S+)\s+\S+\s+\S+\s*/) {
	    $atom_offset = $1;
	    last;
	}
    }
    close(FILE);
    
    open(FILE, "$pdbfile") || die;
    while (<FILE>) {
	if (/SEQADV\s+\S+\s+\S+\s+$chain\s+(\S+)\s+/) {
	    my $val = $1;
	    if ($val < $atom_offset) { $atom_offset = $val; }
	}
    }
    close(FILE);

    return $atom_offset;
}

sub sq_conversion {
    my ($ret_seq) = @_;

    my $seq = $$ret_seq;
    my $new = "";
    my $aa;

    #print "seq\n$seq\n";
    while ($seq) {
	$seq =~ s/^(\S+)\s+//; $aa = $1;
	$new .= aa_conversion($aa); 
    }
    #printf "new $new\n";
    $$ret_seq = $new;
    return length($new);
}

sub aa_conversion {
    my ($aa) = @_;
    my $new;

    if ($aa =~ /^\S$/) { return $aa; } # DNA/RNA
    
    if    ($aa =~ /^ALA$/)  { $new = "A"; }
    elsif ($aa =~ /^CYS$/)  { $new = "C"; }
    elsif ($aa =~ /^ASP$/)  { $new = "D"; }
    elsif ($aa =~ /^GLU$/)  { $new = "E"; }
    elsif ($aa =~ /^BGLU$/) { $new = "E"; }
    elsif ($aa =~ /^PHE$/)  { $new = "F"; }
    elsif ($aa =~ /^GLY$/)  { $new = "G"; }
    elsif ($aa =~ /^HIS$/)  { $new = "H"; }
    elsif ($aa =~ /^ILE$/)  { $new = "I"; }
    elsif ($aa =~ /^LYS$/)  { $new = "K"; }
    elsif ($aa =~ /^LEU$/)  { $new = "L"; }
    elsif ($aa =~ /^MET$/)  { $new = "M"; }
    elsif ($aa =~ /^BMET$/) { $new = "M"; }
    elsif ($aa =~ /^MSE$/)  { $new = "M"; } #selenomethionine is incorporated as methionine
    elsif ($aa =~ /^ASN$/)  { $new = "N"; }
    elsif ($aa =~ /^PRO$/)  { $new = "P"; }
    elsif ($aa =~ /^GLN$/)  { $new = "Q"; }
    elsif ($aa =~ /^ARG$/)  { $new = "R"; }
    elsif ($aa =~ /^SER$/)  { $new = "S"; }
    elsif ($aa =~ /^SEP$/)  { $new = "S"; } # phosphoseronine
    elsif ($aa =~ /^THR$/)  { $new = "T"; }
    elsif ($aa =~ /^TPO$/)  { $new = "T"; } # phosphothereonine
    elsif ($aa =~ /^VAL$/)  { $new = "V"; }
    elsif ($aa =~ /^TRP$/)  { $new = "W"; }
    elsif ($aa =~ /^TYR$/)  { $new = "Y"; }
    elsif ($aa =~ /^DA$/)   { $new = "A"; }
    elsif ($aa =~ /^DC$/)   { $new = "C"; }
    elsif ($aa =~ /^DG$/)   { $new = "G"; }
    elsif ($aa =~ /^DT$/)   { $new = "T"; }
    elsif ($aa =~ /^A$/)    { $new = "A"; }
    elsif ($aa =~ /^C$/)    { $new = "C"; }
    elsif ($aa =~ /^G$/)    { $new = "G"; }
    elsif ($aa =~ /^T$/)    { $new = "T"; }
    elsif ($aa =~ /^U$/)    { $new = "U"; }
    # modified residues
    elsif ($aa =~ /^1MA$/)  { $new = "A"; }
    elsif ($aa =~ /^12A$/)  { $new = "A"; }
    elsif ($aa =~ /^5MC$/)  { $new = "C"; }
    elsif ($aa =~ /^CCC$/)  { $new = "C"; }
    elsif ($aa =~ /^GDP$/)  { $new = "G"; }
    elsif ($aa =~ /^GTP$/)  { $new = "G"; }
    elsif ($aa =~ /^2MG$/)  { $new = "G"; }
    elsif ($aa =~ /^7MG$/)  { $new = "G"; }
    elsif ($aa =~ /^H2U$/)  { $new = "U"; }
    elsif ($aa =~ /^UMP$/)  { $new = "U"; }
    elsif ($aa =~ /^PSU$/)  { $new = "U"; }
    elsif ($aa =~ /^2MU$/)  { $new = "U"; }
    elsif ($aa =~ /^70U$/)  { $new = "U"; }
    elsif ($aa =~ /^AH2U$/) { $new = "U"; }
    elsif ($aa =~ /^BH2U$/) { $new = "U"; }
    else { print "aa_conversion(): uh? $aa\n"; die; }

    return $new;
}


sub get_atoms_coord {
    my ($pdbfile, $atom_offset, $pdbsq_ref, $len, $chain, $res_ref) = @_;

    my $nat = 0;
    my $type;
    my $char = "9";
    my $coor;
    my $x;
    my $y;
    my $z;
    my $l;

    for ($l = 0; $l < $len; $l ++) {
	$res_ref->[$l] = RES->new();
	$res_ref->[$l]->{"RES::nat"}  = 0;
	$res_ref->[$l]->{"RES::coor"} = $l+$atom_offset-1;
	$res_ref->[$l]->{"RES::char"} = "$pdbsq_ref->[$l]";
    }

    open(FILE, "$pdbfile") || die;
    while (<FILE>) {
	if (  /^HETATM\s+\d+\s+(\S+)\s+HOH\s+$chain\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/  ||
	    /^HETATM\s+\d+\s+MN\s+$chain\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/   ) {
	}
	elsif (  /^ATOM\s+\d+\s+(\S+)\s+(\S+)\s+$chain\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/  || 
		 /^ATOM\s+\d+\s+(\S{3})(\S{4})\s+$chain\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/ ||
	    /^HETATM\s+\d+\s+(\S+)\s+(\S+)\s+$chain\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/  ||
	    /^HETATM\s+\d+\s+(\S{3})(\S{4})\s+$chain\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/   ) {

	    $type  = $1;
	    $char  = $2;
	    $l     = $3-$atom_offset;
	    $x     = $4;
	    $y     = $5;
	    $z     = $6;
	    
	    if ($l >= 0 && $l < $len) {
		$nat = $res_ref->[$l]->{"RES::nat"};	    
		${$res_ref->[$l]->{"RES::type"}}[$nat] = $type;
		${$res_ref->[$l]->{"RES::x"}}[$nat]    = $x;
		${$res_ref->[$l]->{"RES::y"}}[$nat]    = $y;
		${$res_ref->[$l]->{"RES::z"}}[$nat]    = $z;
		$res_ref->[$l]->{"RES::char"}          = $char;
		$res_ref->[$l]->{"RES::nat"} ++;
	    }
	}
    }
    close(FILE);
    
    for ($l = 0; $l < $len; $l ++) {
	$nat  = $res_ref->[$l]->{"RES::nat"};	    
	$coor = $res_ref->[$l]->{"RES::coor"};	    
	$char = $res_ref->[$l]->{"RES::char"};	    
	if ($nat == 0) { print "#res $l has not atoms\n"; }
	if (1) { printf "res %d coor %d char %s nat %d\n", $l, $coor, $char, $nat; }
    }
}


sub distance {
    my ($which, $nat1, $type1_ref, $x1_ref, $y1_ref, $z1_ref, $nat2, $type2_ref, $x2_ref, $y2_ref, $z2_ref) = @_;

    my $distance = 0;
    if ($nat1 == 0 || $nat2 == 0) { return $distance; }
    
    if  ($which =~ /^CA$/ || $which =~ /^C$/)  {
	for (my $a1 = 0; $a1 < $nat1; $a1 ++) {
	    for (my $a2 = 0; $a2 < $nat2; $a2 ++) {
		if ($type1_ref->[$a1] =~ /^$which$/ && $type2_ref->[$a2] =~ /^$which$/) {
		    $distance = euclidean_distance($x1_ref->[$a1], $y1_ref->[$a1], $z1_ref->[$a1], $x2_ref->[$a2], $y2_ref->[$a2], $z2_ref->[$a2]);
		}
	    }
	}
    }
    elsif ($which =~ /^ALL$/) {
	$distance = 1e+50;
	for (my $a1 = 0; $a1 < $nat1; $a1 ++) {
	    for (my $a2 = 0; $a2 < $nat2; $a2 ++) {
		my $dist = euclidean_distance($x1_ref->[$a1], $y1_ref->[$a1], $z1_ref->[$a1], $x2_ref->[$a2], $y2_ref->[$a2], $z2_ref->[$a2]);
		if ($dist < $distance) { $distance = $dist; }
	    }
	}
    }
    elsif ($which =~ /^NOH$/) {
	$distance = 1e+50;
	for (my $a1 = 0; $a1 < $nat1; $a1 ++) {
	    for (my $a2 = 0; $a2 < $nat2; $a2 ++) {
		if ($type1_ref->[$a1] =~ /^[^H]/ && $type2_ref->[$a2] =~ /^[^H]/) {
		    my $dist = euclidean_distance($x1_ref->[$a1], $y1_ref->[$a1], $z1_ref->[$a1], $x2_ref->[$a2], $y2_ref->[$a2], $z2_ref->[$a2]);
		    if ($dist < $distance) { $distance = $dist; }
		}
	    }
	}
    }

    return $distance;
}

sub euclidean_distance {
    my ($x1, $y1, $z1, $x2, $y2, $z2) = @_;
   
    my $distance;
    
    $distance  = ($x1-$x2) * ($x1-$x2);
    $distance += ($y1-$y2) * ($y1-$y2);
    $distance += ($z1-$z2) * ($z1-$z2);
    $distance  = sqrt($distance);

    return $distance;
}

sub plot_contact_map {
    my ($mapfile, $len, $seeplots) = @_;
    
    my $psfile = "$mapfile.ps";
    #if ($psfile =~ /\/([^\/]+)\s*$/) { $psfile = "$1"; }
    my $pdffile = $psfile;
    if ($pdffile =~ /^(\S+).ps$/) { $pdffile = "$1.pdf"; }

    my $xlabel = "alignment position";
    my $title  = "$mapfile";
    
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    
    print GP "set terminal postscript color solid 14\n";

    print GP "set output '$psfile'\n";
    print GP "unset key\n";
    print GP "set size ratio -1\n";
    print GP "set size 1,1\n";

    print GP "set nokey\n";
    print GP "set xlabel '$xlabel'\n";
    print GP "set ylabel '$xlabel'\n";
    if ($len > 0) {
	print GP "set xrange [0:$len]\n";
	print GP "set yrange [0:$len]\n";
    }

    #print GP "set title \"$title\\n\\n$key\"\n";
    print GP "set title '$title'\n";

    my $cmd = "'$mapfile' using 1:3  title '' ls 1, '$mapfile' using 3:1  title '' ls 1";
 
    print GP "plot $cmd\n";

    close (GP);

    #system ("/usr/local/bin/ps2pdf $psfile $pdffile\n"); 
    #system("/bin/rm $psfile\n");
    #if ($seeplots) { system ("open $pdffile&\n"); }
}
