#!/usr/bin/perl -w
#pdb_parse.pl

use strict;
use Class::Struct;

use vars qw ($opt_W $opt_D $opt_L $opt_v);  # required if strict used
use Getopt::Std;
getopts ('W:D:L:v');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  pdb_parse.pl [options] <pdbfile> <stofile> \n\n";
        print "options:\n";
 	exit;
}

my $pdbfile = shift;
my $stofile = shift;
my $gnuplot = shift;
use constant GNUPLOT => '$gnuplot';

my $currdir = $ENV{PWD};

my $stoname = $stofile;
if ($stoname =~ /\/([^\/]+)\.\S+$/) { $stoname = $1; }
my $pfamname = $stofile;
if ($pfamname =~ /\/([^\/]+)\.\S+$/) { $pfamname = $1; }

my $hmmer       = "lib/hmmer";
my $hmmbuild    = "$hmmer/src/programs/hmmbuild";
my $hmmersearch = "$hmmer/src/programs/hmmsearch";

my $easel    = "$hmmer/lib/easel";
my $sfetch   = "$easel/miniapps/esl-sfetch";
my $reformat = "$easel/miniapps/esl-reformat";

my $maxD = 8;
if ($opt_D) { $maxD = $opt_D; }

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

my $seeplots = 1;

my $pdbname = parse_pdb($pdbfile, \$nch, \@chname, \@sq);
#print     "# PFAM: $stofile\n";
#print     "# PDB:  $pdbname\n# $nch chains\n";

for (my $n = 0; $n < $nch; $n ++) {
    my $map0file = "$pdbfile.chain$chname[$n].maxD$maxD.map";
    my $map1file = "$pdbfile.chain$chname[$n].maxD$maxD.$pfamname.map";
    my $mapfile  = "$currdir/$stoname.chain$chname[$n].maxD$maxD.map";
    my $corfile  = "$currdir/$stoname.chain$chname[$n].maxD$maxD.cor";

    print "$mapfile\n";
    
    open(COR,  ">$corfile")  || die;
    open(MAP,  ">$mapfile")  || die;
    open(MAP0, ">$map0file") || die;
    open(MAP1, ">$map1file") || die;
    
    #print      "# PDB: $pdbname\n";
    print COR  "# PDB: $pdbname\n";
    print MAP  "# PDB: $pdbname\n";
    print MAP0 "# PDB: $pdbname\n";
    print MAP1 "# PDB: $pdbname\n";
    #print      "# chain $chname[$n]\n# $sq[$n]\n";
    print COR  "# chain $chname[$n]\n";
    print MAP  "# chain $chname[$n]\n";
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
	elsif (/^SEQRES\s+\d+\s+(\S+)\s+(\d+)\s+(\S+.+\S+.+)$/) {
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

sub map_sq {
    my ($stofile, $pdbname, $chname, $sq, $map_ref) = @_;

    my $len = length($sq);
    for (my $l = 0; $l < $len; $l ++) {
	$map_ref->[$l] = -1;
    }

    my $name = "";
    my $asq = "";
    my $i;
    my $j;
    my $exactsq = 1;
    open(STO, "$stofile") || die;
    while (<STO>) {
	if (/^\#\=GS\s+(\S+)\s+DR PDB\;\s+\S*$pdbname\S*\s+$chname\;\s+(\d+)\-(\d+)\;/) {
	    $name = $1;
	    $i    = $2;
	    $j    = $3;
	    #print "# $name $i $j\n";
	}
	elsif (/^\#\=GS\s+(\S+)\s+AC\s+\S*$pdbname\S*\s*$/) {
	    $name = $1;
	    $i    = 1;
	    $j    = -1;
	    #print "# $name\n";
	}
	elsif (/^$name\s+(\S+)\s*$/) {
	    $asq .= uc $1;
	}
    }
    close(STO);
    if ($name =~ /^$/) {
	$exactsq = 0;
	print "could not find $pdbname name chain $chname in sto file.. trying by homology\n"; 
	return 0;
	if (find_pdbsq_in_pfam($stofile, $sq, \$name, \$asq) == 0) {
	    print "could not find $pdbname in sto file\n";
	    return 0;
	}
    }
    
    my @sq  = split(//,$sq);
    my @asq = split(//,$asq);
    my $x = first_pos_in_pdbsq($sq, $asq);    
    my $y = 0;
    while ($x < $len && $y < length($asq)) {
	if ($asq[$y] =~ /^\.$/  || $asq[$y] =~ /^\-$/) { 
	    #print "# gap $y $asq[$y]\n"; 
	    $y ++; 
	}
	elsif ($exactsq && $sq[$x] =~ /^$asq[$y]$/)  { 
	    $map_ref->[$x] = $y;  
	    #print "# match $x $sq[$x] | $y $asq[$y]\n"; 
	    $x ++; $y ++; 
	}
	elsif ($exactsq == 0)               { 
	    $map_ref->[$x] = $y;  
	    #print "# match $x $sq[$x] | $y $asq[$y]\n"; 
	    $x ++; $y ++;
	}

	else {
	    # print "uh? $x $sq[$x] | $y $asq[$y]\n"; die;
	    $map_ref->[$x] = $y;  
	    #print "# mismach $x $sq[$x] | $y $asq[$y]\n"; 
	    $x ++; $y ++;
	} 
    }
    
    return length($asq);
}

sub find_pdbsq_in_pfam {
    my ($stofile, $pdbsq, $ret_name, $ret_pfam_asq) = @_;
 
    my $pfam_asq = "";
    my $name     =  "";
    my $from_pdb;
    my $to_pdb;
    my $from;
    my $to;

    my $pdbsqfile = "$currdir/pdbsqfile";
    
    open(F, ">$pdbsqfile") || die;
    print F ">pdbsq\n$pdbsq\n";
    close(F);
     
    my $hmm    = "$currdir/hmmfile";
    my $hmmout = "$currdir/hmmout";
 
    system("$hmmbuild             $hmm $pdbsqfile >  $hmmout\n");
    system("$hmmersearch -E 1e-15 $hmm $stofile   >  $hmmout\n");
    
    system("more $hmmout\n");

    # take hit with best evalue
    my $pdbasq;
    my $asq;
    parse_hmmout_for_besthit($hmmout, \$name, \$from_pdb, \$to_pdb, \$from, \$to, \$pdbasq, \$asq);
  
    $pfam_asq = get_asq_from_sto($stofile, $name, $from, $to);
    #print "# pfamseq $name $from-$to\n# $pfam_asq";

    
    $$ret_pfam_asq = $pfam_asq;

    system("rm $pdbsqfile\n");
    system("rm $hmm\n");
    system("rm $hmmout\n");

    return length($pfam_asq);
}

sub parse_hmmout_for_besthit {
    my ($hmmout, $ret_name, $ret_from_pdb, $ret_to_pdb, $ret_from, $ret_to, $ret_pdbasq, $ret_asq) = @_;

    my $from_pdb;
    my $to_pdb;
    my $from = -1;
    my $to = -1;
    my $pdbasq = "";
    my $asq    = "";
    my $name   = "";

    $$ret_from_pdb = $from_pdb;
    $$ret_to_pdb   = $to_pdb;
    $$ret_from     = $from;
    $$ret_to       = $to;
    $$ret_name     = $name;
    $$ret_pdbasq   = $pdbasq;
    $$ret_asq      = $asq;
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
    my ($stofile, $name, $from, $to) = @_;
    my $asq = "";

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
	elsif ($found && /^(\S+)$/) {
	    $asq .= $1;
	}
 	elsif ($found == 1 && /^\>/) {
	    $found = 0;
	    last;
	}       
    }
    close(AFA);
    
    return $asq;
}

sub first_pos_in_pdbsq {
    my ($pdbsq, $asq) = @_;
    
    my $pos = 0;
    my $sq = $asq;
    $sq =~ s/\.//g;
    $sq =~ s/\-//g;

    #print "# $pdbsq\n# $asq\n# $sq\n";
 
    my $pdbsqfile = "$currdir/pdbsqfile";
    my $sqfile    = "$currdir/sqfile";
    
    open(F, ">$pdbsqfile") || die;
    print F ">pdbsq\n$pdbsq\n";
    close(F);
    open(F, ">$sqfile") || die;
    print F ">sq\n$sq\n";
    close(F);
    
    my $hmm    = "$currdir/hmmfile";
    my $hmmout = "$currdir/hmmout";
    my $hmmtbl = "$currdir/hmmtbl";

    system("$hmmbuild                        $hmm $pdbsqfile >  $hmmout\n");
    system("$hmmersearch --domtblout $hmmtbl $hmm $sqfile    >> $hmmout\n");

    #system("more $hmmout\n");
    open(H, "$hmmtbl") || die;
    while(<H>) {
	if (/^\#/) {
	    }
	elsif (/^\S+\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\s+\S+\s+\S+\s+\S+\s+\d+\s+\d+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\d+)\s+\d+/) {
	    $pos = $1-1;
	}
    }
    close(H);

    system("/bin/rm $pdbsqfile\n");
    system("/bin/rm $sqfile\n");
    system("/bin/rm $hmm\n");
    system("/bin/rm $hmmout\n");
    system("/bin/rm $hmmtbl\n");
    
    return $pos;
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
    printf MAP0 "# maxD  $maxD\n";
    printf MAP0 "# minL  $minL\n";
    printf MAP1 "# maxD  $maxD\n";
    printf MAP1 "# minL  $minL\n";

    my @map = (); # map sq to the sequence in the alignment  
    my $alen = map_sq($stofile, $pdbname, $chain, $sq, \@map);
    if ($alen == 0) { return $alen; }
    for (my $x = 0; $x < $len; $x ++) {
	printf COR "%d %d\n", $x+1, $map[$x]+1;
    }
 
    my $nct = 0;
    my $nc = 0;
    my $nat1;
    my $nat2;
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
    #print "#atom offset $atom_offset\n";
    
    for ($l1 = 0; $l1 < $len; $l1 ++) {
       	get_atom_coords($pdbfile, $l1+$atom_offset, $chain, \$nat1, \$char1, \@type1, \@x1, \@y1, \@z1);
	if ($nat1 == 0) { #print "atom $l1 not present\n"; 
	}
	else {
	    if (aa_conversion($char1) =~ /^$sq[$l1]$/) {
		#printf "\n#chain$chain;position %d/%d: character $char1\n", $l1+$atom_offset, $len;
	    } 
	    else { printf "#chain$chain;position %d/%d: mismatched character %s should be $sq[$l1]\n", $l1+$atom_offset, $len, aa_conversion($char1); die; }
	}
	
	for ($l2 = $l1+1; $l2 < $len; $l2 ++) {
	    get_atom_coords($pdbfile, $l2+$atom_offset, $chain, \$nat2, \$char2, \@type2, \@x2, \@y2, \@z2);
	    if ($nat2 == 0) { #print "#atom $l2 not present\n"; 
	    }
	    else {
		if (aa_conversion($char2) =~ /^$sq[$l2]$/) {
		    #printf "#chain$chain;position %d/%d: $char2\n", $l2+$atom_offset, $len;
		} 
		else { printf "#chain$chain;position %d/%d: mismatched character %s should be $sq[$l2]\n", $l2+1, $len, aa_conversion($char2); die; }
	    }

	    $distance = distance($which, $nat1, \@type1, \@x1, \@y1, \@z1, $nat2, \@type2, \@x2, \@y2, \@z2);

	    if ($distance > 0) {
		if (0 && $map[$l1] >= 0 && $map[$l2] >= 0) {
		    printf "^^^%d | %d %d %s | %d %d %s | %.2f\n", abs($l1-$l2), $l1+1, $map[$l1]+1, $sq[$l1], $l2+1, $map[$l2]+1, $sq[$l2], $distance;
		}  
		
		if (abs($l1-$l2) >= $minL && $distance <= $maxD) {
		    printf MAP0 "%d %s %d %s %.2f\n", $l1+1, $sq[$l1], $l2+1, $sq[$l2], $distance;
		    $nct ++;
		    
		    if ($map[$l1] >= 0 && $map[$l2] >= 0) {
			$nc ++;
			#printf      "%d %d %s | %d %d %s | %.2f\n", $l1+1, $map[$l1]+1, $sq[$l1], $l2+1, $map[$l2]+1, $sq[$l2], $distance; 
			printf MAP  "%d %s %d %s %.2f\n", $map[$l1]+1, $sq[$l1], $map[$l2]+1, $sq[$l2], $distance; 
			printf MAP1 "%d %s %d %s %.2f\n", $l1+1, $sq[$l1], $l2+1, $sq[$l2], $distance; 
		    }
		}
	    }
	}
    }
    #print      "\# contacts: $nc\n";
    print MAP  "\# contacts: $nc\n";
    print MAP0 "\# contacts: $nct\n";
    print MAP1 "\# contacts: $nc\n";

    return $alen;
}

sub atom_offset {
    my ($pdbfile, $chain) = @_;

    my $atom_offset = 1;
    open(FILE, "$pdbfile") || die;
    while (<FILE>) {
	if (/^DBREF\s+\S+\s+$chain\s+(\S+)\s+\S+\s+PDB/) {
	    $atom_offset = $1;
	}
	elsif (/^DBREF\s+\S+\s+$chain\s+(\S+)\s+\S+\s+UNP/) {
	    $atom_offset = $1;
	}
	elsif (/SEQADV\s+\S+\s+\S+\s+$chain\s+(\S+)\s+/) {
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
    #print "seq\n$new\n";
    $$ret_seq = $new;
    return length($new);
}

sub aa_conversion {
    my ($aa) = @_;
    my $new;
    if    ($aa =~ /^ALA$/) { $new = "A"; }
    elsif ($aa =~ /^CYS$/) { $new = "C"; }
    elsif ($aa =~ /^ASP$/) { $new = "D"; }
    elsif ($aa =~ /^GLU$/) { $new = "E"; }
    elsif ($aa =~ /^PHE$/) { $new = "F"; }
    elsif ($aa =~ /^GLY$/) { $new = "G"; }
    elsif ($aa =~ /^HIS$/) { $new = "H"; }
    elsif ($aa =~ /^ILE$/) { $new = "I"; }
    elsif ($aa =~ /^LYS$/) { $new = "K"; }
    elsif ($aa =~ /^LEU$/) { $new = "L"; }
    elsif ($aa =~ /^MET$/) { $new = "M"; }
    elsif ($aa =~ /^MSE$/) { $new = "M"; } #selenomethionine is incorporated as methionine
    elsif ($aa =~ /^ASN$/) { $new = "N"; }
    elsif ($aa =~ /^PRO$/) { $new = "P"; }
    elsif ($aa =~ /^GLN$/) { $new = "Q"; }
    elsif ($aa =~ /^ARG$/) { $new = "R"; }
    elsif ($aa =~ /^SER$/) { $new = "S"; }
    elsif ($aa =~ /^SEP$/) { $new = "S"; } # phosphoseronine
    elsif ($aa =~ /^THR$/) { $new = "T"; }
    elsif ($aa =~ /^TPO$/) { $new = "T"; } # phosphothereonine
    elsif ($aa =~ /^VAL$/) { $new = "V"; }
    elsif ($aa =~ /^TRP$/) { $new = "W"; }
    elsif ($aa =~ /^TYR$/) { $new = "Y"; }
    elsif ($aa =~ /^DA$/)  { $new = "A"; }
    elsif ($aa =~ /^DC$/)  { $new = "C"; }
    elsif ($aa =~ /^DG$/)  { $new = "G"; }
    elsif ($aa =~ /^DT$/)  { $new = "T"; }
    elsif ($aa =~ /^A$/)   { $new = "A"; }
    elsif ($aa =~ /^C$/)   { $new = "C"; }
    elsif ($aa =~ /^G$/)   { $new = "G"; }
    elsif ($aa =~ /^T$/)   { $new = "T"; }
    elsif ($aa =~ /^U$/)   { $new = "U"; }
    elsif ($aa =~ /^UMP$/) { $new = "U"; }
    else { print "aa_conversion(): uh? $aa\n"; die; }

    return $new;
}


sub get_atom_coords {
    my ($pdbfile, $nn, $chain, $ret_nat, $ret_char, $type_ref, $x_ref, $y_ref, $z_ref) = @_;
    
    my $nat = 0;
    my $type;
    my $char = "9";
    my $x;
    my $y;
    my $z;
    
    my $thisn = 0;
    open(FILE, "$pdbfile") || die;
    while (<FILE>) {
	if (/^ATOM\s+\d+\s+N\s+\S*(\S\S\S)\s+$chain\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/) {
	    $type  = "N";
	    $char  = $1;
	    $x     = $3;
	    $y     = $4;
	    $z     = $5;
	    $thisn = $2;
	    if ($thisn == $nn) { $type_ref->[$nat] = $type; $x_ref->[$nat] = $x; $y_ref->[$nat] = $y; $z_ref->[$nat] = $z; $nat ++; $$ret_char = $char; }
	}
	elsif (/^ATOM\s+\d+\s+(\S+)\s*\S\S\S\s+$chain\s+$thisn\s+(\S+)\s+(\S+)\s+(\S+)\s+/) {
	    $type = $1;
	    $x    = $2;
	    $y    = $3;
	    $z    = $4;
	    if ($thisn == $nn) { $type_ref->[$nat] = $type; $x_ref->[$nat] = $x; $y_ref->[$nat] = $y; $z_ref->[$nat] = $z; $nat ++; }
	}
    }

    $$ret_nat = $nat;
    if (0) { printf "res %d char %s nat %d\n", $nn, $$ret_char, $nat; }
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
 
    #print    "plot $cmd\n";
    print GP "plot $cmd\n";

    close (GP);

    system ("/usr/local/bin/ps2pdf $psfile $pdffile\n"); 
    system("/bin/rm $psfile\n");
    if ($seeplots) { system ("open $pdffile&\n"); }
}
