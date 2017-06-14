package PDBFUNCS;

use strict;
use warnings;
use Class::Struct;

our $VERSION = "1.00";

use constant GNUPLOT => '/usr/local/bin/gnuplot';

struct RES => {
    char     => '$', # number of atoms
    coor     => '$', # coordenate in seqref
    nat      => '$', # number of atoms
    type     => '@', # type for eacth atom
    x        => '@', # x position for eacht atom
    y        => '@', # y position for eacht atom
    z        => '@', # z position for eacht atom
};

struct CNT => {
    i        => '$', # pdb position
    j        => '$', # 
    posi     => '$', # alignment position
    posj     => '$', # 
    chi      => '$', # character
    chj      => '$', # 
    bptype   => '$', # WWc WWt HHc HHt SSc SSt WHc WHt WSc WSt HSc HSt STACKED CONTACT BPNONE
    distance => '$', # euclidian distance
};


struct PDB2MSA => {
    pdbname   => '$', # pdb  name
    stoname   => '$',  # alignment name
    pdblen    => '$', # length of pdb sequence
    msalen    => '$', # length of pfam alignment
    
    which     => '$', # method used to defince "distance", default MIN (minimum eucledian distance between any two atoms)
    
    maxD      => '$', # maximun spacial distance (A) to define a contact
    minL      => '$', # minimum backbone distance (in pdbseq) to define a contact 
    ncnt      => '$', # number of contacts
    cnt       => '@', # a CNT structure for each contact

    nbp       => '$', # (RNA) number of basepairs (any of 12 types)
    nwc       => '$', # (RNA) number of basepairs (WWc only)
    
};


sub pdb2msa {
    my ($gnuplot, $rscapebin, $pdbfile, $stofile, $pdb2msa_ref, $maxD, $minL, $which, $isrna, $seeplots) = @_;

    $$pdb2msa_ref = PDB2MSA->new();
    my $stoname = $stofile;
    if ($stoname =~ /(PF[^\.]+)\./) { $stoname = $1; }
    if ($stoname =~ /(RF[^\.]+)\./) { $stoname = $1; }

    my $pdbname = $pdbfile;
    if ($pdbname =~ /\/([^\/]+)\s*$/) { $pdbname = $1; }

    my $pdblen;
    my $msalen;
    
    my $ncnt;
    my @cnt;
    my $nbp;
    my $nwc;
    
    contacts_from_pdbfile ($gnuplot, $rscapebin, $pdbfile, $stofile, \$msalen, \$ncnt, \@cnt, $maxD, $minL, $which, $isrna, "", "", $seeplots);
    contactlist_print(\*STDOUT, $ncnt, \@cnt, 1);
    contactlist_bpinfo($ncnt, \@cnt, \$nbp, \$nwc);
    contactlist_maxlen($ncnt, \@cnt, \$pdblen);

    $$pdb2msa_ref->{"PDB2MSA::pdbname"}  = $pdbname;
    $$pdb2msa_ref->{"PDB2MSA::stoname"}  = $stoname;
    $$pdb2msa_ref->{"PDB2MSA::pdblen"}   = $pdblen;
    $$pdb2msa_ref->{"PDB2MSA::msalen"}   = $msalen;
    
    $$pdb2msa_ref->{"PDB2MSA::ncnt"}     = $ncnt;
    $$pdb2msa_ref->{"PDB2MSA::nbp"}      = $nbp;
    $$pdb2msa_ref->{"PDB2MSA::nwc"}      = $nwc;
    $$pdb2msa_ref->{"PDB2MSA::maxD"}     = $maxD;
    $$pdb2msa_ref->{"PDB2MSA::minL"}     = $minL;
    @{$$pdb2msa_ref->{"PDB2MSA::cnt"}}   = @cnt;
}

sub contacts_from_pdbfile {
	
    my ($gnuplot, $rscapebin, $pdbfile, $stofile, $ret_msalen, $ret_ncnt_t, $cnt_t_ref, $maxD, $minL, $which, 
	$dornaview, $coorfile, $mapallfile, $seeplots) = @_;

    my $ncnt_t = 0;
    
    my $currdir = $ENV{PWD};

    if ($coorfile) { open(COORF,  ">$coorfile")  || die; }

    my $stoname = $stofile;
    if ($stoname =~ /\/([^\/]+)\.[^\/]+$/) { $stoname = $1; }
    my $pfamname = $stofile;
    if ($pfamname =~ /\/([^\/]+)\.[^\/]+$/) { $pfamname = $1; }
    
    my $hmmer       = "$rscapebin/../lib/hmmer";
    my $hmmbuild    = "$hmmer/src/hmmbuild";
    my $hmmersearch = "$hmmer/src/hmmsearch";
    
    my $easel    = "$hmmer/easel";
    my $sfetch   = "$easel/miniapps/esl-sfetch";
    my $reformat = "$easel/miniapps/esl-reformat";
    
    my $isrna = 0;
    if ($dornaview) { $isrna = 1; }
    
    my $nch = 0;
    my @chname = ();
    my @chsq = ();
    my $len;
    my $alen;
    
    my $resolution;
    my $pdbname = parse_pdb($pdbfile, \$resolution, \$nch, \@chname, \@chsq, $isrna);
    print     "# PFAM:       $stofile\n";
    print     "# PDB:        $pdbname\n";
    print     "# chains:     $nch\n";
    print     "# resolution: $resolution\n";
        
    for (my $n = 0; $n < $nch; $n ++) {
	my $map0file = "$pdbfile.chain$chname[$n].maxD$maxD.map";
	my $map1file = "$pdbfile.chain$chname[$n].maxD$maxD.$pfamname.map";
	my $mapfile  = "$currdir/$stoname.chain$chname[$n].maxD$maxD.map";
	my $corfile  = "$currdir/$stoname.chain$chname[$n].maxD$maxD.cor";
	
	print "\n chain $chname[$n]\n";
	if ($coorfile) {
	    print COORF "$corfile\n";
	}
	
	open(COR,  ">$corfile")  || die;
	open(MAP,  ">$mapfile")  || die;
	open(MAP0, ">$map0file") || die;
	open(MAP1, ">$map1file") || die;
	
	print COR  "# PDB: $pdbname\n";
	#print MAP  "# PDB: $pdbname\n"; 
	#print MAP0 "# PDB: $pdbname\n";
	#print MAP1 "# PDB: $pdbname\n";
	print COR  "# chain $chname[$n]\n";
	#print MAP  "# chain $chname[$n]\n";
	#print MAP0 "# chain $chname[$n]\n";
	#print MAP1 "# chain $chname[$n]\n";
	#print MAP1 "# maps to $pfamname\n";
	
	$len = length($chsq[$n]);
	$alen = parse_pdb_contact_map($rscapebin, $currdir, $pdbfile, $pdbname, $pfamname, \$ncnt_t, $cnt_t_ref, 
				      $stofile, $chname[$n], $chsq[$n], $which, $maxD, $minL, $isrna);
	close(COR);
	close(MAP);
	close(MAP0);
	close(MAP1);
	if ($alen == 0) { next; }

	my $xfield  = 1;
	my $yfield  = 3;
	my $xylabel = "PDB position";
	my $title   = "Contacts in pdb sequence";
	plot_contact_map($map0file, $len,  $xfield, $yfield, $title, $xylabel, $seeplots);

	$xfield  = 1;
	$yfield  = 4;
	$xylabel = "PDB position";
	$title   = "Contacts in pdb sequence that map to the alignment";
	plot_contact_map($map1file, -1,    $xfield, $yfield, $title, $xylabel, $seeplots);
	
	$xfield  = 2;
	$yfield  = 5;
	$xylabel = "Alignment position";
	$title   = "Contacts in the alignment";
	plot_contact_map($mapfile,  $alen, $xfield, $yfield, $title, $xylabel, $seeplots);
    }
    
    if ($coorfile) { close(COORF); }

    if ($mapallfile) { allcontacts_dump($mapallfile, $ncnt_t, $cnt_t_ref, $pdbname, $pfamname, $maxD, $minL); }

    if (1) {
	my $hisfile  = "$currdir/$stoname.$pdbname.chains$nch.maxD$maxD.his";
	allcontacts_histogram($hisfile, $ncnt_t, $cnt_t_ref, $pdbname, $pfamname, $maxD, $minL, 1); 
    }

    $$ret_msalen = $alen;
    $$ret_ncnt_t = $ncnt_t;
}

sub parse_pdb {
    my ($pdbfile, $ret_resolution, $ret_nch, $chname_ref, $chsq_ref, $isrna) = @_;

    my $pdbname;
    my $nch = 0;
    my $resolution;

    my $cur_chain;
    my $prv_chain = "";
    my $sq;
    my @sqlen;

    open(FILE, "$pdbfile") || die;
    while (<FILE>) {
	if (/^HEADER\s+.+\s+(\S+)\s*$/) {
	    $pdbname = $1;
	}
	elsif (/RESOLUTION.\s+(.+)$/) {
	    $resolution = $1;
	}
	elsif (/^SEQRES\s+\d+\s+(\S+)\s+(\d+)\s+(\S+.+)$/) {
	    $cur_chain   = $1;
	    my $sqlen    = $2;
	    $sq          = $3;
	    $sq =~ s/\n//g;
	    
	    if ($prv_chain =~ /^$/) {
		$chsq_ref->[$nch]   = $sq;
		$chname_ref->[$nch] = $cur_chain;
		$sqlen[$nch]        = $sqlen;
	    }
	    elsif ($cur_chain =~ /^$prv_chain$/) {
		$chsq_ref->[$nch] .= $sq;
	    }
	    else { 
		$nch ++; 
		$chsq_ref->[$nch]   = $sq;
		$chname_ref->[$nch] = $cur_chain;
		$sqlen[$nch]        = $sqlen;
	    }
	    $prv_chain = $cur_chain;
	}
    }
    close(FILE);
    $nch ++;

    for (my $n = 0; $n < $nch; $n ++) {
	my $len = sq_conversion(\$chsq_ref->[$n], $isrna);
	if ($len != $sqlen[$n]) { print "parse_pdb(): seq len is $len should be $sqlen[$n]\n"; die; }
    }
    
    $$ret_resolution = $resolution;    
    $$ret_nch        = $nch;
    return $pdbname;
}

# map[0..pdblen-1] taking values in 0..msa_alen-1
sub map_pdbsq {
    my ($rscapebin, $currdir, $stofile, $pdbname, $pfamname, $chname, $pdbsq, $map_ref) = @_;

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
    my $alen = find_pdbsq_in_pfam($rscapebin, $currdir, $stofile, $chname, $pdbname, $pdbsq, 
				  \$pfam_name, \$from_pdb, \$ali_pdb, \$from_pfam, \$ali_pfam, \$pfam_asq);
    if ($alen == 0 || $pfamname =~ //) {
	print "could not find $pdbname chain $chname in sto file\n";
	return 0;
    }

    # x (0..len-1) is coord in the pbdseq of the first aligned position
    my @ali_pdb  = split(//,$ali_pdb);
    my @ali_pfam = split(//,$ali_pfam);
    my @pfam_asq = split(//,$pfam_asq);
    my $x = $from_pdb-1;

    # y (0..pfamalen-1) is the coord in the pfam alignment of the first aligned position
    my $n = 0;
    my $y = 0;
    while ($n < $from_pfam-1) {
	if ($pfam_asq[$y] =~ /^[\.\-]$/) { 
	    #printf "# skip pfam gap $y %s \n", $pfam_asq[$y]; 
	    $y ++;
	}
	else {
	    $n ++; $y ++;
 	}
    }
    print "^^1st in pdb $x | 1st in pfamali $y\n";
    
    my $pos = 0;
    while ($pos < $alen) {
	my $pos_pdb  = uc($ali_pdb[$pos]);
	my $pos_pfam = uc($ali_pfam[$pos]);
	if ($pos_pdb =~ /^[\.\-]$/  && $pos_pfam =~ /^[\.\-]$/  ) { 
	    #printf "# double gap at pos %d\n", $pos; 
	}
	elsif ($pos_pdb =~ /^[\.\-]$/)  { 
	    #printf "# pdb gap | move pfam %d $pos_pfam\n", $y;
	    while($pfam_asq[$y] =~ /^[\.\-]$/) { 
		printf "# skip pfam gap $y %s \n", $pfam_asq[$y]; 
		$y ++; 
	    }
	    $y ++;	
	}
	elsif ($pos_pfam =~ /^[\.\-]$/)  { 
	    #printf "# pfam gap | move pdb $x $pos_pdb \n"; 
	    $x ++; 
	}
	elsif ($pos_pfam =~ /^$pos_pdb$/)  { 
	    while($pfam_asq[$y] =~ /^[\.\-]$/) { 
		#printf "# skip pfam gap $y %s \n", $pfam_asq[$y]; 
		$y ++; 
	    }
	    $map_ref->[$x] = $y;  
	    #printf "# match $x $pos_pdb | %d $pos_pfam\n", $y; 
	    $x ++; $y ++; 
	}
	else {
	    while($pfam_asq[$y] =~ /^[\.\-]$/) { 
		#printf "# skip pfam gap $y %s \n", $pfam_asq[$y]; 
		$y ++; 
	    }
	    $map_ref->[$x] = $y;  
	    #printf "# mismach $x $pos_pdb | %d $pos_pfam\n", $y; 
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
    my ($rscapebin, $currdir, $stofile, $chain, $pdbname, $pdbsq, $ret_pfamsqname, 
	$ret_from_pdb, $ret_ali_pdb, $ret_from_pfam, $ret_ali_pfam, $ret_pfam_asq) = @_;
 
    my $ali_pdb  = "";
    my $ali_pfam = "";
    my $from_pdb;
    my $from_pfam;
    my $pfam_asq = "";
    my $eval = 1;

    my $pfamsqname =  $$ret_pfamsqname;
     if ($pfamsqname eq "") {    
	print "could not find $pdbname chain $chain in sto file... trying by homology\n";
    }

    my $pdbsqfile = "$currdir/$pdbname";

    my $hmmer       = "$rscapebin/../lib/hmmer";
    my $hmmbuild    = "$hmmer/src/hmmbuild";
    my $hmmersearch = "$hmmer/src/hmmsearch";
    
    my $easel    = "$hmmer/easel";
    my $sfetch   = "$easel/miniapps/esl-sfetch";
    my $reformat = "$easel/miniapps/esl-reformat";

    open(F, ">$pdbsqfile") || die;
    print F ">$pdbname\n$pdbsq\n";
    close(F);
     
    my $hmm       = "$currdir/hmmfile";
    my $hmmout    = "$currdir/hmmout";
 
    system("$hmmbuild                       $hmm  $pdbsqfile   >  /dev/null\n");
    system("$hmmersearch -E $eval           $hmm  $stofile     >  $hmmout\n");
    system("/bin/echo $hmmersearch -E $eval $hmm  $stofile \n");
    #system("/usr/bin/more $hmmout\n");

    # take hit with best evalue
    my $pdbasq;
    my $asq;
    parse_hmmout_for_besthit($hmmout, $pdbname, \$pfamsqname, \$from_pdb, \$ali_pdb, \$from_pfam, \$ali_pfam);
    if ($pfamsqname eq "") { print "could not find best hit for chain $chain\n"; }
    
    $pfam_asq = get_asq_from_sto($reformat, $stofile, $pfamsqname, 0);

    if ($pfamsqname ne "") {
	printf "^^>$pdbname len=%d\n$pdbsq\n", length($pdbsq);
	print "^^>$pfamsqname\n$pfam_asq\n";
	print "^^$ali_pdb\n";
	print "^^$ali_pfam\n";
    }
    
    $$ret_pfamsqname = $pfamsqname;
    $$ret_ali_pdb    = $ali_pdb;
    $$ret_ali_pfam   = $ali_pfam;
    $$ret_from_pdb   = $from_pdb;
    $$ret_from_pfam  = $from_pfam;
    $$ret_pfam_asq   = $pfam_asq;

    system("/bin/rm $pdbsqfile\n");
    system("/bin/rm $hmm\n");
    system("/bin/rm $hmmout\n");

    return length($ali_pfam);
}

    
sub parse_hmmout_for_besthit {
    my ($hmmout, $pdbname, $ret_pfamsqname, $ret_from_pdb, $ret_ali_pdb, $ret_from_pfam, $ret_ali_pfam) = @_;

    my $from_pdb   = -1;
    my $to_pdb     = -1;
    my $from_pfam  = -1;
    my $to_pfam    = -1;
    my $ali_pdb    = "";
    my $ali_pfam   = "";
    my $pfamsqname = "";

    my $n = 0;
    open(HF, "$hmmout") || die;
    while(<HF>) {
	if (/#/) {
	}
	elsif (/\>\>\s+(\S+)\s*$/ && $n == 0) {
	    $pfamsqname = $1;
	    $ali_pfam   = "";
	    $ali_pdb    = "";
	    $from_pfam  = 123456789;
	    $from_pdb   = 123456789;
	    $n ++;
	}
	elsif (/\>\>\s+(\S+)/ && $n > 0) {
	    last;
	}
	elsif ($n==1 && /^\s*$pdbname\s+(\d+)\s+(\S+)\s+(\d+)\s*$/) {
	    my $i     = $1;
	    $ali_pdb .= $2;
	    $to_pdb   = $3;
	    $from_pdb = ($i < $from_pdb)? $i:$from_pdb;
	}
	elsif ($n==1 && /^\s*$pfamsqname\s+(\d+)\s+(\S+)\s+(\d+)\s*$/) {
	    my $i       = $1;
	    $ali_pfam  .= $2;
	    $to_pfam    = $3;
	    $from_pfam  = ($i < $from_pfam)? $i:$from_pfam;
	}
    }
    close(HF);

    $$ret_from_pdb   = $from_pdb;
    $$ret_from_pfam  = $from_pfam;
    $$ret_pfamsqname = ($from_pfam<0)? "":$pfamsqname;
    $$ret_ali_pdb    = $ali_pdb;
    $$ret_ali_pfam   = $ali_pfam;
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
    my ($reformat, $stofile, $name, $from) = @_;
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

sub get_first_asq_from_sto {
    my ($reformat, $stofile, $from) = @_;
    my $asq = "";

    my $afafile = "$stofile";
    if    ($afafile =~ /^(\S+).txt$/) { $afafile = "$1.afa"; }
    elsif ($afafile =~ /^(\S+).sto$/) { $afafile = "$1.afa"; }
    if (-f $afafile) { } else { system("$reformat afa $stofile > $afafile\n"); }

    my $n = 0;
    open(AFA, "$afafile") || die;
    while(<AFA>) {
	if (/^\>(\S+)/) {
	    $n ++;
	}
	elsif ($n == 1 && /^\>/) {
	    last;
	}       
	elsif ($n == 1 && /^(\S+)$/) {
	    $asq .= $1;
	}
     }
    close(AFA);
    

    return substr($asq, $from);
}


sub parse_pdb_contact_map {
    my ($rscapebin, $currdir, $pdbfile, $pdbname, $pfamname, $ret_ncnt_t, $cnt_t_ref, $stofile, $chain, $chsq, $which, $maxD, $minL, $isrna) = @_;

    my $len  = length($chsq);
    my @chsq = split(//,$chsq);
    
    printf COR  "# maxD  $maxD\n";
    printf COR  "# minL  $minL\n";
    #printf MAP  "# maxD  $maxD\n";
    #printf MAP  "# minL  $minL\n";
    #printf MAP0 "# maxD  $maxD\n";
    #printf MAP0 "# minL  $minL\n";
    #printf MAP1 "# maxD  $maxD\n";
    #printf MAP1 "# minL  $minL\n";
    
    my @map = (); # map sq to the sequence in the alignment  
    my $alen = map_pdbsq($rscapebin, $currdir, $stofile, $pdbname, $pfamname, $chain, $chsq, \@map);
    if ($alen == 0) { return $alen; }
    for (my $x = 0; $x < $len; $x ++) {
	printf COR "%d %d\n", $x+1, $map[$x]+1;
	#printf     "%d %d\n", $x+1, $map[$x]+1;
    }

    my $ncnt = 0; # the contact for this chain
    my @cnt;
    
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
    get_atoms_coord($pdbfile, $atom_offset, \@chsq, $len, $chain, \@res);
    
    for ($l1 = 0; $l1 < $len; $l1 ++) {
	$nat1  = $res[$l1]->{"RES::nat"};
	$coor1 = $res[$l1]->{"RES::coor"};
	$char1 = $res[$l1]->{"RES::char"};
	@type1 = @{$res[$l1]->{"RES::type"}};
	@x1    = @{$res[$l1]->{"RES::x"}};
	@y1    = @{$res[$l1]->{"RES::y"}};
	@z1    = @{$res[$l1]->{"RES::z"}};

	if ($char1 ne $chsq[$l1]) {
	    printf "#chain$chain;position %d/%d %d/%d: mismatched character %s should be $chsq[$coor1]\n", 
	    ($l1+1), $len, $coor1+1, $len+$atom_offset-1, $char1; 
	    die;
	}
       
	for ($l2 = $l1+1; $l2 < $len; $l2 ++) {
	    $nat2  = $res[$l2]->{"RES::nat"};
	    $coor2 = $res[$l2]->{"RES::coor"};
	    $char2 = $res[$l2]->{"RES::char"};
	    @type2 = @{$res[$l2]->{"RES::type"}};
	    @x2    = @{$res[$l2]->{"RES::x"}};
	    @y2    = @{$res[$l2]->{"RES::y"}};
	    @z2    = @{$res[$l2]->{"RES::z"}};

	    if ($char2 ne $chsq[$l2]) {
		printf "#chain$chain;position %d/%d %d/%d: mismatched character %s should be $chsq[$coor2]\n", 
		($l2+1), $len, $coor2+1, $len+$atom_offset-1, $char2; 
		die; 
	    }

	    $distance = distance($which, $char1, $nat1, \@type1, \@x1, \@y1, \@z1, $char2, $nat2, \@type2, \@x2, \@y2, \@z2);

	    if ($distance > 0) {
		if (abs($l1-$l2) >= $minL && $distance <= $maxD) {
		    printf MAP0 "%d %s %d %s %.2f\n", $l1+1, $chsq[$l1], $l2+1, $chsq[$l2], $distance;

		
		    if ($map[$l1] >= 0 && $map[$l2] >= 0) {

			$cnt[$ncnt] = CNT->new();
			$cnt[$ncnt]->{"CNT::i"}        = $l1+1;
			$cnt[$ncnt]->{"CNT::j"}        = $l2+1;
			$cnt[$ncnt]->{"CNT::posi"}     = $map[$l1]+1;
			$cnt[$ncnt]->{"CNT::posj"}     = $map[$l2]+1;
			$cnt[$ncnt]->{"CNT::chri"}     = $chsq[$l1];
			$cnt[$ncnt]->{"CNT::chrj"}     = $chsq[$l2];
			$cnt[$ncnt]->{"CNT::bptype"}   = "CONTACT";
			$cnt[$ncnt]->{"CNT::distance"} = $distance;
			
			$ncnt ++;
		    }
		}
	    }
	}
    }

    ## If RNA, run rnaview to extract the bptypes
    if ($isrna) { run_rnaview($rscapebin, $currdir, $pdbfile, $stofile, $pdbname, $pfamname, $chain, $minL, \$ncnt, \@cnt, $isrna); }
    
    # add all contacts from this chain to the list of total contacts if new
    for (my $c = 0; $c < $ncnt; $c ++) {
	addcontact($ret_ncnt_t, $cnt_t_ref, $cnt[$c]);
    }

    contactlist_print(\*STDOUT, $ncnt, \@cnt, 1);
    contactlist_print(\*MAP,    $ncnt, \@cnt, 0);
    contactlist_print(\*MAP1,   $ncnt, \@cnt, 0);
    
    return $alen;
}

sub contactlist_print {
    my ($fp, $ncnt, $cnt_ref, $addcomments) = @_;
    
    my $nbp = 0;
    my $nwc = 0;
    
    for (my $c = 0; $c < $ncnt; $c ++) {

	my $bptype = $cnt_ref->[$c]->{"CNT::bptype"};
	if    ($bptype =~ /^WWc$/)                           { $nbp ++; $nwc ++; }
	elsif ($bptype ne "STACKED" && $bptype ne "CONTACT") { $nbp ++; };
	
    	printf $fp "%d\t%d\t%s\t%d\t%d\t%s\t%s\t%.2f\n", 
		$cnt_ref->[$c]->{"CNT::i"}, $cnt_ref->[$c]->{"CNT::posi"}, $cnt_ref->[$c]->{"CNT::chri"}, 
		$cnt_ref->[$c]->{"CNT::j"}, $cnt_ref->[$c]->{"CNT::posj"}, $cnt_ref->[$c]->{"CNT::chrj"}, 
		$cnt_ref->[$c]->{"CNT::bptype"}, $cnt_ref->[$c]->{"CNT::distance"};
    }
    if ($addcomments) {
	printf $fp "# contacts: $ncnt\n";
	printf $fp "# bpairs   %d \n", $nbp;
	printf $fp "# WWc      %d \n", $nwc;
    }
}

sub contactlist_bpinfo {
    my ($ncnt, $cnt_ref, $ret_nbp, $ret_nwc) = @_;
    
    my $nbp = 0;
    my $nwc = 0;
    
    for (my $c = 0; $c < $ncnt; $c ++) {

	my $bptype = $cnt_ref->[$c]->{"CNT::bptype"};
	if    ($bptype =~ /^WWc$/)                           { $nbp ++; $nwc ++; }
	elsif ($bptype ne "STACKED" && $bptype ne "CONTACT") { $nbp ++; };
    }
    $$ret_nbp = $nbp;
    $$ret_nwc = $nwc;
}
sub contactlist_maxlen {
    my ($ncnt, $cnt_ref, $ret_maxlen) = @_;
    
    my $maxlen = 0;
     
    for (my $c = 0; $c < $ncnt; $c ++) {

	my $i = $cnt_ref->[$c]->{"CNT::i"};
	my $j = $cnt_ref->[$c]->{"CNT::j"};
	my $len = $j - $i + 1;
	if ($len > $maxlen) { $maxlen = $len; }
    }
    $$ret_maxlen = $maxlen;
}

sub found_alicoords_in_contactlist {
    my ($posi, $posj, $minL, $ncnt, $cnt_ref, $ret_type, $ret_pdbi, $ret_pdbj) = @_;

    my $found = 0;
    my $type = 14; # not a contact 

    for (my $c = 0; $c < $ncnt; $c ++) {
	if ($cnt_ref->[$c]->{"CNT::j"}-$cnt_ref->[$c]->{"CNT::i"}+1 < $minL) { continue; } # too close in pdbseq distance
	
	if ($posi == $cnt_ref->[$c]->{"CNT::posi"} && $posj == $cnt_ref->[$c]->{"CNT::posj"}) {
	    $found = 1;
	    
	    my $bptype = $cnt_ref->[$c]->{"CNT::bptype"};
	    if    ($bptype =~ /^WWc$/)     { $type = 0;  }
	    elsif ($bptype =~ /^WWt$/)     { $type = 1;  }
	    elsif ($bptype =~ /^HHc$/)     { $type = 2;  }
	    elsif ($bptype =~ /^HHt$/)     { $type = 3;  }
	    elsif ($bptype =~ /^SSc$/)     { $type = 4;  }
	    elsif ($bptype =~ /^SSt$/)     { $type = 5;  }
	    elsif ($bptype =~ /^WHc$/)     { $type = 6;  }
	    elsif ($bptype =~ /^WHt$/)     { $type = 7;  }
	    elsif ($bptype =~ /^WSc$/)     { $type = 8;  }
	    elsif ($bptype =~ /^WSt$/)     { $type = 9;  }
	    elsif ($bptype =~ /^HSc$/)     { $type = 10; }
	    elsif ($bptype =~ /^HSt$/)     { $type = 11; }
	    elsif ($bptype =~ /^STACKED$/) { $type = 12; }
	    elsif ($bptype =~ /^CONTACT$/) { $type = 13; }
	    else { print "uh? bptype = $bptype\n"; die;  }
	    
	    $$ret_type = $type;
	    $$ret_pdbi = $cnt_ref->[$c]->{"CNT::i"};
	    $$ret_pdbj = $cnt_ref->[$c]->{"CNT::j"};
	    return $found;
	}
    }
    
    $$ret_type = $type;
    return $found;
}
sub found_pdbcoords_in_contactlist {
    my ($i, $j, $minL, $ncnt, $cnt_ref, $ret_type, $ret_posi, $ret_posj) = @_;

    my $found = 0;
    my $type = 14; # not a contact 

    for (my $c = 0; $c < $ncnt; $c ++) {
	if ($cnt_ref->[$c]->{"CNT::j"}-$cnt_ref->[$c]->{"CNT::i"}+1 < $minL) { continue; } # too close in pdbseq distance
	
	if ($i == $cnt_ref->[$c]->{"CNT::i"} && $j == $cnt_ref->[$c]->{"CNT::j"}) {
	    $found = 1;
	    
	    my $bptype = $cnt_ref->[$c]->{"CNT::bptype"};
	    if    ($bptype =~ /^WWc$/)     { $type = 0;  }
	    elsif ($bptype =~ /^WWt$/)     { $type = 1;  }
	    elsif ($bptype =~ /^HHc$/)     { $type = 2;  }
	    elsif ($bptype =~ /^HHt$/)     { $type = 3;  }
	    elsif ($bptype =~ /^SSc$/)     { $type = 4;  }
	    elsif ($bptype =~ /^SSt$/)     { $type = 5;  }
	    elsif ($bptype =~ /^WHc$/)     { $type = 6;  }
	    elsif ($bptype =~ /^WHt$/)     { $type = 7;  }
	    elsif ($bptype =~ /^WSc$/)     { $type = 8;  }
	    elsif ($bptype =~ /^WSt$/)     { $type = 9;  }
	    elsif ($bptype =~ /^HSc$/)     { $type = 10; }
	    elsif ($bptype =~ /^HSt$/)     { $type = 11; }
	    elsif ($bptype =~ /^STACKED$/) { $type = 12; }
	    elsif ($bptype =~ /^CONTACT$/) { $type = 13; }
	    else { print "uh? bptype = $bptype\n"; die;  }
	    
	    $$ret_type = $type;
	    $$ret_posi = $cnt_ref->[$c]->{"CNT::posi"};
	    $$ret_posj = $cnt_ref->[$c]->{"CNT::posj"};
	    return $found;
	}
    }
    
    $$ret_type = $type;
    return $found;
}

sub addcontact {
    my ($ret_ncnt, $cnt_ref, $new_ref) = @_;

    my $ncnt = $$ret_ncnt;

    my $i        = $new_ref->{"CNT::i"};
    my $j        = $new_ref->{"CNT::j"};
    my $posi     = $new_ref->{"CNT::posi"};
    my $posj     = $new_ref->{"CNT::posj"};
    my $chri     = $new_ref->{"CNT::chri"};
    my $chrj     = $new_ref->{"CNT::chrj"};
    my $bptype   = $new_ref->{"CNT::bptype"};
    my $distance = $new_ref->{"CNT::distance"};

    my $c;
    for ($c = 0; $c < $ncnt; $c ++) {
	if ($posi == $cnt_ref->[$c]->{"CNT::posi"} && $posj == $cnt_ref->[$c]->{"CNT::posj"}) { last; }
    }

    if ($c == $ncnt) {
	$cnt_ref->[$ncnt] = CNT->new();
	$cnt_ref->[$ncnt]->{"CNT::i"}        = $i;
	$cnt_ref->[$ncnt]->{"CNT::j"}        = $j;
	$cnt_ref->[$ncnt]->{"CNT::posi"}     = $posi;
	$cnt_ref->[$ncnt]->{"CNT::posj"}     = $posj;
	$cnt_ref->[$ncnt]->{"CNT::chri"}     = $chri;
	$cnt_ref->[$ncnt]->{"CNT::chrj"}     = $chrj;
	$cnt_ref->[$ncnt]->{"CNT::bptype"}   = $bptype;
	$cnt_ref->[$ncnt]->{"CNT::distance"} = $distance;
	$$ret_ncnt ++;
    }

}

 sub allcontacts_dump {
    my ($file, $ncnt, $cnt_ref, $pdbname, $pfamname, $maxD, $minL) = @_;
    
    if ($file =~ //) { return; }

    my $nbp = 0;
    my $nwc = 0;
    open(FILE,  ">$file")  || die;
    print FILE  "# PDB   $pdbname\n"; 
    print FILE  "# MSA   $pfamname\n"; 
    print FILE  "# maxD  $maxD\n";
    print FILE  "# minL  $minL\n";
    for (my $c = 0; $c < $ncnt; $c ++) {
	if    ($cnt_ref->[$c]->{"CNT::bptype"} =~ /^WWc$/) { $nwc ++; $nbp ++; }
	elsif ($cnt_ref->[$c]->{"CNT::bptype"} =~ /^WWt$/) { $nbp ++; }
	elsif ($cnt_ref->[$c]->{"CNT::bptype"} =~ /^HH/)   { $nbp ++; }
	elsif ($cnt_ref->[$c]->{"CNT::bptype"} =~ /^SS/)   { $nbp ++; }
	elsif ($cnt_ref->[$c]->{"CNT::bptype"} =~ /^WH/)   { $nbp ++; }
	elsif ($cnt_ref->[$c]->{"CNT::bptype"} =~ /^WS/)   { $nbp ++; }
	elsif ($cnt_ref->[$c]->{"CNT::bptype"} =~ /^HS/)   { $nbp ++; }
	
	printf FILE  "%d\t%s\t%d\t%s\t%s\t%.2f\n", 
	$cnt_ref->[$c]->{"CNT::posi"}, $cnt_ref->[$c]->{"CNT::chri"}, 
	$cnt_ref->[$c]->{"CNT::posj"}, $cnt_ref->[$c]->{"CNT::chrj"}, 
	$cnt_ref->[$c]->{"CNT::bptype"}, 
	$cnt_ref->[$c]->{"CNT::distance"}; 
    }
    print FILE  "\# contacts: $ncnt\n";
    print FILE  "\# bpairs:   $nbp\n";
    print FILE  "\# wc:       $nwc\n";
    
    close (FILE);   
}

sub allcontacts_histogram {
    my ($hfile, $ncnt, $cnt_ref, $pdbname, $pfamname, $maxD, $minL, $seeplots) = @_;
    
    my @his;
    my $N = 50;
    my $k = 100;
    my $shift = 0;
    FUNCS::init_histo_array($N, $k, \@his);
    for (my $c = 0; $c < $ncnt; $c ++) {
	my $distance = $cnt_ref->[$c]->{"CNT::distance"};
	FUNCS::fill_histo_array(1, $distance, $N, $k, $shift, \@his);
    }

    FUNCS::write_histogram($N, $k, $shift, \@his, 1, $hfile, 0);
    
    my $title = "PDB: $pdbname   MSA: $pfamname   maxD: $maxD   minL: $minL\n";
    my $xlabel = "euclidian minimum distance (Angstroms)";
    my $ylabel = "number of contacts";
    my $key = "";
    my $psfile = "$hfile.ps";
    my $xleft = 0;
    my $xright = $maxD;
    my $ymax = -1;
    my $xfield = 1;
    my $yfield = 2;
    FUNCS::gnuplot_histo($hfile, $xfield, $yfield, $psfile, $title, $xlabel, $ylabel, $key, 0, $seeplots, $xleft, $xright, $ymax);
}

sub run_rnaview {
    my ($rscapebin, $currdir, $pdbfile, $stofile, $pdbname, $pfamname, $chain, $minL, $ret_ncnt, $cnt_ref, $isrna) = @_;

    my $ncnt = $$ret_ncnt;
    my $rnaviewfile = "rnaviewf";
    my $sq;
    my @map = (); # map sq to the sequence in the alignment

    my $rnaview = "$rscapebin/rnaview";

    system("$rnaview -c $chain $pdbfile > $rnaviewfile\n");
    system("/usr/bin/more $rnaviewfile\n");

    open(RF, "$rnaviewfile") || die;
    while(<RF>) {
	if (/\#\s+seq_$chain\s+(\S+)\s*$/) { 
	    $sq = $1;
	    my $alen = map_pdbsq($rscapebin, $currdir, $stofile, $pdbname, $pfamname, $chain, $sq, \@map);
	    if ($alen == 0) { return; } # cannot find this chain in the msa
	}
	elsif (/^(\d+)\s+(\d+)\s+$chain\s+\d+\s+(\S)\s+(\S)\s+\d+\s+$chain\s+(\S+)/) {
	    my $i      = $1;
	    my $j      = $2;
	    print "^^^^ i $i $3 j $j $4\n";
	    my $posi   = $map[$i-1]+1;
	    my $posj   = $map[$j-1]+1;
	    my $chri   = aa_conversion($3, $isrna);
	    my $chrj   = aa_conversion($4, $isrna);
 	    print "^^^^ i $i $3 $chri j $j $4 $chrj\n";
	    my $bptype = $5;

	    if ($posi > 0 && $posj > 0 && ($j-$i+1) >= $minL) {
		rnaview2list($posi, $chri, $posj, $chrj, $bptype, \$ncnt, $cnt_ref);
	    }
	}
    }
    close(RF);   
    
    system("/bin/rm $rnaviewfile\n");
    system("/bin/rm base_pair_statistics.out\n");
    
    $$ret_ncnt = $ncnt;
}

sub rnaview2list {
    my ($posi, $chri, $posj, $chrj, $bptype, $ret_ncnt, $cnt_ref) = @_;

    my $ncnt = $$ret_ncnt;
    my $found = 0;
    for (my $c = 0; $c < $ncnt; $c ++) {
	if ($posi == $cnt_ref->[$c]->{"CNT::posi"} && 
	    $posj == $cnt_ref->[$c]->{"CNT::posj"} ) {
	    if ($chri =~ /^$cnt_ref->[$c]->{"CNT::chri"}$/ &&
		$chrj =~ /^$cnt_ref->[$c]->{"CNT::chrj"}$/ ) {
		$found = 1;
		$cnt_ref->[$c]->{"CNT::bptype"} = $bptype;
	    }
	    else { 
		printf "i %d %s %s %d\n", $posi, $chri, $cnt_ref->[$c]->{"CNT::chri"}, $cnt_ref->[$c]->{"CNT::i"} ; 
		printf "j %d %s %s %d\n", $posj, $chrj, $cnt_ref->[$c]->{"CNT::chrj"}, $cnt_ref->[$c]->{"CNT::j"} ; 
		print "bad rnaview correspondence at $c/$ncnt\n"; 
		die; 
	    }
	}
	if ($found == 1) { last; }
    }

    if ($found == 0) {
	printf "new bpair: i %d %s j %d %s\n", $posi, $chri, $posj, $chrj; 
	
	$cnt_ref->[$ncnt] = CNT->new();
	$cnt_ref->[$ncnt]->{"CNT::i"}        = 0;
	$cnt_ref->[$ncnt]->{"CNT::j"}        = 0;
	$cnt_ref->[$ncnt]->{"CNT::posi"}     = $posi;
	$cnt_ref->[$ncnt]->{"CNT::posj"}     = $posj;
	$cnt_ref->[$ncnt]->{"CNT::chri"}     = $chri;
	$cnt_ref->[$ncnt]->{"CNT::chrj"}     = $chrj;
	$cnt_ref->[$ncnt]->{"CNT::bptype"}   = $bptype;
	$cnt_ref->[$ncnt]->{"CNT::distance"} = 0;
	$ncnt ++;
    }
    $$ret_ncnt = $ncnt;
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
    my ($ret_seq, $isrna) = @_;

    my $seq = $$ret_seq;
    my $new = "";
    my $aa;

    #print "seq\n$seq\n";
    while ($seq) {
	$seq =~ s/^(\S+)\s+//; $aa = $1;
	$new .= aa_conversion($aa, $isrna); 
    }
    #printf "new $new\n";
    $$ret_seq = $new;
    return length($new);
}

sub aa_conversion {
    my ($aa, $isrna) = @_;
    my $new;

    my $AA = uc($aa);

    if ($isrna) {
	if    ($AA =~ /^A$/)    { $new = "A"; return $new; }
	elsif ($AA =~ /^I$/)    { $new = "A"; return $new; }
	elsif ($AA =~ /^C$/)    { $new = "C"; return $new; }
	elsif ($AA =~ /^G$/)    { $new = "G"; return $new; }
	elsif ($AA =~ /^T$/)    { $new = "T"; return $new; }
	elsif ($AA =~ /^U$/)    { $new = "U"; return $new; }
	elsif ($AA =~ /^P$/)    { $new = "U"; return $new; }
	# modified residues
	elsif ($AA =~ /^1MA$/)  { $new = "A"; return $new; }
	elsif ($AA =~ /^12A$/)  { $new = "A"; return $new; }
	elsif ($AA =~ /^5MC$/)  { $new = "C"; return $new; }
	elsif ($AA =~ /^CCC$/)  { $new = "C"; return $new; }
	elsif ($AA =~ /^OMC$/)  { $new = "C"; return $new; }
	elsif ($AA =~ /^M2G$/)  { $new = "G"; return $new; }
	elsif ($AA =~ /^OMG$/)  { $new = "G"; return $new; }
	elsif ($AA =~ /^YYG$/)  { $new = "G"; return $new; }
	elsif ($AA =~ /^GDP$/)  { $new = "G"; return $new; }
	elsif ($AA =~ /^GTP$/)  { $new = "G"; return $new; }
	elsif ($AA =~ /^2MG$/)  { $new = "G"; return $new; }
	elsif ($AA =~ /^7MG$/)  { $new = "G"; return $new; }
	elsif ($AA =~ /^H2U$/)  { $new = "U"; return $new; }
	elsif ($AA =~ /^UMP$/)  { $new = "U"; return $new; }
	elsif ($AA =~ /^PSU$/)  { $new = "U"; return $new; }
	elsif ($AA =~ /^2MU$/)  { $new = "U"; return $new; }
	elsif ($AA =~ /^70U$/)  { $new = "U"; return $new; }
	elsif ($AA =~ /^5MU$/)  { $new = "U"; return $new; }
	elsif ($AA =~ /^AH2U$/) { $new = "U"; return $new; }
	elsif ($AA =~ /^BH2U$/) { $new = "U"; return $new; }
   }
    
    
    if    ($AA =~ /^ALA$/)  { $new = "A"; }
    elsif ($AA =~ /^CYS$/)  { $new = "C"; }
    elsif ($AA =~ /^ASP$/)  { $new = "D"; }
    elsif ($AA =~ /^GLU$/)  { $new = "E"; }
    elsif ($AA =~ /^BGLU$/) { $new = "E"; }
    elsif ($AA =~ /^PHE$/)  { $new = "F"; }
    elsif ($AA =~ /^GLY$/)  { $new = "G"; }
    elsif ($AA =~ /^HIS$/)  { $new = "H"; }
    elsif ($AA =~ /^ILE$/)  { $new = "I"; }
    elsif ($AA =~ /^LYS$/)  { $new = "K"; }
    elsif ($AA =~ /^LEU$/)  { $new = "L"; }
    elsif ($AA =~ /^MET$/)  { $new = "M"; }
    elsif ($AA =~ /^BMET$/) { $new = "M"; }
    elsif ($AA =~ /^MSE$/)  { $new = "M"; } #selenomethionine is incorporated as methionine
    elsif ($AA =~ /^ASN$/)  { $new = "N"; }
    elsif ($AA =~ /^PRO$/)  { $new = "P"; }
    elsif ($AA =~ /^GLN$/)  { $new = "Q"; }
    elsif ($AA =~ /^ARG$/)  { $new = "R"; }
    elsif ($AA =~ /^SER$/)  { $new = "S"; }
    elsif ($AA =~ /^SEP$/)  { $new = "S"; } # phosphoseronine
    elsif ($AA =~ /^THR$/)  { $new = "T"; }
    elsif ($AA =~ /^TPO$/)  { $new = "T"; } # phosphothereonine
    elsif ($AA =~ /^VAL$/)  { $new = "V"; }
    elsif ($AA =~ /^TRP$/)  { $new = "W"; }
    elsif ($AA =~ /^TYR$/)  { $new = "Y"; }
    elsif ($AA =~ /^DA$/)   { $new = "A"; }
    elsif ($AA =~ /^DC$/)   { $new = "C"; }
    elsif ($AA =~ /^DG$/)   { $new = "G"; }
    elsif ($AA =~ /^DT$/)   { $new = "T"; }
    elsif ($AA =~ /^\S$/)   { $new = $AA; } # if AA is already given in the 1 letter code
    else { print "aa_conversion(): uh? |$AA|\n"; die; }

    return $new;
}


# The numbering of the sequence in SEQRES does not have to agree with the numbers
# given in the ATOM line as "resSeq"
#
# Two reasons for that non-correspondence
#
# (1) Sometimes in ATOM a number is just not used
#
#       example 3hl2.pdb chain E: it goes from 16 to 18 missing 17, but 17 ISNOT a missing residue
#
#                    ATOM  14473  C5 B  U E  16      46.093 -21.393   9.929  0.50100.66           C  
#                    ATOM  14474  C6 A  U E  16      40.956 -19.338  30.396  0.50 99.67           C  
#                    ATOM  14475  C6 B  U E  16      46.175 -22.317   8.966  0.50 99.69           C  
#                    ATOM  14476  P  A  G E  18      38.337 -22.615  34.967  0.50101.62           P  
#                    ATOM  14477  P  B  G E  18      44.644 -26.224   4.396  0.50105.38           P  
#
#        the U(16) G(18) are contiguous in SEQRES.
#
# (2) There is something called "insertion residues". PDB documentations says:
#
#          Alphabet letters are commonly used for insertion code. The
#          insertion code is used when two residues have the same
#          numbering. The combination of residue numbering and
#          insertion code defines the unique residue.
#
#      example 3hl2.pdb chain E: there are 3 insertion residues:  5A, 5B, 5C
#
#                         ATOM  13879  C6 B  C E   4      28.577  -2.586  10.917  0.50 90.30           C  
#                         ATOM  13880  P  A  G E   5A     62.899 -27.100  31.234  0.50 94.87           P  
#                         .
#                         .
#                         .
#                         ATOM  13924  C4 A  G E   5A     62.453 -21.295  28.899  0.50101.28           C  
#                         ATOM  13925  C4 B  G E   5A     33.721  -4.688  10.483  0.50101.54           C  
#                         ATOM  13926  P  A  G E   5B     57.209 -23.191  31.258  0.50114.69           P  
#                         .
#                         .
#                         .
#                         ATOM  13971  C4 B  G E   5B     38.398  -6.815  12.263  0.50104.40           C  
#                         ATOM  13972  P  A  A E   5C     53.245 -19.845  29.978  0.50103.68           P  
#                         ATOM  13973  P  B  A E   5C     39.584 -11.933   9.395  0.50103.81           P  
#                         .
#                         .
#                         .
#                         ATOM  14015  C4 B  A E   5C     38.787  -9.494  15.853  0.50102.79           C  
#                         ATOM  14016  P  A  U E   6      49.823 -18.749  24.136  0.50107.74           P  
#                         ATOM  14017  P  B  U E   6      42.253 -14.348  15.233  0.50103.96           P  
#
# IN SEQRES those insertion residues are all there as CGGAU
#                                                     
#    for the 3hl2.pdb chain E the actual correspondence between
#          the SEQRES (1-90)
#          and
#          the atom resSeq numbering (1-76)
#     is as follows
#
# SEQRES_E G C C C G G A U G A U C  C  U  C  A  G  U  G  G  U  C  U  G  G  G  G  U  G  C  A  G  G
# resSeq_E 1 2 3 4 5 5 5 6 7 8 9 10 11 12 13 14 15 16 18 19 20 20 21 22 23 24 25 26 27 28 29 30 31
#                  * * *                             *      *  *
#
# SEQRES_E C U U C A A A  C  C  U  G  U  A  G  C  U  G  U  C  U  A  G  C  G  A  C  A  G  A  G  U  G  G
# resSeq_E 0 0 0 0 0 0 0  39 40 41 42 43 44 45 46 46 46 46 46 46 46 46 46 46 46 46 47 48 49 50 51 52 53
#                                               *  *  *  *  *  *  *  *  *  *  *  *
#
#          * * * * * * * -> these are documented "missing residues" but those do not affect the ordering
#
# SEQRES_E U  U  C  A  A  U  U  C  C  A  C  C  U  U  U  C  G  G  G  C  G  C  C  A
# resSeq_E 54 55 56 57 58 59 60 61 62 63 64 65 66 67 67 68 69 70 71 72 73 74 75 0 
#                                              *  *
#
#                                                                               * -> another documented "missing residue"
#
#
sub get_atoms_coord {
    my ($pdbfile, $atom_offset, $seqres_ref, $len, $chain, $res_ref) = @_;

    my $type;
    my $char;
    my $coor;
    my $nat;
    my $x;
    my $y;
    my $z;
    my $l = 0;
    my $nn = 0;
    my $recording = 0;
    my $respos_first;
    my $respos_prv = -1;
    my $id;
    my $id_prv = -1;
    
    
    for ($l = 0; $l < $len; $l ++) {
	$res_ref->[$l] = RES->new();
	$res_ref->[$l]->{"RES::nat"}  = 0;
	$res_ref->[$l]->{"RES::coor"} = -1;
	$res_ref->[$l]->{"RES::char"} = aa_conversion($seqres_ref->[$l]);
    }

    my @ismissing;
    my $remarknum = 0;
    my $from;
    my $to;
    open(FILE, "$pdbfile") || die;
    while (<FILE>) {
	if (/DBREF\s+\S+\s+$chain\s+(\d+)\s+(\d+)\s+/) {
	    $from = $1;
	    $to   = $2;
	    for (my $r = $from-1; $r < $to; $r++) {
		$ismissing[$r] = 0;
	    }
	}
    }
    close(FILE);
    
    open(FILE, "$pdbfile") || die;
    while (<FILE>) {
	if (/^REMARK\s+(\d+)\s+MISSING\s+RESIDUES/) {
	    $remarknum = $1;
	}
	elsif (/^REMARK\s+$remarknum\s+\S+\s+$chain\s+(\d+)\s*$/) {
	    my $pos = $1;
	    if ($pos < $from) {
		$ismissing[$pos-1] = 1;
		for (my $x = $pos; $x < $from; $x ++) { $ismissing[$x] = 0; }
	    }
	    elsif ($pos > $to) {  
		$ismissing[$pos-1] = 1;
		for (my $x = $to; $x < $pos-1; $x ++) { $ismissing[$x] = 0; }
	    }
	    else { $ismissing[$pos-1] = 1; }
	}
	
    }
    close(FILE);
    
    # ATOM  17182  C2'A  C E  75      91.905 -22.497  17.826  0.50 94.20           C  
    #
    #
    $l = 0;
    open(FILE, "$pdbfile") || die;
    while (<FILE>) {
	my $line = $_;

	if ($line =~ /^ATOM/ || $line =~ /^HETATM/) {
	    my $atom     = substr($line, 0,  6); if ($atom     =~ /^\s*(\S+)\s*$/) { $atom = $1; }
	    my $serial   = substr($line, 6,  7); if ($serial   =~ /^\s*(\S+)\s*$/) { $serial = $1; }
	    my $atomname = substr($line, 12, 4); if ($atomname =~ /^\s*(\S+)\s*$/) { $atomname = $1; }
	    my $altloc   = substr($line, 16, 1); if ($altloc   =~ /^\s*(\S*)\s*$/) { $altloc = $1; }
	    my $resname  = substr($line, 17, 3); if ($resname  =~ /^\s*(\S+)\s*$/) { $resname = aa_conversion($1); }
	    my $chainid  = substr($line, 21, 1); if ($chainid  =~ /^\s*(\S*)\s*$/) { $chainid = $1; }
	    my $respos   = substr($line, 22, 4); if ($respos   =~ /^\s*(\S+)\s*$/) { $respos = $1; }
	    my $icode    = substr($line, 26, 1); if ($icode    =~ /^\s*(\S)\s*$/)  { $icode = $1; }
	    my $x        = substr($line, 30, 8); if ($x        =~ /^\s*(\S+)\s*$/) { $x = $1; }
	    my $y        = substr($line, 38, 8); if ($y        =~ /^\s*(\S+)\s*$/) { $y = $1; }
	    my $z        = substr($line, 46, 8); if ($z        =~ /^\s*(\S+)\s*$/) { $z = $1; }

	    # Look for the target chain
	    if ($chainid ne $chain) { next; }
	    
	    # Ignore HOH WAT MN atoms
	    if ($atom =~ /^HETATM$/ && ( $resname =~ /^HOH$/ || $resname =~ /^WAT$/ || $resname =~ /^MN$/) ) { next; }
	    
	    # An atom to record
	    $id = "$respos"."$icode";
	    if ($recording == 0) { $respos_first = $respos; }
	    $recording = 1;
	    
	    if ($respos - $respos_first > $len) { next; }
	    if ($respos < $respos_first)        { next; }

	    if ($nn == 0 || $id =~ /^$id_prv$/) {
		$nat = $res_ref->[$l]->{"RES::nat"};
		${$res_ref->[$l]->{"RES::type"}}[$nat] = $atomname;
		${$res_ref->[$l]->{"RES::x"}}[$nat]    = $x;
		${$res_ref->[$l]->{"RES::y"}}[$nat]    = $y;
		${$res_ref->[$l]->{"RES::z"}}[$nat]    = $z;
		$res_ref->[$l]->{"RES::nat"} ++;
	    }
	    else {
		$l ++;
		if ($l > 0 && $respos != $respos_prv+1) {
		    for (my $p = $respos_prv+1; $p < $respos; $p ++) {
			if ($ismissing[$p-1]) { $l ++; }
		    }

		}
		
	    }
	    if ($l >= $len) { printf("l %d >= len %d\n", $l, $len); die; }
	    $res_ref->[$l]->{"RES::coor"} = $respos;
 
	    #printf "^^at %d> |$atom|\t|$serial|\t|$atomname|\t|$altloc|\t|$resname|$icode|\t|$chainid|\t|$respos|\t|$icode|\t|$x|\t|$y|\t|$z|\n",  $l+1;
	    if ($l < 0)     { print "bad lenght\n"; die; }
 	    if ($l >= $len) { print "too long?\n";  die; }

	    $id_prv     = $id;
	    $respos_prv = $respos;
	    $nn ++;
	}
	# How to terminate chain
	elsif ($recording && $line =~ /TER/) { last; }
   }
    close(FILE);
    
    for ($l = 0; $l < $len; $l ++) {
	$nat  = $res_ref->[$l]->{"RES::nat"};	    
	$coor = $res_ref->[$l]->{"RES::coor"};	    
	$char = $res_ref->[$l]->{"RES::char"};	    
	#if ($nat == 0) { print "#res $l has not atoms\n"; }
	if (0) { printf "res %d coor %d char %s nat %d\n", $l+1, $coor, $char, $nat; }
    }
}


sub distance {
    my ($which, $char1, $nat1, $type1_ref, $x1_ref, $y1_ref, $z1_ref, $char2, $nat2, $type2_ref, $x2_ref, $y2_ref, $z2_ref) = @_;

    my $distance = 1e+50;
    if ($nat1 == 0 || $nat2 == 0) { return $distance; }
    
    if  ($which =~ /^CA$/ || $which =~ /^C$/ || $which eq "C1'")  {
	for (my $a1 = 0; $a1 < $nat1; $a1 ++) {
	    if ($type1_ref->[$a1] =~ /^$which$/) {
		for (my $a2 = 0; $a2 < $nat2; $a2 ++) {
		    if ($type2_ref->[$a2] =~ /^$which$/) {
			$distance = euclidean_distance($x1_ref->[$a1], $y1_ref->[$a1], $z1_ref->[$a1], $x2_ref->[$a2], $y2_ref->[$a2], $z2_ref->[$a2]);
		    }
		}
	    }
	}
    }
    elsif  ($which =~ /^CB$/)  { # use CA if a GLY
	for (my $a1 = 0; $a1 < $nat1; $a1 ++) {
	    if ($type1_ref->[$a1] =~ /^$which$/ || ($type1_ref->[$a1] =~ /^CA$/ &&  $char1 =~ /^G$/) ) {
		for (my $a2 = 0; $a2 < $nat2; $a2 ++) {
		    if ($type2_ref->[$a2] =~ /^$which$/ || ($type2_ref->[$a2] =~ /^CA$/ &&  $char2 =~ /^G$/) ) {
			$distance = euclidean_distance($x1_ref->[$a1], $y1_ref->[$a1], $z1_ref->[$a1], $x2_ref->[$a2], $y2_ref->[$a2], $z2_ref->[$a2]);
		    }
		}
	    }
	}
    }
    elsif ($which =~ /^MIN$/) { # minimum atom distance
	$distance = 1e+50;
	for (my $a1 = 0; $a1 < $nat1; $a1 ++) {
	    for (my $a2 = 0; $a2 < $nat2; $a2 ++) {
		my $dist = euclidean_distance($x1_ref->[$a1], $y1_ref->[$a1], $z1_ref->[$a1], $x2_ref->[$a2], $y2_ref->[$a2], $z2_ref->[$a2]);
		if ($dist < $distance) { $distance = $dist; }
	    }
	}
    }
    elsif ($which =~ /^AVG$/) { # average
	$distance = 0;
	my $n = 0;
	for (my $a1 = 0; $a1 < $nat1; $a1 ++) {
	    for (my $a2 = 0; $a2 < $nat2; $a2 ++) {
		my $dist = euclidean_distance($x1_ref->[$a1], $y1_ref->[$a1], $z1_ref->[$a1], $x2_ref->[$a2], $y2_ref->[$a2], $z2_ref->[$a2]);
		if  ($dist > 0) {
		    $distance += $dist;
		    $n ++;
		}
	    }
	}
	$distance /= $n;
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
    my ($mapfile, $len, $xfield, $yfield, $title, $xylabel, $seeplots) = @_;

    my $key = $mapfile;
    if ($key =~ /([^\/]+)\s*$/) { $key = $1; }
    
    my $psfile = "$mapfile.ps";
    if ($psfile =~ /\/([^\/]+)\s*$/) { $psfile = "$1"; }
    #my $pdffile = $psfile;
    #if ($pdffile =~ /^(\S+).ps$/) { $pdffile = "$1.pdf"; }
    
    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    
    print GP "set terminal postscript color solid 14\n";

    print GP "set output '$psfile'\n";
    #print GP "unset key\n";
    print GP "set size ratio -1\n";
    print GP "set size 1,1\n";

    print GP "set nokey\n";
    print GP "set xlabel '$xylabel'\n";
    print GP "set ylabel '$xylabel'\n";
    if ($len > 0) {
	print GP "set xrange [0:$len]\n";
	print GP "set yrange [0:$len]\n";
    }

    print GP "set title \"$title\\n\\n$key\"\n";
    #print GP "set title '$title'\n";

    my $cmd = "x title '' with lines, '$mapfile' using $xfield:$yfield  title '' ls 3, '$mapfile' using $yfield:$xfield  title '' ls 3";
 
    print GP "plot $cmd\n";
    close (GP);

    #system ("/usr/local/bin/ps2pdf $psfile $pdffile\n"); 
    #system("/bin/rm $psfile\n");
    if ($seeplots) { system ("open $psfile&\n"); }
}

1
