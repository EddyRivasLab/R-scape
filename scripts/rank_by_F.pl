#!/usr/bin/perl -w
# 
#  Rank a R-scape output by the F measure
#

use strict;
use Class::Struct;

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
my $file = shift;

my $file3d;
if ($opt_F) { $file3d = "$opt_F"; }

my $n3d = 0;
my @fam3d;
my @name3d;
if ($file3d) { 
    parse_3dfile($file3d, \$n3d, \@fam3d, \@name3d);
    print "3D fam $n3d\n";
    for (my $f = 0; $f < $n3d; $f ++) {
	printf "%d |$fam3d[$f]| |$name3d[$f]|\n", $f+1;
    }
}

# Method Target_E-val [cov_min,conv_max] [FP | TP True Found | Sen PPV F] 
# GTp    0.05         [-9.85,91.47]     [0 | 5 6 5 | 83.33 100.00 90.91] 

my @fam;
my %fam_id;
my %fam_alen;
my %fam_nseq;
my %fam_avgsub;
my %fam_tp;
my %fam_true;
my %fam_found;
my %fam_S;
my %fam_P;
my %fam_F;
my $nf = 0;
my $fam;
my $usefam;
my $acc;
open (FILE, "$file") || die;
while(<FILE>) {
    # MSA RF00001_5S_rRNA nseq 712 (712) alen 119 (230) avgid 56.09 (55.86) nbpairs 34 (34)
    if (/^# MSA\s+(\S+)\s+nseq\s+(\S+)\s+\(.+alen\s+(\S+)\s+\(.+avgid\s+(\S+)\s+/) {
	$fam      = $1;
  	my $nseq  = $2;
	my $alen  = $3;
	my $avgid = $4;
	#print "^^ nseq $nseq alen $alen id $avgid\n";

	$acc = $fam;
	if ($acc =~ /^(RF\d\d\d\d\d)/) { $acc = $1; }
	
	$usefam = 0;
	if ($n3d == 0 || ($n3d > 0 && is_3dfam($acc, $n3d, \@fam3d)) ) { $usefam = 1; }	    

 
	if ($usefam) {
	    $fam[$nf]       = $fam;
	    $fam_nseq{$fam} = $nseq;
	    $fam_alen{$fam} = $alen;
	    $fam_id{$fam}   = $avgid;
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
	    
	    $fam_tp{$fam}    = $tp;
	    $fam_true{$fam}  = $true;
	    $fam_found{$fam} = $found;
	    
	    $fam_S{$fam}     = $sen;
	    $fam_P{$fam}     = $ppv;
	    $fam_F{$fam}     = $F;
	    $nf ++;
	}
    }
    # BPAIRS 20
    # avg substitutions per BP 15.8
    # BPAIRS expected covary 7.1
    elsif (/^#\s+avg substitutions per BP\s+(\S+)/) {
	if ($usefam) {
	    my $avgsub = $1;
	    $fam_avgsub{$fam} = $avgsub;
	}
    }
    
}
close (FILE);

my @S_order = sort { $fam_S{$b} <=> $fam_S{$a} } keys(%fam_S);
my @P_order = sort { $fam_P{$b} <=> $fam_P{$a} } keys(%fam_P);
my @F_order = sort { $fam_F{$b} <=> $fam_F{$a} } keys(%fam_F);

print "\n# nfam = $nf\n";
print "# family\t\t\t\t\tF\tSEN\tPPV\tTRUE\tFOUND\tTP\tavgsub\tavgid\talen\tnseq\n";
for (my $f = 0; $f < $nf; $f ++) {
    my $fam = $F_order[$f];
    printf "%-40s\t$fam_F{$fam}\t$fam_S{$fam}\t$fam_P{$fam}\t$fam_true{$fam}\t$fam_found{$fam}\t$fam_tp{$fam}\t$fam_avgsub{$fam}\t$fam_id{$fam}\t$fam_alen{$fam}\t$fam_nseq{$fam}\n", $fam;   
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
