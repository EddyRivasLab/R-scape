#!/usr/bin/perl -w
#pdb_get.pl

use strict;
use Class::Struct;
#use lib '/Users/erivas/projects/R-scape/scripts';
use lib '/Users/rivase/projects/R-scape/scripts';
use FUNCS;
#use constant GNUPLOT => '/usr/bin/gnuplot';
use constant GNUPLOT => '/opt/local/bin/gnuplot';
use LWP::Simple;

use vars qw ($opt_W $opt_D $opt_L $opt_v);  # required if strict used
use Getopt::Std;
getopts ('W:D:L:v');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  pdb_get.pl [options] <pdblist>  \n\n";
        print "options:\n";
 	exit;
}

my $pdblist = shift;

my $pdb_url = 'http://files.rcsb.org/view';

my $nf;
my @pdbname;
my @pdbf;
parse_pdblist($pdblist, \$nf, \@pdbname, \@pdbf);
print "nf = $nf\n";
for (my $f = 0; $f < $nf; $f ++) {
    printf "%d>$pdbname[$f]\n", $f+1;
    get_file_from_url($pdbname[$f], $pdbf[$f]);
}


sub parse_pdblist {
    my ($pdblist, $ret_nf, $pdbname_ref, $pdbf_ref) = @_;

    my $nf = 0;
    
    open(LIST, "$pdblist") || die;
    while (<LIST>) {
	if (/\#/) {
	}
	elsif (/^(\S+)([A_Z])\s+/) {
	    $pdbname_ref->[$nf] = $1;
	    $pdbf_ref->[$nf]    = "data/PDB/$1.$2.pdb";
	    $nf ++;
	}
	elsif (/^(\S+)\s+/) {
	    $pdbname_ref->[$nf] = $1;
	    $pdbf_ref->[$nf]    = "data/PDB/$1.pdb";
	    $nf ++;
	}
    }
    close(LIST);

    $$ret_nf = $nf;
}


 sub get_file_from_url {
     my ($pdbname, $pdbf) = @_;

     my $url = "$pdb_url/$pdbname.pdb";
     print "$url\n";
     print "$pdbf\n";
     
     open(PDB, ">$pdbf") || die;
     my $content = get($url);
     die "Couldn't get $url" unless defined $content;
     print PDB $content;
     close(PDB);
}
