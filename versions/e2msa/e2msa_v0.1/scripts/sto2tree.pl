#!/usr/bin/perl -w
#sto2tree.pl

use strict;
use Class::Struct;
use lib '/groups/eddy/home/rivase/projects/evohmm/scripts';
use lib '/Users/rivase/projects/evohmm/scripts';
use FUNCS;

use vars qw ($opt_f $opt_q $opt_v $opt_w);  # required if strict used
use Getopt::Std;
getopts ('fqvw');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  sto2tree.pl [options] <result_dir> <msafile>\n\n";
        print "options:\n";
        print "-f    :  msa is in afa format [default is stockholm format]\n";
        print "-v    :  be verbose\n";
        print "-w    :  Whelan-And-Goldman 2001 model (amino acid alignments only)\n";
        exit;
}
my $wrkdir   = shift;
my $msafile  = shift;
if (! -e $wrkdir)  { die "$wrkdir doesn't exist"; }

my $msaname = $msafile;
if ($msaname =~ /^\S+\/([^\/]+)$/) { $msaname = $1; }
if ($msaname =~ /^(\S+)\.sto$/)    { $msaname = $1; }
if ($msaname =~ /^(\S+)\.stk$/)    { $msaname = $1; }
if ($msaname =~ /^(\S+)\.fa$/)     { $msaname = $1; }
my $treefile  = "$wrkdir/$msaname.nh";

#options
my $afa          = 0;
my $do_quicktree = 0;
my $verbose      = 0;
my $wag          = 0;
if ($opt_f) { $afa          = 1; }
if ($opt_q) { $do_quicktree = 1; }
if ($opt_v) { $verbose      = 1; }
if ($opt_w) { $wag          = 1; }


# tree programs
my $fastreedir   = "/groups/eddy/home/rivase/alien-src/FastTree";
my $fastree      = "$fastreedir/src/FastTree";
my $quicktreedir = "/groups/eddy/home/rivase/alien-src/quicktree/quicktree_1.1";
my $quicktree    = "$quicktreedir/bin/quicktree";


if ($do_quicktree) {
    if (! -d $quicktreedir) { die "didn't find FastTree directory $quicktreedir"; }
    run_quicktree($quicktree, $wrkdir, $afa, $msaname, $msafile, $treefile, $verbose);
}
else {
    if (! -d $fastreedir) { die "didn't find FastTree directory $fastreedir"; }
    run_FastTree($fastree, $wrkdir, $afa, $msaname, $msafile, $treefile, $verbose);
}


# subroutines
sub run_FastTree {
    my ($fastree, $wrkdir, $is_afa, $msaname, $msafile, $treefile, $verbose) = @_;

    if (! -x $fastree)    { die "didn't find executable $fastree"; }

    # alignment needs to be in afa format
    my $infile;
    if (!$is_afa) { $infile = FUNCS::sto2afa("$wrkdir/$msaname", $msafile); }
    else          { $infile = $msafile; }

    my $cmd = "$fastree ";
    if ($wag) { $cmd .= "-wag "; }
    $cmd .= "$infile "; 

    system("echo $cmd > $treefile\n");
    system("$cmd > $treefile\n");
}

sub run_quicktree {
    my ($quicktree, $wrkdir, $is_afa, $msaname, $msafile, $treefile, $verbose) = @_;

    if (! -x $quicktree) { die "didn't find executable $quicktree"; }

    # alignment needs to be in stockholm format
    my $infile;
    if ($is_afa) { $infile = FUNCS::afa2sto("$wrkdir/$msaname", $msafile); }
    else         { $infile = $msafile; }

    my $cmd  = "$quicktree ";
    $cmd    .= "-out t "; # output is a tree in New Hampshire format
    $cmd    .= "$infile "; 

    system("$cmd > $treefile\n");

}
