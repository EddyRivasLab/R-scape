#!/usr/bin/perl -w
# msa2sto.pl
# msa = the format required by gremlin
# name1 aseq1
# name2 aseq2


use strict;
use Class::Struct;
use lib '/groups/eddy/home/rivase/projects/evohmm/scripts';
use lib '/Users/rivase/projects/evohmm/scripts';
use FUNCS;

use vars qw ($opt_v );  # required if strict used
use Getopt::Std;
getopts ('v');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  msa2sto.pl [options] <msa_file> \n\n";
        print "options:\n";
 	print "-v    :  be verbose\n";
 	exit;
}
my $msafile = shift;
my $afafile = "$msafile";
my $stofile = "$msafile";
if ($afafile =~ /^(\S+).msa$/) { $afafile = "$1.afa"; }
if ($stofile =~ /^(\S+).msa$/) { $stofile = "$1.sto"; }

my $eslreformat = "/n/eddyfs01/home/erivas/src/src/mysource/versions/rscape/rscape_v0.6.1/lib/hmmer/easel/miniapps/esl-reformat ";

my $name = "";
my $asq = "";

open (AFA, ">$afafile") || die;
open (MSA, "$msafile") || die;
while(<MSA>) {
    if (/^(\S+)\s+(\S+)$/) {
	$name = $1;
	$asq  = $2;

	print AFA ">$name\n$asq\n";
    }
}
close (MSA);

close (AFA);

system("$eslreformat STOCKHOLM $afafile > $stofile\n");
