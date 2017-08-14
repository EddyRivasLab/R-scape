#!/usr/bin/perl -w
#afa2msa.pl
# the format required by gremlin
# name1 aseq1
# name2 aseq2


use strict;
use Class::Struct;

use vars qw ($opt_v );  # required if strict used
use Getopt::Std;
getopts ('v');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  afa2msa.pl [options] <afa_file> \n\n";
        print "options:\n";
 	print "-v    :  be verbose\n";
 	exit;
}
my $afafile = shift;
my $msafile = "$afafile";
if ($msafile =~ /^(\S+).afa$/) { $msafile = "$1.msa"; }
if ($msafile =~ /^(\S+).fas$/) { $msafile = "$1.msa"; }

my $nsq = 0;
my $name = "";
my $curname;
my $asq = "";

open (MSA, ">$msafile") || die;
open (AFA, "$afafile") || die;
while(<AFA>) {
    if (/^\>(\S+)\s*/) {
	$curname = $1;

	if ($nsq > 0) {
	    print MSA "$name $asq\n";
	}
	
	$asq = "";
	$name = $curname;
	$nsq ++;
    }
    elsif(/^(\S+)$/) {
	$asq .= "$1";
    }
    
}
close (AFA);

#last case
print MSA "$name $asq\n";

close (MSA);
