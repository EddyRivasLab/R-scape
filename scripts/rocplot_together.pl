#!/usr/bin/perl -w
# rocplot_together.pl

use strict;
use Class::Struct;

use vars qw ($opt_L $opt_v);  # required if strict used
use Getopt::Std;
getopts ('L:v');


# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  rocplot_together.pl [options] <string_name> <string_type> <string_fam> \n\n";
        print "options:\n";
 	exit;
}

my $string_name = shift;
my $string_type = shift;
my $string_fam  = shift;

my @type = split(/\s+/, $string_type);
my $M = $#type+1;
print "NTYPE $M\n";
for (my $m = 0; $m < $M; $m++)
{
    $type[$m] =~ s/ //g;
    print "$type[$m]\n";
}

my @family = split(/\s+/, $string_fam);
my $F = $#family+1;
print "\nNFAM $F\n";
for (my $f = 0; $f < $F; $f++)
{
    $family[$f] =~ s/ //g;
    print "$family[$f]\n";
}


my $seeplots = 1;
my $verbose  = 0;
