#!/usr/bin/perl -w
#create_markplot.pl

use strict;
use Class::Struct;


use vars qw ($opt_x $opt_X $opt_y $opt_Y $opt_t);  # required if strict used
use Getopt::Std;
getopts ('x:X:y:Y:t:');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage: create_markplot.pl [options] <n> <models_string>\n\n";
        print "options:\n";
	exit;
}

my $n      = shift;
my $string = shift;

my $outfile = "figure$n.dat";

my $M = 0;
my @method = ();
while ($string) { $string =~ s/^\s*(\S+)//; $method[$M++] = $1; }

print "Methods = $M\n";
for (my $m = 0; $m < $M; $m ++) { print "$method[$m]\n"; }

my @dataset = ();
for (my $m = 0; $m < $M; $m ++) { $dataset[$m] = "$method[$m]-mono.dat$n"; }

for (my $m = 0; $m < $M; $m ++) {   
    open (OUT, ">>$outfile") || die;
    printf OUT "# $method[$m]\n";
    close(OUT);
    system("cat $dataset[$m] >> $outfile\n");
}
