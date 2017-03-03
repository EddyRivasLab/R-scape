#!/usr/bin/perl -w
#train.pl
 
use strict;
use Class::Struct;
use lib '/groups/eddy/home/rivase/projects/evohmm/scripts';
use lib '/Users/rivase/projects/evohmm/scripts';
use FUNCS;
#use constant GNUPLOT => '/usr/bin/gnuplot';
use constant GNUPLOT => '/opt/local/bin/gnuplot';
#use constant GNUPLOT => '/sw/bin/gnuplot';

use vars qw ($opt_i $opt_I $opt_f $opt_F);  # required if strict used
use Getopt::Std;
getopts ('i:I:f:F:t:');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  plot_prodist.pl [options] <file>  \n\n";
        print "options:\n";
	exit;
}

my $file = shift;
my $fmin = -1.0;
my $fmax = -1.0;
if ($opt_f) { $fmin = $opt_f; if ($fmin < 0.0)  { print "ilegal fmin $fmin\n"; die; } }
if ($opt_F) { $fmax = $opt_F; if ($fmax < 0.0)  { print "ilegal fmax $fmax\n"; die; } }
if ($fmin >= 0 && $fmax >= 0 && $fmin >= $fmax) { print "ilegal fmin > fmax \n"; die; }


my $minid = 0.;
my $maxid = 100.;
if ($opt_i) { $minid = $opt_i; if ($minid < 0.0) { print "ilegal minid $minid\n"; die; } }
if ($opt_I) { $maxid = $opt_I; if ($maxid < 0.0) { print "ilegal maxid $maxid\n"; die; } }

train($file, $minid, $maxid, $fmin, $fmax);



sub train {

    my ($file, $minid, $maxid, $fmin, $fmax) = @_;

    my $nf  = 0;
    my $nft = 0;

    my $id;
    my $match;
    my $alen;
    my $nins;
    my $itlen;
    my $ilen;
    
    my $avg_id    = 0.0;
    my $std_id    = 0.0;
    my $avg_match = 0.0;
    my $std_match = 0.0;
    my $avg_alen  = 0.0;
    my $std_alen  = 0.0;
    my $avg_nins  = 0.0;
    my $std_nins  = 0.0;
    my $avg_ilen  = 0.0;
    my $std_ilen  = 0.0;
    
    open(FILE, "$file")|| die;
    while (<FILE>) {
	if (/^\#\s+id\s+(\S+)\s+match\s+(\S+)\s+alen\s+(\S+)\s+inserts\s+(\S+)\s+ilen\s+(\S+)\s*/) {
	    $id    = $1;
	    $match = $2;
	    $alen  = $3;
	    $nins  = $4;
	    $itlen = $5;

	    $ilen = ($nins > 0)? $itlen/$nins : 0;

	    if ($id >= $minid && $id <= $maxid) { 
		$nf ++; 

		$avg_id    += $id;
		$std_id    += $id*$id;

		$avg_match += $match;
		$std_match += $match*$match;

		$avg_alen  += $alen;
		$std_alen  += $alen*$alen;

 		$avg_nins  += $nins;
		$std_nins  += $nins*$nins;

		$avg_ilen  += $ilen;
		$std_ilen  += $ilen*$ilen;
 	    }
	    $nft ++;
	}
    }
    close(FILE);

    FUNCS::calculate_averages (\$avg_id,    \$std_id,    $nf);
    FUNCS::calculate_averages (\$avg_match, \$std_match, $nf);
    FUNCS::calculate_averages (\$avg_alen,  \$std_alen,  $nf);
    FUNCS::calculate_averages (\$avg_nins,  \$std_nins,  $nf);
    FUNCS::calculate_averages (\$avg_ilen,  \$std_ilen,  $nf);
    
    printf "NF %d/%d\n", $nf, $nft;
    printf "ID    %10.2f +\- %10.2f\n", $avg_id,    $std_id;
    printf "MATCH %10.2f +\- %10.2f\n", $avg_match, $std_match;
    printf "ALEN  %10.2f +\- %10.2f\n", $avg_alen,  $std_alen;
    printf "NINS  %10.2f +\- %10.2f\n", $avg_nins,  $std_nins;
    printf "ILEN  %10.2f +\- %10.2f\n", $avg_ilen,  $std_ilen;

}


