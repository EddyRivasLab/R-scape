# bash GeneratePseudobonds.bash 0  DCApairs.map filename_prefix A
# Input: offset ($1), mapfile ($2), protein prefix ($3), chain ($4)  

# offset - an offset when using C-alpha models that have a different residue indexing. 
# mapfile - is a file containing the protein residue pairs
# protein prefix - a user defined string to be used as part of the output files 
# chain - the chain ID where we want to display contacts 

#  output are two files: 
# 1) a text file with the pseudo-bonds to be read in the pseudobond reader panel in Chimera and 2) a simple Chimera script to display such bonds.  

awk -f map2pseudobond.awk -v offset=$1  -v chain=$4  $2 > pseudobond_$3_dist$4_pairs.dat

awk 'BEGIN{FS="@"}{print $1,$2}' pseudobond_$3_dist$4_pairs.dat | awk '{print "display",$1,$3}' > pseudobond_$3_dist$4_pairs.cmd


