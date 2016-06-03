This is a package containing scripts for computation of coupled residue pairs
in protein families. It contains:

1. dca.m - matlab script for direct coupling analysis of MSA
2. plotDCAmap.m - matlab script for contact map visualization
3. GeneratePseudobond.bash -  generate Chimera files for contact visuzalition in 3D structures
4. map2pseudobond.awk - an awk script that is a dependency of GeneratePseudobond.bash

Some of the sample input files used by the visualization scripts:

5.  BPTI_5pti_top80.pairs - contains the top 80 DCA pairs to be displayed as a contact map by plotDCAmap.m
6.  plotDCAmapVariables.mat - a matlab file containing the sample variables needed to reproduce the example described in the book chapter below. It has the top 80 pairs in a matlab vector as well as the native contact map for protein PDB 5PTI
7.  5pti_dcaTop30.map - top 30 DCA contacts to be displayed in Chimera 3D protein model of protein PDB: 5PTI. This is used by the script GeneratePseudobond.bash
8.  pseudobond_5pti_distA_pairs.dat/pseudobond_5pti_distA_pairs.cmd - are two sample output files of the script GeneratePseudobond.bash. These files can be used directly in Chimera with the PDB 5PTI. Further instructions are described in the reference below.


For more information on how to use these tools pleae refer to:

   Morcos et al. DCA for Contact Prediction. Methods in Molecular Biology, Protein Structure Prediction 3rd Edition. 2013

Please cite the following article in case of a publication resulting from the use of this software:

F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander,R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), 
Direct-coupling analysis of residue co-evolution captures native contacts across many protein families, 
Proc. Natl. Acad. Sci. 108:E1293-1301.

Further updates will be posted in  http://dca.ucsd.edu or http://dca.rice.edu or http://dca.upmc.fr

