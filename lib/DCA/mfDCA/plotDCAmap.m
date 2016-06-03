function [Mat] = plotDCAmap(dca_pairs,native_pairs, pRange, ranking, mirror)
% Direct Coupling Analysis (DCA)
%
% function dca(inputfile , outputfile)
% 
% INPUTS: 
%   dca_pairs - two column vector with residue index couplets computed from DCA
%
%   native_pairs (optional) - two column vector with native contacts. Used for comparison
%
%   pRange - this is the range of the protein in the format [initial_resID final_resID]. 
%            This is needed when the Pfam domain does not cover the complete protein. 
%
%   ranking -  this is a 0/1 flag that instructs the function to color the residue 
%              pairs according to the DI rank. A color bar will appear to reflect 
%              the rank color mapping.
%   
%   mirror ? this is a 0/1 flag to decide if we want to display the contact map also
%            in the upper triangular part of the map.
%
% OUTPUT:
%
%   Mat - a matrix with the contact map combining the dca pairs, and if
%   selected, native pairs
%  
%   
% 
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied. All use is entirely at the user's own risk.
%
% Any publication resulting from applications of DCA should cite:
%
%     F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander, R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), 
%     Direct-coupling analysis of residue co-evolution captures native contacts across many protein families, 
%     Proc. Natl. Acad. Sci. 108:E1293-1301.
%
%     For more details on how to use this tool, please refer to:
%     Morcos et al. DCA for Contact Prediction. Methods in Molecular Biology, 
%     Protein Structure Prediction 3rd Edition. 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add path to Customizable heatmap that can be downloaded here:
% http://www.mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps

addpath('~/heatmaps/')  % Modify this path to comply to user system


% make sure first column is larger than second

if(~isempty(native_pairs))
     reslist = native_pairs;
     for i=1:length(reslist)
            if( reslist(i,1) > reslist(i,2) )
                d = reslist(i,2);
                reslist(i,2)= reslist(i,1);
                reslist(i,1)=d;
            end
     end
end

% make sure first column is smaller than second

 for i=1:length(dca_pairs)
        if( dca_pairs(i,1) > dca_pairs(i,2) )
            d = dca_pairs(i,2);
            dca_pairs(i,2)= dca_pairs(i,1);
            dca_pairs(i,1)=d;
        end
 end
 
 
% Visualize map
 
 Mat=zeros(pRange(1,2)- pRange(1,1)+1, pRange(1,2)- pRange(1,1)+1);
 
 rank=0;
 for i=1:length(dca_pairs)

     if (ranking ==1)
         rank= rank +1;
     else
         rank=1;
     end

     if(mirror==1 && isempty(native_pairs)) % Activate if we want maps on both sides of diagonals
         Mat(dca_pairs(i,1), dca_pairs(i,2))=rank;
     end         
      Mat( dca_pairs(i,2), dca_pairs(i,1))=rank;
 end
 
% Fill up native upper diagonal

 if( ~isempty(native_pairs))
     rank =length(dca_pairs)+10;
     for i=1:length(native_pairs)
             Mat(reslist(i,1), reslist(i,2))=rank;
     end
 end
 
  
 figure;

 if(ranking==1)
    heatmap(Mat, pRange(1,1):pRange(1,2), pRange(1,1):pRange(1,2), [], 'Colormap', 'summer')
     colorbar
     map= colormap;
     map(1,:)=[1,1,1];
     map(end,:)=[0,0,0.8];
     colormap(map);
 else
    spy(Mat, 8)
 end
 

 axis image 
 grid on
 set(gca, 'GridLineStyle', ':')
 xlabel('DCA contacts')
 if( ~isempty(native_pairs))
     title('Native contacts')
 end
 
