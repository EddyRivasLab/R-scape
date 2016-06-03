function dca(inputfile, outputfile)
% Direct Coupling Analysis (DCA)
%
% function dca(inputfile , outputfile)
% 
% INPUTS: 
%   inputfile  - file containing the FASTA alignment
%   outputfile - file for dca results. The file is composed by N(N-1)/2 
%                (N = length of the sequences) rows and 4 columns: 
%                residue i (column 1), residue j (column 2),
%                MI(i,j) (Mutual Information between i and j), and 
%                DI(i,j) (Direct Information between i and j).
%                Note: all insert columns are removed from the alignment.
%
% SOME RELEVANT VARIABLES:
%   N        number of residues in each sequence (no insert)
%   M        number of sequences in the alignment
%   Meff     effective number of sequences after reweighting
%   q        equal to 21 (20 aminoacids + 1 gap)
%   align    M x N matrix containing the alignmnent
%   Pij_true N x N x q x q matrix containing the reweigthed frequency
%            counts.
%   Pij      N x N x q x q matrix containing the reweighted frequency 
%            counts with pseudo counts.
%   C        N(q-1) x N(q-1) matrix containing the covariance matrix.
%
%
% Copyright for this implementation: 
%             2011/12 - Andrea Pagnani and Martin Weigt
%                       andrea.pagnani@gmail.com 
%                       martin.weigt@upmc.fr
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
%     F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander, 
%     R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
%     analysis of residue co-evolution captures native contacts across 
%     many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    pseudocount_weight = 0.5; % relative weight of pseudo count   
    theta = 0.2;              % threshold for sequence id in reweighting

    
    [N,M,q,align] = return_alignment(inputfile); 
    [Pij_true,Pi_true, Meff]=Compute_True_Frequencies(align,M,N,q, theta);
    
    fprintf('### N = %d M = %d Meff = %.2f q = %d\n', N,M,Meff,q); 
    [Pij,Pi] = with_pc(Pij_true,Pi_true,pseudocount_weight,N,q);
    C = Compute_C(Pij,Pi,N,q);
    invC = inv(C);
    
    fp = fopen(outputfile,'w');
    Compute_Results(Pij, Pi,Pij_true, Pi_true, invC, N, q,fp);
    fclose(fp); 
end

function [N,M,q,Z] = return_alignment(inputfile)
% reads alignment from inputfile, removes inserts and converts into numbers

    align_full = fastaread(inputfile);
    M = length(align_full);
    ind = align_full(1).Sequence ~= '.' & align_full(1).Sequence == upper( align_full(1).Sequence );
    N = sum(ind);
    Z = zeros(M,N);

    for i=1:M
        counter = 0;
        for j=1:length(ind)
            if( ind(j) )
                counter = counter + 1;
                Z(i,counter)=letter2number( align_full(i).Sequence(j) );
            end
        end
    end
    q = max(max(Z));
end

function Compute_Results(Pij,Pi,Pij_true,Pi_true,invC, N,q, fp)
% computes and prints the mutual and direct informations

    for i=1:(N-1)
        for j=(i+1):N
            % mutual information
            [MI_true,si_true,sj_true] = calculate_mi(i,j,Pij_true,Pi_true,q);
            
            % direct information from mean-field
            W_mf = ReturnW(invC,i,j,q); 
            DI_mf_pc = bp_link(i,j,W_mf,Pi,q);
            fprintf(fp,'%d %d %g %g\n', i, j, MI_true, DI_mf_pc);
        end
    end
end

function [Pij_true,Pi_true,Meff] = Compute_True_Frequencies(align,M,N,q,theta)
% computes reweighted frequency counts

    W = ones(1,M);
    if( theta > 0.0 )
        W = (1./(1+sum(squareform(pdist(align,'hamm')<theta))));
    end
    Meff=sum(W);

    Pij_true = zeros(N,N,q,q);
    Pi_true = zeros(N,q);

    for j=1:M
        for i=1:N
            Pi_true(i,align(j,i)) = Pi_true(i,align(j,i)) + W(j);
        end
    end
    Pi_true = Pi_true/Meff;

    for l=1:M
        for i=1:N-1
            for j=i+1:N
                Pij_true(i,j,align(l,i),align(l,j)) = Pij_true(i,j,align(l,i),align(l,j)) + W(l);
                Pij_true(j,i,align(l,j),align(l,i)) = Pij_true(i,j,align(l,i),align(l,j));
            end
        end
    end
    Pij_true = Pij_true/Meff;

    scra = eye(q,q);
    for i=1:N
        for alpha=1:q
            for beta=1:q
                Pij_true(i,i,alpha,beta) = Pi_true(i,alpha) * scra(alpha,beta);
            end
        end
    end
end

function x=letter2number(a)
    switch(a)
        % full AA alphabet
        case '-'
             x=1;
        case 'A'    
            x=2;    
        case 'C'    
            x=3;
        case 'D'
            x=4;
        case 'E'  
            x=5;
        case 'F'
            x=6;
        case 'G'  
            x=7;
        case 'H'
            x=8;
        case 'I'  
            x=9;
        case 'K'
            x=10;
        case 'L'  
            x=11;
        case 'M'
            x=12;
        case 'N'  
            x=13;
        case 'P'
            x=14;
        case 'Q'
            x=15;
        case 'R'
            x=16;
        case 'S'  
            x=17;
        case 'T'
            x=18;
        case 'V'
            x=19;
        case 'W'
            x=20;
        case 'Y'
            x=21;
        otherwise
            x=1;
    end
end

function [Pij,Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight,N,q)
% adds pseudocount

    Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*ones(N,N,q,q);
    Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*ones(N,q);

    scra = eye(q);

    for i=1:N
        for alpha = 1:q
            for beta = 1:q
               Pij(i,i,alpha,beta) =  (1.-pseudocount_weight)*Pij_true(i,i,alpha,beta) + pseudocount_weight/q*scra(alpha,beta);
            end
        end
    end 
end

function C = Compute_C(Pij,Pi,N,q)
% computes correlation matrix

    C=zeros(N*(q-1),N*(q-1));
    for i=1:N
        for j=1:N
            for alpha=1:q-1
                for beta=1:q-1
                     C(mapkey(i,alpha,q),mapkey(j,beta,q)) = Pij(i,j,alpha,beta) - Pi(i,alpha)*Pi(j,beta);
                end
            end
        end
    end
end

function A=mapkey(i,alpha,q)
    A = (q-1)*(i-1)+alpha;
end

function [M,s1,s2] = calculate_mi(i,j,P2,P1,q)
% computes mutual information between columns i and j

    M = 0.;
    for alpha=1:q
        for beta = 1:q
             if( P2(i,j,alpha,beta)>0 )
                M = M + P2(i,j,alpha, beta)*log(P2(i,j, alpha, beta) / P1(i,alpha)/P1(j,beta));
            end
        end
    end

    s1=0.;
    s2=0.;
    for alpha=1:q
        if( P1(i,alpha)>0 )
            s1 = s1 - P1(i,alpha) * log(P1(i,alpha));
        end
        if( P1(j,alpha)>0 )
            s2 = s2 - P1(j,alpha) * log(P1(j,alpha));
        end
    end

end

function W=ReturnW(C,i,j,q)
% extracts coupling matrix for columns i and j

    W = ones(q,q);
    W(1:q-1,1:q-1) = exp( -C(mapkey(i,1:q-1,q),mapkey(j,1:q-1,q)) );

end

function DI = bp_link(i,j,W,P1,q)
% computes direct information

    [mu1, mu2] = compute_mu(i,j,W,P1,q);
    DI = compute_di(i,j,W, mu1,mu2,P1);
    return;
end

function [mu1,mu2] = compute_mu(i,j,W,P1,q)

    epsilon=1e-4;
    diff =1.0;
    mu1 = ones(1,q)/q;
    mu2 = ones(1,q)/q;
    pi = P1(i,:);
    pj = P1(j,:);

    while ( diff > epsilon )

        scra1 = mu2 * W';
        scra2 = mu1 * W;

        new1 = pi./scra1;
        new1 = new1/sum(new1);

        new2 = pj./scra2;
        new2 = new2/sum(new2);

        diff = max( max( abs( new1-mu1 ), abs( new2-mu2 ) ) );

        mu1 = new1;
        mu2 = new2;

    end
end

function DI = compute_di(i,j,W, mu1,mu2, Pia)
% computes direct information

    tiny = 1.0e-100;

    Pdir = W.*(mu1'*mu2);
    Pdir = Pdir / sum(sum(Pdir));

    Pfac = Pia(i,:)' * Pia(j,:);

    DI = trace( Pdir' * log( (Pdir+tiny)./(Pfac+tiny) ) );

end
