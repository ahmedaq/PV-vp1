function [sim_mat,G,msa20]=sim_seq(algn)
% usage: [sim_mat,G]=sim_seq(algn)
%
% Returns matrices of correlation and similarity between sequences:
%
% sim_mat is a matrix of similarity which gives the fraction of common
% amino acids between pairs of sequences (counting gaps in the
% normalization).
%
% G is a matrix of covariance between sequences where each sequence of 
% length L is treated as a vector of size 20L with 0 and 1.
%
% These 2 matrices capture essentially the same information, and differ
% only by a scaling factor.
%
% Variations of these measures of correlation treating differently the gaps
% are conceivable, but they should give similar results when the alignment 
% has been properly truncated to eliminate positions with a large fraction 
% of gaps.
%
% Author: Olivier Rivoire
% 2/2010
% revised 9/2010 R. Ranganathan
%**************************************************************************
alg=lett2num(algn);

[N_seq,N_pos]=size(alg); 
msa20=zeros(N_seq,20*N_pos);
for i=1:N_pos 
    for a=1:20 
        msa20(:,20*(i-1)+a)=(alg(:,i)==a); 
    end; 
end;
G=cov(msa20',1);
sim_mat=msa20*msa20'/N_pos; 