function [no_of_gaps_per_seq,no_of_X_per_seq] = plot_no_of_gaps_per_seq(msa,heading,flag)

[Nseq,Npos] = size(msa);
no_of_gaps_per_seq = [];
no_of_stars_per_seq= [];
no_of_B_per_seq= [];
no_of_X_per_seq= [];
for k = 1:Nseq
    no_of_gaps_per_seq(k) = length(find(msa(k,:) == '-'));    
    no_of_stars_per_seq(k) = length(find(msa(k,:) == '*'));
    no_of_B_per_seq(k) = length(find(msa(k,:) == 'B'));
    no_of_X_per_seq(k) = length(find(msa(k,:) == 'X'));
end

if flag == 1
    figure;bar(1:Nseq,no_of_gaps_per_seq)
    title(sprintf('Number of gaps per seq in %s',heading))
end