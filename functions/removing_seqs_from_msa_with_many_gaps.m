function [msa_filtered,seqs_with_many_gaps] = removing_seqs_from_msa_with_many_gaps(msa,thresh,heading)

[Nseq,Npos] = size(msa);

% [no_of_gaps_per_seq,no_of_X_per_seq] = plot_no_of_gaps_per_seq(msa,heading);
no_of_gaps_per_seq = plot_no_of_gaps_per_seq(msa,heading,0);

seqs_with_many_gaps = find(no_of_gaps_per_seq>=(thresh*Npos));

msa_filtered = msa;
msa_filtered(seqs_with_many_gaps,:) = [];
