function [nuc_mahoney,aa_mahoney,aa_variant,aa_index] = convert_nuc_index_to_aa_index_mahoney(nuc_index,new_nuc)

% new_nuc = 'A'
% nuc_index = 2749
aa_index = floor((nuc_index-740)/3);
[nuc_indices,nuc] = convert_aa_index_to_nuc_index_mahoney(aa_index);
a = nuc_indices==nuc_index;
nuc_mahoney = nuc(a);
aa_mahoney = nt2aa(nuc);
nuc_variant = nuc;
nuc_variant(find(a)) = new_nuc;
aa_variant = nt2aa(nuc_variant);