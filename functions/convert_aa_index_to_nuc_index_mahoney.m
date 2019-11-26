function [index_nuc,nuc] = convert_aa_index_to_nuc_index_mahoney(aa_index)

%range of aa_index = 1-2209

[h_nuc,s_nuc] = fastaread('mahoney_v01149_T2133C.fa');
[h_aa,s_aa] = fastaread('Mahoney_strain_Human_poliovirus1.fasta');

s_nuc_frame_shift = s_nuc(2:end);
s_nuc_to_aa = nt2aa(s_nuc_frame_shift); %frame shift by 1

[score,alignment,index] = swalign(s_nuc_to_aa,s_aa);
% showalignment(alignment)

s_nuc_to_aa_correct_index = s_nuc_to_aa(index(1):end);

[score2,alignment2,index2] = swalign(s_nuc_to_aa_correct_index,s_aa);
% showalignment(alignment2)


% aa_index = 750;

% amino_acid_in_mahoney = s_aa(aa_index)

index_nuc = ((aa_index+index(1)-1)*3) -2  +1:((aa_index+index(1)-1)*3)+1;

nuc = s_nuc(index_nuc);

% amino_acid_obtained_from_nuc_indices = nt2aa(nuc)

