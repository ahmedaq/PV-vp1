function [indices_conserved_sites,no_of_conserved_sites] = find_conserved_sites(msa)

[ns,ls]=size(msa); %ns=no.of seq, ls=length of seq
profile_seq = seqprofile(msa);
[freq_wt,pos_wt] = max(profile_seq);
Cseq = int2aa(pos_wt);
% Cseq = seqconsensus(msa);
Cseq_mtrx = repmat(Cseq,ns,1);
msa_binary = double(msa==Cseq_mtrx); %replacing mutation by 0 and no mutation by 1
freq_single_mutation = mean(msa_binary);
site_freq_no_mutation = find(freq_single_mutation==1);

no_of_mutations_per_site = sum(double(~msa_binary),1);
site_with_all_blank_mutations = [];
for kk = 1:ls
    no_of_Bs(kk) = length(find(msa(:,kk)=='B'));
    no_of_Js(kk) = length(find(msa(:,kk)=='J'));
    no_of_Os(kk) = length(find(msa(:,kk)=='O'));
    no_of_Us(kk) = length(find(msa(:,kk)=='U'));
    no_of_Xs(kk) = length(find(msa(:,kk)=='X'));
    no_of_Zs(kk) = length(find(msa(:,kk)=='Z'));
    if ((no_of_Bs(kk) ~= 0) && (no_of_Bs(kk) == no_of_mutations_per_site(kk))) ...
       || ((no_of_Js(kk) ~= 0) && (no_of_Js(kk) == no_of_mutations_per_site(kk))) ...     
       || ((no_of_Os(kk) ~= 0) && (no_of_Os(kk) == no_of_mutations_per_site(kk))) ...
       || ((no_of_Us(kk) ~= 0) && (no_of_Us(kk) == no_of_mutations_per_site(kk))) ...
       || ((no_of_Xs(kk) ~= 0) && (no_of_Xs(kk) == no_of_mutations_per_site(kk))) ...
       || ((no_of_Zs(kk) ~= 0) && (no_of_Zs(kk) == no_of_mutations_per_site(kk)))
       
            site_with_all_blank_mutations = [site_with_all_blank_mutations kk]; 
    end
end

% site_with_all_blank_mutations

site_with_blanks_mutations_gt_12pc = [];
for kk = 1:ls
    per_no_of_Bs(kk) = length(find(msa(:,kk)=='B'))/ns*100;
    if per_no_of_Bs(kk) > 12.5
       site_with_blanks_mutations_gt_12pc = [site_with_blanks_mutations_gt_12pc kk];
    end
end

% site_with_blanks_mutations_gt_12pc
%%%
% if strcmp(remove_sites_with_all_blank_mutations,'yes')
indices_conserved_sites = [site_freq_no_mutation site_with_all_blank_mutations site_with_blanks_mutations_gt_12pc];
% end
% msa(:,[site_freq_no_mutation])=[];
% msa_binary(:,[site_freq_no_mutation])=[]; 
% freq_single_mutation_after_processing = mean(msa_binary);

% true_indices = setdiff([1:ls],site_freq_no_mutation);
% percentage_conserved_sites = length(site_freq_no_mutation)/ls*100

no_of_conserved_sites = length(indices_conserved_sites);