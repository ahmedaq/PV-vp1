function random_start_seqs = generate_random_start_seqs(nrs,M,...
    msa,msa_aa_ex,phi_curr,phi_cumulative,mutant_order)

[ns_ex,ls_ex] = size(msa_aa_ex);
[ns,ls] = size(msa);

msa_ex_unique = unique(msa_aa_ex,'rows');
[ns_unique,~] = size(msa_ex_unique);

% nrs = 1e5; %no of random seqs
% we want to generate nrs seqs in the M hamming distance neighborhood of msa

% M = 5; %max number of mutations in a randomly selected MSA seq

phi_cum = [0 phi_cumulative]; %for proper indexing in loop

% matlabpool open
parfor kk = 1:nrs
    %select a random sequence from observed msa
    sel_seq_from_msa = msa_ex_unique(randi(ns_unique),:);
    %randomly find the number of mutations to make in it
    no_mutations = randi(M); %mutations can vary from 1 to M
    %randomly select the sites at which mutations should happen
    sites_mutations = randi(ls,1,no_mutations);
    
    for mm = 1:no_mutations
        if phi_curr(sites_mutations(mm))==1
            sel_seq_from_msa(phi_cum(sites_mutations(mm))+1:phi_cum(sites_mutations(mm)+1)) ...
                = ~(sel_seq_from_msa(phi_cum(sites_mutations(mm))+1:phi_cum(sites_mutations(mm)+1)));
        elseif phi_curr(sites_mutations(mm))~=1
            possibilities = [zeros(1,phi_curr(sites_mutations(mm))) ; eye(phi_curr(sites_mutations(mm)))];
            sel_aa = randi(phi_curr(sites_mutations(mm))+1); %randomly select an amino acid at this position
            sel_seq_from_msa(phi_cum(sites_mutations(mm))+1:phi_cum(sites_mutations(mm)+1))...
                = possibilities(sel_aa,:);       
        end
        
%         elseif phi_curr(sites_mutations(mm))~=2
%             possibilities = [zeros(1,phi_curr(mm)) ; eye(phi_curr(mm))];
%             sel_aa = randi(3); %randomly select an amino acid at this position
%             sel_seq_from_msa(phi_cum(sites_mutations(mm))+1:phi_cum(sites_mutations(mm)+1))...
%                 = possibilities(sel_aa,:);
%         elseif phi_curr(sites_mutations(mm))==3
%             possibilities = [0 0 0; 0 0 1; 0 1 0; 1 0 0];
%             sel_aa = randi(4); %randomly select an amino acid at this position
%             sel_seq_from_msa(phi_cum(sites_mutations(mm))+1:phi_cum(sites_mutations(mm)+1))...
%                 = possibilities(sel_aa,:);
%         elseif phi_curr(sites_mutations(mm))==4
%             possibilities = [0 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0; 1 0 0 0];
%             sel_aa = randi(5); %randomly select an amino acid at this position
%             sel_seq_from_msa(phi_cum(sites_mutations(mm))+1:phi_cum(sites_mutations(mm)+1))...
%                 = possibilities(sel_aa,:);
%         else 
%             sprintf('Error! Code not written for phi_curr (number of mutants at a site) > 4')
%             break
%         end
        
    end
    random_start_seqs(kk,:) = sel_seq_from_msa;    
end
% matlabpool close

% save random_start_seqs random_start_seqs




