function [energy_matrix, ak, ak_log_energy] = correlation_calculation(nrs,M,no_of_steps,...
    ls_msa,msa_aa_ex,phi_curr,phi_cumulative,mutant_order,H)


[ns_ex,ls_ex] = size(msa_aa_ex);
%[ns,ls] = size(msa);

msa_ex_unique = unique(msa_aa_ex,'rows');
[ns_unique,~] = size(msa_ex_unique);


%% Generate random sequences

random_start_seqs = generate_random_start_seqs(nrs,M,...
    ls_msa,msa_aa_ex,phi_curr,phi_cumulative,mutant_order);


%% correlation length calculation

phi_cum = [0 phi_cumulative]; %for proper indexing in loop
energy_matrix = zeros(nrs,no_of_steps+1);
%first column is the fitness of randomly selected sequences

%calculating the energy of the random start seqs
for kk = 1:nrs
    energy_matrix(kk,1) = calcSeqEnergy(random_start_seqs(kk,:),H);
end

tic
for k = 2:no_of_steps+1
    for kk = 1:nrs
        %select a site to mutate
        mm = randi(ls_msa);
%         new_seq = random_start_seqs(kk,:); %not saving seqs
        
        if phi_curr(mm) == 1
            random_start_seqs(kk,phi_cum(mm)+1:phi_cum(mm+1)) = ~random_start_seqs(kk,phi_cum(mm)+1:phi_cum(mm+1));
            energy_matrix(kk,k) = calcSeqEnergy(random_start_seqs(kk,:),H);
        
        elseif phi_curr(mm) ~= 1
            possibilities = [zeros(1,phi_curr(mm)) ; eye(phi_curr(mm))];
            current_pos = find(sum(~xor(repmat(random_start_seqs(1,phi_cum(mm)+1:phi_cum(mm+1)),phi_curr(mm)+1,1),...
                    possibilities),2)==phi_curr(mm));
            remaining_pos = setdiff(1:phi_curr(mm)+1,current_pos);
            sel_aa = randi(phi_curr(mm));
            random_start_seqs(kk,phi_cum(mm)+1:phi_cum(mm+1)) = possibilities(remaining_pos(sel_aa),:);
            energy_matrix(kk,k) = calcSeqEnergy(random_start_seqs(kk,:),H);
            
        end
    end
end
toc

%%%%%%%%%%%%%%%%%%%%
%changing energy of a seq = 0 to a very small value to avoid problem in
%taking log
%%%%%%%%%%%%%%%%%%%%
energy_matrix(find(energy_matrix==0))=1e-10;

%finding correlation coefficient for each step
ak = zeros(no_of_steps,1);
ak_log_energy = zeros(no_of_steps,1);
for k = 1:no_of_steps
    ak_log_energy(k) = corr(log10(energy_matrix(:,1)),log10(energy_matrix(:,k+1)));
    ak(k) = corr(energy_matrix(:,1),energy_matrix(:,k+1));
end
