
function mean_hamming_distance_max_epsilon = compute_neutrality_landscape(data_protein)

load(data_protein)

[ns_ex,ls_ex] = size(msa_aa_ex);
[ns,ls] = size(msa);

%% Generate random sequences
nrs = 1e5;
M = 5;

random_start_seqs = generate_random_start_seqs(nrs,M,...
    msa,msa_aa_ex,phi_curr,phi_cumulative,mutant_order);

%% Neutrality

L = 1000;  %number of steps ;500
epsilon_values = [1e-01 3e-01 5e-01 8e-01 3]; %threshold

phi_cum = [0 phi_cumulative]; %for proper indexing in loop

for ep = 1:length(epsilon_values)
    
    epsilon = epsilon_values(ep);
    
    parfor kk = 1:nrs
        
        init_seq = random_start_seqs(kk,:);
        initial_seq = random_start_seqs(kk,:);
        hamming_distance = [];
        
        for ss = 1:L
            
            energy_initial_seq = calcSeqEnergy(initial_seq,H);
            
            mm = randi(ls);%select a site to mutate
            
            if phi_curr(mm) == 1
                new_seq = initial_seq;
                new_seq(phi_cum(mm)+1:phi_cum(mm+1)) = ~initial_seq(phi_cum(mm)+1:phi_cum(mm+1));
                energy_seq = calcSeqEnergy(new_seq,H);
                
            elseif phi_curr(mm) ~= 1
                new_seq = initial_seq;
                possibilities = [zeros(1,phi_curr(mm)) ; eye(phi_curr(mm))];
                current_pos = find(sum(~xor(repmat(new_seq(1,phi_cum(mm)+1:phi_cum(mm+1)),phi_curr(mm)+1,1),...
                    possibilities),2)==phi_curr(mm));
                remaining_pos = setdiff(1:phi_curr(mm)+1,current_pos);
                sel_aa = randi(phi_curr(mm));
                new_seq(phi_cum(mm)+1:phi_cum(mm+1)) = possibilities(remaining_pos(sel_aa),:);
                energy_seq = calcSeqEnergy(new_seq,H);
                
            end
            
            if abs((energy_seq) - (energy_initial_seq)) <= epsilon
                hamming_distance = [hamming_distance; sum((new_seq-init_seq)~=0)];
                initial_seq = new_seq;
            else
                hamming_distance = 0;
            end
        end
        hamming_distance_max(kk) = max(hamming_distance);
        
    end
    hamming_distance_max_epsilon{ep} = hamming_distance_max;
    mean_hamming_distance_max_epsilon(ep) = mean(hamming_distance_max);
    %     ep
end
