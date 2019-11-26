%% Determining number of peaks using a steepest-ascent method

function determine_local_peaks(input_data_file,H)

load(input_data_file)

[ns_ex,ls_ex] = size(msa_aa_ex);
[ns,ls] = size(msa); 

phi_cum = [0 phi_cumulative]; %for proper indexing in loop
minima_seq = [];

tic
parfor kk = 1:ns_ex
    
    initial_seq = msa_aa_ex(kk,:); %select each sequence in the sma
    local_minima = 0;
    pp = 0;
    ppp = 0;
    while local_minima == 0
               
        %loop to flip each sequence and obtain the corresponding energy
        %value
        nn = 1;
        new_seq_store = cell(1,ls_ex);
        f_new = zeros(1,ls_ex);
        for mm = 1:ls
            
            new_seq = initial_seq;
            
            f_curr = calcSeqEnergy(new_seq,H); %energy of current sequence
            
            if phi_curr(mm) == 1            
%                 new_seq(1,phi_cum(mm+1)) = ~initial_seq(1,phi_cum(mm+1));
                new_seq(1,nn) = ~initial_seq(1,nn);
                f_new(nn) = calcSeqEnergy(new_seq,H);
                new_seq_store{nn} = new_seq;
                nn = nn+1;
            elseif phi_curr(mm) ~= 1%== 2
                possibilities = [zeros(1,phi_curr(mm)) ; eye(phi_curr(mm))];
                current_pos = find(sum(~xor(repmat(new_seq(1,phi_cum(mm)+1:phi_cum(mm+1)),phi_curr(mm)+1,1),...
                    possibilities),2)==phi_curr(mm));
                remaining_pos = setdiff(1:phi_curr(mm)+1,current_pos);
                for rr = 1:phi_curr(mm)
                    new_seq1 = new_seq;
                    new_seq1(1,phi_cum(mm)+1:phi_cum(mm+1)) = possibilities(remaining_pos(rr),:);
                    f_new(nn) = calcSeqEnergy(new_seq1,H);
                    new_seq_store{nn} = new_seq1;
                    nn = nn+1;
                end
            end
        end
        
        [min_f_new,indx_f_new] = min(f_new); %max value of energy in mutated sequences
	       
        if min_f_new < f_curr %if a mutated seq has energy lower than prev
            initial_seq = new_seq_store{indx_f_new};
            pp = pp+1;  
            no_of_times_a_better_seq_available(kk)=pp;
	    
	    if sum(ismember(f_new,min_f_new))>1 %multiple options available?
		ppp = ppp+1;
		multiple_minima(kk) = ppp;
       	    end
        else %if the prev seq is the local minima
            minima_seq(kk,:) = initial_seq; %store peak sequence
            energy_min_seq(kk) = f_curr; %store energy value
            local_minima = 1; %this will terminate while loop
        end        
        
    end

end
toc

avg_no_steps_at_each_seq_to_reach_local_minima = mean(no_of_times_a_better_seq_available);
std_no_steps_at_each_seq_to_reach_local_minima = std(no_of_times_a_better_seq_available);

save vp1_peaks minima_seq energy_min_seq no_of_times_a_better_seq_available multiple_minima
