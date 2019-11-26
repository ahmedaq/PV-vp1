
function perc = ZeroTempMCSim(input_data_file,data_peaks_rank)

load(input_data_file)
load(data_peaks_rank)

[Nseq,total_length] = size(msa_aa_ex);
[Nseq,protein_length_aa] = size(msa);

% protein_length_aa = N;
% total_length = N;
strMethod = ' ';
strLargeOrSmall = 'large'; % 'small';

thin = 1;1e1;%1
burnin = 0;1e5;%0
nosim = 5e5;%1e7 %5*1e7; % 1e8;

MCsweepLength = 1e4; % 2*1e4; % 5*1e4; % 1e5;
numParallel = 1; %4;%number of parallel chains
betaArray = inf; %[1 0.9 0.8 0.7]; % [1 0.85 0.7 0.55]; %Temperature

J_MINFLOW_mat = convertJohnFLParamsToRayParamsFormat(H);

% seedSeq = double(rand(1,protein_length_aa) > 0.1);
% seedSeq = zeros(1,protein_length_aa);

param_verifyMCMC_TC_PT = cell(1,11);
param_verifyMCMC_TC_PT{1} = J_MINFLOW_mat;
param_verifyMCMC_TC_PT{2} = total_length;
% param_verifyMCMC_TC_PT{3} = msa_aa_ex; % comment out this line
param_verifyMCMC_TC_PT{4} = phi_curr;
param_verifyMCMC_TC_PT{5} = phi_cumulative;
param_verifyMCMC_TC_PT{6} = protein_length_aa;
param_verifyMCMC_TC_PT{7} = strMethod;
param_verifyMCMC_TC_PT{8} = strLargeOrSmall;
param_verifyMCMC_TC_PT{9} = [thin burnin nosim];
param_verifyMCMC_TC_PT{10} = [MCsweepLength numParallel betaArray];


N_valleys = length(indices_seqs_that_converged_to_peak_seq); %number of valleys

ITER = 600; %number of MCMC runs starting from a sequence in a valley

for valley = 1:N_valleys
    
    Nseq_valley = length(indices_seqs_that_converged_to_peak_seq{valley})
    fraction_paths_converging_to_valley_seq = zeros(1,Nseq_valley);
    
    
    for kk = 1:Nseq_valley %kk is the index of a sequence in a valley
        
        %Defining the seed sequence as one sequence in the valley
        seedSeq = msa_aa_ex(indices_seqs_that_converged_to_peak_seq{valley}(kk),:);
        
        param_verifyMCMC_TC_PT{11} = seedSeq;
        
        last_seqs_all_similar_to_valley_seq = zeros(1,ITER);
        
        parfor iter = 1:ITER
        	%MCMC sampler
        	[~, samples_MCMC] = verifyMCMC_TC_PT_noMSA(param_verifyMCMC_TC_PT);
        
	        %checking the last 100*protein_length_aa samples
        
        	samples_MCMC_last_100 = samples_MCMC(end-100*protein_length_aa:end,:);
	        a{iter} = unique(samples_MCMC_last_100,'rows');
        
        	if size(a{iter},1) ~= 1
	            sprintf('Last 100xN samples are not same')
        	    no_of_unique_samples = size(a{iter},1)
		%            break
	        else
        	    %sprintf('Last 100xN samples are same')
	            if mean(a{iter}==unique_minima_seqs_sorted(valley,:))==1
        	    	last_seqs_all_similar_to_valley_seq(iter) = 1;
                    end
	        end
        end
        
        fraction_paths_converging_to_valley_seq(kk) = mean(last_seqs_all_similar_to_valley_seq);
        
        %Display
        %fraction_paths_converging_to_valley_seq(1:kk)
        
        if fraction_paths_converging_to_valley_seq(kk) ~= 1
		seqs_not_converged = find(last_seqs_all_similar_to_valley_seq==0);
		aa = [];
			for abc = 1:length(seqs_not_converged)
				aa = [aa;a{seqs_not_converged(abc)}];
			end
		other_valleys{valley}{kk} = aa;
        end
        
    end
    
    mean_over_all_valley_seqs_converging_to_valley_seq(valley) = ...
        mean(fraction_paths_converging_to_valley_seq)    
    
end

%%
no_studied_valleys = 10;

%% quantifying how many sequences went to which valleys

hist_valleys = zeros(no_studied_valleys,N_valleys+1);

for valley = 1:no_studied_valleys
    
    Nseq_valley = length(other_valleys{valley});
    
    for kk = 1:Nseq_valley %a sequence in specific valley
        
        other_valleys_that_this_seq_went_to = other_valleys{valley}{kk};
        no_other_valleys = size(other_valleys_that_this_seq_went_to,1);
        other_valley_that_this_seq_went_to = zeros(1,no_other_valleys);

        if isempty(other_valleys_that_this_seq_went_to)==0 %if the sequence goes to some other valley
            
            for mm = 1:no_other_valleys
                
                %other valley that this seq went to
                temp = find(mean(repmat(other_valleys_that_this_seq_went_to(mm,:),N_valleys,1)==unique_minima_seqs_sorted,2)==1);
                if isempty(temp)==1 %valley is not in steepest descent found valleys
                    other_valley_that_this_seq_went_to(mm) = N_valleys+1; %collect all such valleys in a new valley number N_valleys+1
                else                     
                    other_valley_that_this_seq_went_to(mm) = temp; %else save the valley number to which it went
                end
                                                
            end
            
        end
        
        %keeping track of the following vector (size 1 x N_valleys) for each sequence
        %whose i-th element indicates the number of times this
        %sequence fell in i-th valley
        hist_valleys(valley,:) = hist_valleys(valley,:) + histcounts(other_valley_that_this_seq_went_to,0.5:1:N_valleys+1+0.5); 
        
    end
    
end

clear perc
figure;
for kk = 1:no_studied_valleys
    subplot(2,5,kk)
    perc(kk,:) = (hist_valleys(kk,1:no_studied_valleys)./(ITER*length(indices_seqs_that_converged_to_peak_seq{kk})));
    perc(kk,:) = perc(kk,:)/sum(perc(kk,:))*100;
    bar(1:no_studied_valleys,perc(kk,:),0.3)
%     if kk == 13
    xlabel('Valley')
%     elseif kk == 6
    ylabel('Percentage')
%     end
    title(sprintf('Valley %d',kk))
    xlim([0 no_studied_valleys+1])
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'off'      , ...
        'YMinorTick'  , 'on'      , ...
        'YGrid'       , 'off'      , ...
        'XColor'      , [.2 .2 .2], ...
        'YColor'      , [.2 .2 .2], ...        
        'LineWidth'   , 1        );
end
