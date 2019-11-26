function plot_rank_peaks(data_peaks,H,input_data_file)

load(data_peaks)
load(input_data_file)

[ns_ex,ls_ex] = size(msa_aa_ex); %extended 

no_of_times_a_better_seq_available_vp1 = no_of_times_a_better_seq_available;
energy_min_seq_vp1 = energy_min_seq;
minima_seq_vp1 = minima_seq;
H_vp1 = H;

[unique_minima_seqs_vp1,b_vp1,c_vp1] = unique(minima_seq_vp1,'rows');

for kk = 1:length(b_vp1)
    no_seqs_copies_vp1(kk) = sum(c_vp1==kk);
end

[sort_no_seqs_copies_vp1,indx_sort_seqs_vp1] = sort(no_seqs_copies_vp1,'descend');

unique_minima_seqs_sorted_vp1 = unique_minima_seqs_vp1(indx_sort_seqs_vp1,:);

for kk = 1:size(unique_minima_seqs_sorted_vp1,1)
    energy_min_seq_sorted_vp1(kk) = calcSeqEnergy(unique_minima_seqs_sorted_vp1(kk,:),H_vp1);
end

for kk = 1:size(unique_minima_seqs_sorted_vp1,1)
    indices_seqs_that_converged_to_peak_seq_vp1{kk} = ...
        find(round2dp(energy_min_seq_vp1,5)==round2dp(energy_min_seq_sorted_vp1(kk),5));
end

no_of_peak_seqs_vp1 = length(sort_no_seqs_copies_vp1)

%% saving data

save vp1_peaks_rank

%% Plot

figure;
loglog(1:no_of_peak_seqs_vp1,sort_no_seqs_copies_vp1/length(c_vp1),'o',...
    'MarkerFaceColor',blue,'MarkerEdgeColor','w','MarkerSize',8)
% hold on
% loglog(1:no_of_peak_seqs_p24,sort_no_seqs_copies_p24/length(c_p24),'o',...
%     'MarkerFaceColor',red,'MarkerEdgeColor',red,'MarkerSize',6)
% loglog(1:no_of_peak_seqs_nef,sort_no_seqs_copies_nef/length(c_nef),'o',...
%     'MarkerFaceColor',green,'MarkerEdgeColor',green,'MarkerSize',6)
% loglog(1:no_of_peak_seqs_vp1,sort_no_seqs_copies_vp1/length(c_vp1),'o',...
%     'MarkerFaceColor',blue,'MarkerEdgeColor',blue,'MarkerSize',6)
% legend('vp1','p24','nef')
% legend boxoff
xlabel('Peak rank')
ylabel('Fraction of sequences in peak')
xlim([0.9 1e2])
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'XTick'       , [1 1e1 1e2 1e3], ...
    'YTick'       , [1e-3 1e-2 1e-1 1], ...
    'LineWidth'   , 1        );