function plot_one_and_two_point_correlations(no_of_mutating_sites,freqs_data,freqs_sampler,color,markersize)

%Inputs: 1. no_of_mutating_sites
%        2. freqs_data -> excel file containing point mutations of the
%        data. Obtained by pasting data from p file in excel file
%        3. freqs_sampler -> excel file containing point mutations of
%        the sampler. Obtained by pasting data from p file in excel file
%Outputs: Two plots

%Data
% numdata_data = xlsread('vp1-sec1-data.xlsx');
numdata_data = xlsread(freqs_data);


one_point_freq_data = numdata_data(1:no_of_mutating_sites,:);
% two_point_freq_data = numdata_data(no_of_mutating_sites+1:end,:);

one_point_freq_all_data = one_point_freq_data(:);
one_point_freq_all_data(isnan(one_point_freq_all_data))=[];

% two_point_freq_all_data = two_point_freq_data(:);
% two_point_freq_all_data(isnan(two_point_freq_all_data))=[];

% no_two_point_freq_all_data_zeros = sum(two_point_freq_all_data==0);
% perc_sparsity_two_points_data = no_two_point_freq_all_data_zeros/length(two_point_freq_all_data)*100;

% no_two_point_freq_all_data_thresh = sum(two_point_freq_all_data<=1e-5);
% perc_sparsity_two_points_data_thresh = no_two_point_freq_all_data_thresh/length(two_point_freq_all_data)*100



%Sampler
% numdata = xlsread('vp1-sec1-sampler.xlsx');
numdata = xlsread(freqs_sampler);

one_point_freq_sampler = numdata(1:no_of_mutating_sites,:);
% two_point_freq_sampler = numdata(no_of_mutating_sites+1:end,:);

one_point_freq_all_sampler = one_point_freq_sampler(:);
one_point_freq_all_sampler(isnan(one_point_freq_all_sampler))=[];

% two_point_freq_all_sampler = two_point_freq_sampler(:);
% two_point_freq_all_sampler(isnan(two_point_freq_all_sampler))=[];

%%

    
    figure
    plot(one_point_freq_all_data,one_point_freq_all_sampler,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',markersize)
    hold on;
    max_value = max(max(one_point_freq_all_data,one_point_freq_all_sampler));
    plot(0:max_value/20:max_value,0:max_value/20:max_value,'k')
    xlabel('One-point correlations (MSA)')
    ylabel('One-point correlations ')

    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'      , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'LineWidth'   , .5        );
    set(gcf, 'PaperPositionMode', 'auto');
       

%%
%     figure
%     plot(two_point_freq_all_data,two_point_freq_all_sampler,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',markersize)
%     hold on
%     max_value = max(max(two_point_freq_all_data,two_point_freq_all_sampler));
%     plot(0:max_value/20:max_value,0:max_value/20:max_value,'k')
%     xlabel('Double mutant probability (MSA)')
%     ylabel('Double mutant probability')
%     set(gca, ...
%         'Box'         , 'off'     , ...
%         'TickDir'     , 'out'     , ...
%         'TickLength'  , [.02 .02] , ...
%         'XMinorTick'  , 'on'      , ...
%         'YMinorTick'  , 'on'      , ...
%         'XColor'      , [.1 .1 .1], ...
%         'YColor'      , [.1 .1 .1], ...
%         'LineWidth'   , .5        );
%     set(gcf, 'PaperPositionMode', 'auto');