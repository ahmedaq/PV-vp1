function peaks_antigenic_analysis(input_data_file,data_peaks_rank)

load(input_data_file)
load(data_peaks_rank)

run startup.m

%% finding the mutated sites in each peak seq

antigenic_sites_VP1 = [93 95:101 254 168 221:223 224 226];  %\Hogle1989 (only mutating)

for kk = 1:no_of_peak_seqs_vp1
    indices_mutation = find(unique_minima_seqs_sorted_vp1(kk,:)~=0);
    indx_aa_mutated = [];
    for mm = 1:length(indices_mutation)
        temp = find(phi_cumulative>=indices_mutation(mm));
        indx_aa_mutated(mm) = temp(1);
    end
    indx_aa_mutated_seq{kk} = (indx_aa_mutated);
    indx_aa_mutated_seq_true{kk} = true_indices(indx_aa_mutated);
    indx_aa_mutated_seq_true_in_antigenic_sites{kk} = ...
        indx_aa_mutated_seq_true{kk}(ismember(indx_aa_mutated_seq_true{kk},antigenic_sites_VP1));
    no_of_antigenic_sites_in_mutated_seq(kk) = ...
        sum(ismember(indx_aa_mutated_seq_true{kk},antigenic_sites_VP1));
end

%% Finding the amino acids at the mutant antigenic sites

% antigenic_sites_VP1 = [221:223 224 226];  %\Minor1990 
antigenic_sites_VP1_mut_index = rev_translation_indices(antigenic_sites_VP1,true_indices);
phi_cum = [0 phi_cumulative]; %for proper indexing in loop

% temp = [antigenic_sites_VP1_mut_index antigenic_sites_VP1_mut_index(end)+1] ;
for kk = 1:length(antigenic_sites_VP1_mut_index)
    indices_in_msa_aa_ex{kk} = phi_cum(antigenic_sites_VP1_mut_index(kk))+1:phi_cum(antigenic_sites_VP1_mut_index(kk)+1);
end

for kk = 1:no_of_peak_seqs_vp1
    for mm = 1:length(indices_in_msa_aa_ex)
        if length(indices_in_msa_aa_ex{mm})==1
            if unique_minima_seqs_sorted_vp1(kk,indices_in_msa_aa_ex{mm})==0
                seq_aa_antigenic_sites_in_peak_seq(kk,mm)=mutant_order{antigenic_sites_VP1_mut_index(mm)}(1);
            else
                seq_aa_antigenic_sites_in_peak_seq(kk,mm)=mutant_order{antigenic_sites_VP1_mut_index(mm)}(2);
            end
        else
            if sum(unique_minima_seqs_sorted_vp1(kk,indices_in_msa_aa_ex{mm})==[0 0])==2
                seq_aa_antigenic_sites_in_peak_seq(kk,mm)=mutant_order{antigenic_sites_VP1_mut_index(mm)}(1);
            elseif sum(unique_minima_seqs_sorted_vp1(kk,indices_in_msa_aa_ex{mm})==[0 1])==2
                seq_aa_antigenic_sites_in_peak_seq(kk,mm)=mutant_order{antigenic_sites_VP1_mut_index(mm)}(2);
            else
                seq_aa_antigenic_sites_in_peak_seq(kk,mm)=mutant_order{antigenic_sites_VP1_mut_index(mm)}(3);
            end
        end
    end
end

unique_aa_combo_at_antigenic_sites = unique(seq_aa_antigenic_sites_in_peak_seq,'rows');
no_of_unique_aa_combo_at_antigenic_sites = size(unique_aa_combo_at_antigenic_sites,1);

%% 

list_indx_aa_mutated_seq_true = [];
for kk =1:no_of_peak_seqs_vp1
    list_indx_aa_mutated_seq_true = [list_indx_aa_mutated_seq_true  indx_aa_mutated_seq_true{kk}];
end

list_indx_aa_mutated_seq_true = unique(list_indx_aa_mutated_seq_true);

per_peak_seqs_with_mutations_at_antigenic_sites = ...
    length(find(no_of_antigenic_sites_in_mutated_seq~=0))/no_of_peak_seqs_vp1*100;

% finding enrichment and p values

for kk = 1:no_of_peak_seqs_vp1
    frac_antigenic_sites_in_peak_seq(kk) = length(indx_aa_mutated_seq_true_in_antigenic_sites{kk})/...
        length(indx_aa_mutated_seq_true{kk});
    pvalue(kk) = pvalue_counting(302,length(antigenic_sites_VP1),length(indx_aa_mutated_seq_true{kk}),...
        length(indx_aa_mutated_seq_true_in_antigenic_sites{kk}));
end


%%
pvalue_vec = [];
figure;
plot(frac_antigenic_sites_in_peak_seq,-log10(pvalue),'o','MarkerFaceColor',blue,'MarkerEdgeColor','w','MarkerSize',8)
hold on

a = 0:0.1:(max(frac_antigenic_sites_in_peak_seq))+0.1;
plot(a,-log10(0.05)*ones(1,length(a)),'k--','LineWidth',1)
% xlabel('Enrichment in antigenic mutations')
xlabel('Enrichment')
xlim([0 0.8])
ylim([0 3.5])
% ylabel('$$-\log_{10} (p)$$','interpreter','latex')
ylabel('-log_{10}(p-value)','interpreter','tex')
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...    
    'YTick'       , [-0.5:0.5:3.5],...
    'XTick'       , 0:0.2:0.8,...
    'LineWidth'   , 1        );

%%

figure; 
barh(1:no_of_peak_seqs_vp1,no_of_antigenic_sites_in_mutated_seq/length(antigenic_sites_VP1)*100,0.5,'FaceColor',blue,'EdgeColor',blue)
ylabel('Peak rank')
xlabel('%Mutation on antigenic sites w.r.t. consensus sequence','interpreter','tex')
xlim([0 80])
ylim([0 27])
hold on

for kk = 1:no_of_peak_seqs_vp1
    x_axis_coord = 31;
    [~,alignment] = nwalign(seq_aa_antigenic_sites_in_peak_seq(1,:),seq_aa_antigenic_sites_in_peak_seq(kk,:));
    indices_mut = find(alignment(2,:) ~= '|');
    for mm = 1:size(seq_aa_antigenic_sites_in_peak_seq,2)
        if ismember(mm,indices_mut)==1
%             gray = [0.6 0.6 0.6];
            text(x_axis_coord,kk+0.2,seq_aa_antigenic_sites_in_peak_seq(kk,mm),'FontSize',7,'FontWeight','bold')
%             text(x_axis_coord,kk+0.2,seq_aa_antigenic_sites_in_peak_seq(kk,mm),'FontSize',7,'FontWeight','bold','Color',white,'BackgroundColor',[0.6 0.6 0.6])
        else
            text(x_axis_coord,kk+0.2,seq_aa_antigenic_sites_in_peak_seq(kk,mm),'FontSize',7,'color',blue)
        end
        x_axis_coord = x_axis_coord + 3;
    end
end

x_axis_coord = 31;
for mm = 1:size(seq_aa_antigenic_sites_in_peak_seq,2)
    text(x_axis_coord,26+0.2,int2str(antigenic_sites_VP1(mm)),'FontSize',6)
    x_axis_coord = x_axis_coord + 3;
end

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'on'      , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...    
    'YTick'       , 0:5:30,...
    'LineWidth'   , 1        );


% %% Finding the antigenic mutations in each sequence falling in the valley
% for kk = 1:no_of_peak_seqs_vp1
%     seqs_in_peak = indices_seqs_that_converged_to_peak_seq_vp1{kk};
%     for nn = 1:length(seqs_in_peak)
%         indices_mutation = find(msa_aa_ex(seqs_in_peak(nn),:)~=0);
%         indx_aa_mutated = [];
%         for mm = 1:length(indices_mutation)
%             temp = find(phi_cumulative>=indices_mutation(mm));
%             indx_aa_mutated(mm) = temp(1);
%         end
%         indx_aa_mutated_seq2{kk}{nn} = (indx_aa_mutated);
%         indx_aa_mutated_seq_true2{kk}{nn} = true_indices(indx_aa_mutated);
%         indx_aa_mutated_seq_true_in_antigenic_sites2{kk}{nn} = ...
%             indx_aa_mutated_seq_true2{kk}{nn}(ismember(indx_aa_mutated_seq_true2{kk}{nn},antigenic_sites_VP1));
%         %Calculate Hamming distance of each sequence from the respective
%         %valley sequence in the antigenic region only...
%         antigenic_HD_of_seq_from_its_valleySeq{kk}(nn) = ...
%             length(union(indx_aa_mutated_seq_true_in_antigenic_sites{kk},indx_aa_mutated_seq_true_in_antigenic_sites2{kk}{nn}))-...
%         length(intersect(indx_aa_mutated_seq_true_in_antigenic_sites{kk},indx_aa_mutated_seq_true_in_antigenic_sites2{kk}{nn}));
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Boxplot figure
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% antigenic_HD = [antigenic_HD_of_seq_from_its_valleySeq{1} antigenic_HD_of_seq_from_its_valleySeq{2} ...
%     antigenic_HD_of_seq_from_its_valleySeq{3} antigenic_HD_of_seq_from_its_valleySeq{4} ...
%     antigenic_HD_of_seq_from_its_valleySeq{5}];
% 
% grp = [zeros(1,length(antigenic_HD_of_seq_from_its_valleySeq{1})),ones(1,length(antigenic_HD_of_seq_from_its_valleySeq{2}))...
%     2*ones(1,length(antigenic_HD_of_seq_from_its_valleySeq{3})),3*ones(1,length(antigenic_HD_of_seq_from_its_valleySeq{4}))...
%     4*ones(1,length(antigenic_HD_of_seq_from_its_valleySeq{5})),...
%     ];
% 
% 
% jitter_value = 0.4;
% whisker_value = 1.5;
% outlier_size = 8;
% plot_style = 'traditional';%'compact
% % colors_vec = 'bbbbbb'; %'bmrkyg';
% widths_box_value = 0.2;
% 
% 
% color_scheme2(1,:) = [0 0 0];
% color_scheme2(2:5,:) = color_scheme(1:4,:);
% color_scheme2(6,:) = [150 150 150]./255; %darkgray;
% color_scheme2(7:10,:) = color_scheme(5:8,:);
% 
% 
% % figure;
% % subplot(211)
% figure
% bh = boxplot(antigenic_HD,grp,'Orientation','vertical','PlotStyle',plot_style,...   
%     'labels',1:5,...
%     'whisker',whisker_value, 'jitter', jitter_value,'symbol','.',...
%     'color',color_scheme2,'outliersize',outlier_size,...
%     'widths',widths_box_value);
% set(bh,'linewidth',.5);
% % h=findobj(gca,'tag','Outliers');delete(h);%deleting outlier
% title('Antigenic makeup comparison between sequences associated with a peak and the respective peak sequence','FontWeight', 'Normal','FontSize',11)
% ylabel('Hamming distance')
% xlabel('Peak')
% 
% h=findobj(gca,'tag','Median');
% for kk=1:5
%     h(kk).LineWidth = 3;
% end
% 
% set(gca, ...
%   'Box'         , 'on'     , ...
%   'TickDir'     , 'in'     , ...
%   'TickLength'  , [.01 .01] , ...
%   'XMinorTick'  , 'off'      , ...
%   'YMinorTick'  , 'off'      , ...
%   'YGrid'       , 'off'      , ...
%   'XColor'      , [0 0 0], ...
%   'YColor'      , 'k', ...  
%   'LineWidth'   , .5        );
% 
% 
% %% Other valleys
% 
% %HD of each sequence in valley 1 from other valley sequences in the antigenic region only...
% for other_valley = 2:5
%     for nn = 1:length(indices_seqs_that_converged_to_peak_seq_vp1{1})
%         antigenic_HD_of_seqs_in_valley1_from_valleySeq{other_valley}(nn) = ...
%             length(union(indx_aa_mutated_seq_true_in_antigenic_sites{other_valley},indx_aa_mutated_seq_true_in_antigenic_sites2{1}{nn}))-...
%             length(intersect(indx_aa_mutated_seq_true_in_antigenic_sites{other_valley},indx_aa_mutated_seq_true_in_antigenic_sites2{1}{nn}));
%     end
% end
% 
% %Boxplot figure
% 
% antigenic_HD_other = [antigenic_HD_of_seqs_in_valley1_from_valleySeq{2} antigenic_HD_of_seqs_in_valley1_from_valleySeq{3} ...
%     antigenic_HD_of_seqs_in_valley1_from_valleySeq{4} antigenic_HD_of_seqs_in_valley1_from_valleySeq{5}];
% 
% grp = [zeros(1,length(antigenic_HD_of_seqs_in_valley1_from_valleySeq{2})),ones(1,length(antigenic_HD_of_seqs_in_valley1_from_valleySeq{3}))...
%     2*ones(1,length(antigenic_HD_of_seqs_in_valley1_from_valleySeq{4})),3*ones(1,length(antigenic_HD_of_seqs_in_valley1_from_valleySeq{5}))];
% 
% jitter_value = 0.4;
% whisker_value = 1.5;
% outlier_size = 8;
% plot_style = 'traditional';%'compact
% % colors_vec = 'bbbbbb'; %'bmrkyg';
% widths_box_value = 0.2;
% 
% 
% color_scheme2(1,:) = [0 0 0];
% color_scheme2(2:5,:) = color_scheme(1:4,:);
% color_scheme2(6,:) = [150 150 150]./255; %darkgray;
% color_scheme2(7:10,:) = color_scheme(5:8,:);
% 
% 
% % figure;
% % subplot(211)
% figure
% bh = boxplot(antigenic_HD_other,grp,'Orientation','vertical','PlotStyle',plot_style,...   
%     'labels',2:5,...
%     'whisker',whisker_value, 'jitter', jitter_value,'symbol','.',...
%     'color',color_scheme2(2:end,:),'outliersize',outlier_size,...
%     'widths',widths_box_value);
% set(bh,'linewidth',.5);
% % h=findobj(gca,'tag','Outliers');delete(h);%deleting outlier
% title('Antigenic makeup comparison between sequences associated with peak 1 and other peak sequences','FontWeight', 'Normal','FontSize',11)
% ylabel('Hamming distance')
% xlabel('Peak')
% 
% h=findobj(gca,'tag','Median');
% for kk=1:4
%     h(kk).LineWidth = 3;
% end
% 
% set(gca, ...
%   'Box'         , 'on'     , ...
%   'TickDir'     , 'in'     , ...
%   'TickLength'  , [.01 .01] , ...
%   'XMinorTick'  , 'off'      , ...
%   'YMinorTick'  , 'off'      , ...
%   'YGrid'       , 'off'      , ...
%   'XColor'      , [0 0 0], ...
%   'YColor'      , 'k', ...  
%   'LineWidth'   , .5        );
% 
% 
% 
