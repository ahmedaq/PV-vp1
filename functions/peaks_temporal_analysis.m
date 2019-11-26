function year_seqs_peak_seq = peaks_temporal_analysis(input_data_file,data_peaks_rank,data_country_vp1,data_year_vp1)

load(data_peaks_rank)
load(data_country_vp1)
load(data_year_vp1)
load(input_data_file)

% run startup.m

%%
for kk = 1:no_of_peak_seqs_vp1
    year_seqs_peak_seq{kk} = year_vp1(indices_seqs_that_converged_to_peak_seq_vp1{kk});
end

%Adding years of valley 10 manually based on the associated report
year_seqs_peak_seq{10}(1:17) = 1992;
year_seqs_peak_seq{10}(18:35) = 1993;


%% Plot

color_scheme2(1,:) = [0 0 0];
color_scheme2(2:5,:) = color_scheme(1:4,:);
color_scheme2(6,:) = [150 150 150]./255; %darkgray;
color_scheme2(7:10,:) = color_scheme(5:8,:);

figure('units','normalized','outerposition',[0 0 1 1])
for kk = 1:10
    subplot(5,2,kk)
    title(sprintf('Peak %d',kk))
    indices_year_info_available = find(year_seqs_peak_seq{kk}~=0);
    if ~isempty(indices_year_info_available)        
        unique_years = unique(year_seqs_peak_seq{kk}(indices_year_info_available));
        x = histc(year_seqs_peak_seq{kk}(indices_year_info_available),min(unique_years):max(unique_years));
        perc_x = x/sum(x)*100;
        perc_x(perc_x>0 & perc_x<=5) = 5;
        h = area(min(unique_years):max(unique_years),(perc_x),'LineWidth',0.5);
        h(1).FaceColor = color_scheme2(kk,:);
        alpha 0.8
        xlim([1960 2015])
        ylim([0 70])
        if kk >= 9
            xlabel('Year','FontSize',10)
        end
        if kk == 5
            ylabel('Percentage','FontSize',10)
        end  
    end
    post_proc_fig
    set(gca,...
        'box','off',...
        'TickLength',[.02 .02],...
        'XTick',1960:10:2010,...
        'YTick',[0 70],...
        'LineWidth'   , .5        );
    
end