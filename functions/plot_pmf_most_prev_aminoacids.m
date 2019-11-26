function plot_pmf_most_prev_aminoacids(input_configurations_sampler,input_data_file)

A = xlsread(input_configurations_sampler);

load(input_data_file)

[ns,ls] = size(msa);

%making the modified msa according to mutant_order
msa_new = construct_msa_aa_after_entropy_compression(msa,mutant_order);


%% calculating pmf of first/sec/third dominant mutants...

% close all

dom_mutant_index_data = 2; %first dom mutant = 2; second dom mutant = 3 and so on...
msa_bin_data_dom_mutant = zeros(ns,ls);
msa_bin_sample_dom_mutant = zeros(size(A));
% indx_dom_mutant = zeros(ns,ls);
for kk = 1:ls
    if length(mutant_order{kk})>= dom_mutant_index_data
        
        indx_dom_mutant_data = find(msa_new(:,kk)==mutant_order{kk}(dom_mutant_index_data));
        msa_bin_data_dom_mutant(indx_dom_mutant_data,kk) = 1;
        
        indx_dom_mutant_sample = find(A(:,kk)==dom_mutant_index_data-2);
        msa_bin_sample_dom_mutant(indx_dom_mutant_sample,kk) = 1;
    end
end

sum_mut_data = sum(msa_bin_data_dom_mutant,2);
sum_mut_sample = sum(msa_bin_sample_dom_mutant,2);

x = 0:1:ls;

pmf_mut_data = histc(sum_mut_data,x);
pmf_mut_sample = histc(sum_mut_sample,x);

pmf_mut_data = pmf_mut_data./sum(pmf_mut_data);
pmf_mut_sample = pmf_mut_sample./sum(pmf_mut_sample);

figure;
semilogy(x,smooth(pmf_mut_data),'Color',blue,'LineWidth',1)
xlabel('Number of first most-prevalent mutant amino acids, \itr_1')
ylabel('Probability of {\itr_1}')
hold on
semilogy(x,smooth(pmf_mut_sample),'--','Color',blue,'LineWidth',1)
legend('Data','Model')
legend boxoff
xlim([0 20])
ylim([1e-4 1])

% [r1,pval1] = corr(pmf_mut_data,pmf_mut_sample,'type','spearman','tail','right')
% text(1,5e-5,'$$\rho_p = 0.7129$$','interpreter','latex')
% text(1,2.5e-5,'$$p \; \;= 3.56 \times 10^{-32}$$','interpreter','latex')

% [number,exponent] = convert_exp_power10(pval1);
% text(1,1e-3,sprintf('$$\\rho_s = %0.4f$$',r1),'interpreter','latex')
% text(1,4e-4,sprintf('$$p = %0.1f \\times 10^{%d}$$',number,exponent),'interpreter','latex')

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'YTick'       , [1e-4 1e-3 1e-2 1e-1 1],...
    'LineWidth'   , .5        );
set(gcf, 'PaperPositionMode', 'auto');

% % set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5.25 4])
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 4])
% % print('fig10','-depsc','-r300')
% print('fig10','-dpng','-r300')

%%
dom_mutant_index_data = 3; %first dom mutant = 2; second dom mutant = 3 and so on...
msa_bin_data_dom_mutant = zeros(ns,ls);
msa_bin_sample_dom_mutant = zeros(size(A));
% indx_dom_mutant = zeros(ns,ls);
for kk = 1:ls
    if length(mutant_order{kk})>= dom_mutant_index_data
        
        indx_dom_mutant_data = find(msa_new(:,kk)==mutant_order{kk}(dom_mutant_index_data));
        msa_bin_data_dom_mutant(indx_dom_mutant_data,kk) = 1;
        
        indx_dom_mutant_sample = find(A(:,kk)==dom_mutant_index_data-2);
        msa_bin_sample_dom_mutant(indx_dom_mutant_sample,kk) = 1;
    end
end

sum_mut_data = sum(msa_bin_data_dom_mutant,2);
sum_mut_sample = sum(msa_bin_sample_dom_mutant,2);

x = 0:1:ls;

pmf_mut_data = histc(sum_mut_data,x);
pmf_mut_sample = histc(sum_mut_sample,x);

pmf_mut_data = pmf_mut_data./sum(pmf_mut_data);
pmf_mut_sample = pmf_mut_sample./sum(pmf_mut_sample);

figure;
semilogy(x,smooth(pmf_mut_data),'Color',blue,'LineWidth',1)
xlabel('Number of first most-prevalent mutant amino acids, \itr_2')
ylabel('Probability of {\itr_2}')
hold on
semilogy(x,smooth(pmf_mut_sample),'--','Color',blue,'LineWidth',1)
legend('Data','Model')
legend boxoff
xlim([0 20])
ylim([1e-4 1])

% [r2,pval2] = corr(pmf_mut_data,pmf_mut_sample,'type','spearman','tail','right')
% [number,exponent] = convert_exp_power10(pval2);
% text(1,1e-3,sprintf('$$\\rho_s = %0.4f$$',r2),'interpreter','latex')
% text(1,4e-4,sprintf('$$p = %0.1f \\times 10^{%d}$$',number,exponent),'interpreter','latex')

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'YTick'       , [1e-4 1e-3 1e-2 1e-1 1],...
    'LineWidth'   , .5        );
% 
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 4])
% % set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5.25 4])
% %     print('fig11','-depsc','-r300')
% print('fig11','-dpng','-r300')

%%
% dom_mutant_index_data = 4; %first dom mutant = 2; second dom mutant = 3 and so on...
% msa_bin_data_dom_mutant = zeros(ns,ls);
% msa_bin_sample_dom_mutant = zeros(size(A));
% % indx_dom_mutant = zeros(ns,ls);
% for kk = 1:ls
%     if length(mutant_order{kk})>= dom_mutant_index_data
%         
%         indx_dom_mutant_data = find(msa_new(:,kk)==mutant_order{kk}(dom_mutant_index_data));
%         msa_bin_data_dom_mutant(indx_dom_mutant_data,kk) = 1;
%         
%         indx_dom_mutant_sample = find(A(:,kk)==dom_mutant_index_data-2);
%         msa_bin_sample_dom_mutant(indx_dom_mutant_sample,kk) = 1;
%     end
% end
% 
% sum_mut_data = sum(msa_bin_data_dom_mutant,2);
% sum_mut_sample = sum(msa_bin_sample_dom_mutant,2);
% 
% x = 0:1:ls;
% 
% pmf_mut_data = histc(sum_mut_data,x);
% pmf_mut_sample = histc(sum_mut_sample,x);
% 
% pmf_mut_data = pmf_mut_data./sum(pmf_mut_data);
% pmf_mut_sample = pmf_mut_sample./sum(pmf_mut_sample);
% 
% figure;
% semilogy(x,smooth(pmf_mut_data),'Color',blue,'LineWidth',1)
% xlabel('$$r_3$$','interpreter','latex')
% ylabel('$$p(r_3)$$','interpreter','latex')
% hold on
% semilogy(x,smooth(pmf_mut_sample),'--','Color',blue,'LineWidth',1)
% % legend('Data','Model')
% legend boxoff
% xlim([0 15])
% 
% [r3,pval3] = corr(pmf_mut_data,pmf_mut_sample,'type','spearman','tail','right')
% [number,exponent] = convert_exp_power10(pval3);
% text(5,1e-3,sprintf('$$\\rho_s = %0.4f$$',r3),'interpreter','latex')
% text(5,2.5e-4,sprintf('$$p = %0.3f \\times 10^{%d}$$',number,exponent),'interpreter','latex')
% 
% % text(7.5,1e-3,sprintf('$$\\rho_s = %.4f$$',r3),'interpreter','latex')
% % text(7.5,2.5e-4,'$$p = 2 \times 10^{-82}$$','interpreter','latex')
% 
% set(gca, ...
%     'Box'         , 'off'     , ...
%     'TickDir'     , 'out'     , ...
%     'TickLength'  , [.02 .02] , ...
%     'XMinorTick'  , 'on'      , ...
%     'YMinorTick'  , 'on'      , ...
%     'XColor'      , [.1 .1 .1], ...
%     'YColor'      , [.1 .1 .1], ...
%     'YTick'       , [1e-4 1e-3 1e-2 1e-1 1], ...
%     'LineWidth'   , .5        );
% 
% 
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 4])
% set(gcf,'renderer','painters')
% % print('fig12','-depsc','-r300')
% print('fig12','-dpng','-r300')
