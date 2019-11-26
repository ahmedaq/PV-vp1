function plot_mutdist(mutdist_data,mutdist_sampler,color)

numdata_data = xlsread(mutdist_data);
numdata_sampler = xlsread(mutdist_sampler);

%Normalization
numdata_data = numdata_data/sum(numdata_data);
numdata_sampler = numdata_sampler/sum(numdata_sampler);

%Linear Plot
% figure;
% plot(0:length(numdata_data)-1,numdata_data,color,'LineWidth',3)
% hold on
% plot(0:length(numdata_sampler)-1,numdata_sampler,sprintf('%s--',color),'LineWidth',3)
% % title(sprintf('Distribtion of mutations in %s',heading))
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$p(r)$$','interpreter','latex')
% legend('Data','Model')

%Semilogx Plot
figure;
semilogy(0:length(numdata_data)-1,(numdata_data),'Color',color,'LineWidth',1)
hold on
semilogy(0:length(numdata_sampler)-1,(numdata_sampler),'--','Color',color,'LineWidth',1)
% title(sprintf('Distribtion of mutations in %s',heading))
h = legend('Data','Model');
h.Position = [0.650 0.8524 0.1080 0.0488];
legend boxoff
xlabel('Number of mutations, \itr')
ylabel('Probability of {\itr} mutations')
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$p(r)$$','interpreter','latex')
% xlim([-0.1 no_of_nonzero_mutations(end)])
% ylim([min(numdata_data(numdata_data~=0)) 1])
xlim([0 25])
ylim([1e-4 1])

% no_of_nonzero_mutations = find(numdata_data ~= 0);
% figure;
% semilogy(0:length(numdata_data)-1,(numdata_data),'k','LineWidth',1)
% % set(gca,'YScale','log')
% hold on
% bar(0:length(numdata_sampler)-1,(numdata_sampler),1,'FaceColor',color,'EdgeColor',color)
% alpha 0.8
% % title(sprintf('Distribtion of mutations in %s',heading))
% h = legend('Data','Model');
% h.Position = [0.650 0.8524 0.1080 0.0488];
% legend boxoff
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$p(r)$$','interpreter','latex')
% % no_of_nonzero_mutations = find(numdata_data ~= 0);
% % xlim([-0.1 no_of_nonzero_mutations(end)])
% % ylim([min(numdata_data(numdata_data~=0)) 1])
% xlim([-0.5 24])
% ylim([6e-4 1])




%[r,pval] = corr(numdata_data(1:30)',numdata_sampler(1:30)','type','spearman','tail','right')
% text(2.5,3e-3,'$$\rho_s = 0.97$$','interpreter','latex')

% [number,exponent] = convert_exp_power10(pval);
% text(.6,1e-3,sprintf('$$\\rho_s = %0.4f$$',r),'interpreter','latex')
% text(.6,4e-4,sprintf('$$p = %0.1f \\times 10^{%d}$$',number,exponent),'interpreter','latex')
set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'      , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'YTick'       , [1e-4 1e-3 1e-2 1e-1 1], ...
        'LineWidth'   , .5        );
%     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5.25 4])

% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 4])
% % print(figname,'-depsc','-r300')
% print(figname,'-dpng','-r300')
