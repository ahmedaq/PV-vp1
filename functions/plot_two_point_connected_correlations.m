function plot_two_point_connected_correlations(p2_data,p2_sampler,color,markersize)

thresh_zero = 0;

figure;
plot(p2_data(abs(p2_data)>thresh_zero),p2_sampler(abs(p2_data)>thresh_zero),'o','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',markersize)
hold on;
max_value = max(max(p2_data,p2_sampler));
min_value = min(min(p2_data,p2_sampler));
plot((min_value-0.1):max_value/20:(max_value+0.1),(min_value-0.1):max_value/20:(max_value+0.1),'k')
xlabel('Two-point connected correlations (MSA)')
ylabel('Two-point connected correlations')

% [r,pval] = corr(p2_data(abs(p2_data)>thresh_zero).',p2_sampler(abs(p2_data)>thresh_zero).','type','pearson','tail','right')

xlim([-0.1 0.2])
ylim([-0.1 0.2])
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...
  'YTick'       , [-0.15:0.05:0.3],...
  'XTick'       , [-0.15:0.05:0.3],...
  'LineWidth'   , .5        );
