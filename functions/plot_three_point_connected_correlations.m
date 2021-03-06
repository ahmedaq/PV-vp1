function plot_three_point_connected_correlations(p3_data,p3_sampler,color,markersize)

thresh_zero = 0;

figure;
plot(p3_data(abs(p3_data)>thresh_zero),p3_sampler(abs(p3_data)>thresh_zero),'o','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',markersize)
hold on;
max_value = max(max(p3_data,p3_sampler));
min_value = min(min(p3_data,p3_sampler));
plot((min_value-0.1):max_value/20:(max_value+0.1),(min_value-0.1):max_value/20:(max_value+0.1),'k')
xlabel('Three-point connected correlations (MSA)')
ylabel('Three-point connected correlations')

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
