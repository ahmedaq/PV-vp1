function [CriterionValuesSil, CriterionValuesCH, CriterionValuesDB, ...
    OptimalSil, OptimalCH, OptimalDB] ...
    = evaClusters(A,labels)

eva = evalclusters(A,labels,'silhouette');%,'Distance','correlation');
CriterionValuesSil = eva.CriterionValues;
OptimalSil = eva.OptimalK;
figure;
subplot(311)
plot(1:size(labels,2),eva.CriterionValues)
title(sprintf('silhouette, optimal number of clusters = %d',eva.OptimalK))
xlabel('Number of clusters')
ylabel('Value')

eva = evalclusters(A,labels,'CalinskiHarabasz');
OptimalCH = eva.OptimalK;
CriterionValuesCH = eva.CriterionValues;
subplot(312)
plot(1:size(labels,2),eva.CriterionValues)
title(sprintf('CalinskiHarabasz, optimal number of clusters = %d',eva.OptimalK))
xlabel('Number of clusters')
ylabel('Value')

eva = evalclusters(A,labels,'DaviesBouldin');
OptimalDB = eva.OptimalK;
CriterionValuesDB = eva.CriterionValues;
subplot(313)
plot(1:size(labels,2),eva.CriterionValues)
title(sprintf('DaviesBouldin, optimal number of clusters = %d',eva.OptimalK))
xlabel('Number of clusters')
ylabel('Value')