function pvalue = pvalue_counting(total_no_of_sites,no_of_sites_in_sector,sites_test,sites_test_in_sector)
warning off
% total_no_of_sites = 631;
% no_of_sites_in_sector = 246;
% 
% sites_test = [105;60;23;61;10;70];
% sites_test_in_sector = [49;30;15;33;4;44];

pvalue = zeros(1,length(sites_test));

for j = 1:length(sites_test)    
    for k = sites_test_in_sector(j):1:min(sites_test(j),no_of_sites_in_sector)
        pvalue(j) = pvalue(j) + ...
            (nchoosek(sites_test(j),k)*...
            nchoosek(total_no_of_sites-sites_test(j),no_of_sites_in_sector-k))/nchoosek(total_no_of_sites,no_of_sites_in_sector);
    end
end