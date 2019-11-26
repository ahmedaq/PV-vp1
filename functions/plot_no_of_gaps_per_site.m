function no_of_gaps_per_site = plot_no_of_gaps_per_site(msa,heading,flag)

[Nseq,Npos] = size(msa);
no_of_gaps_per_site = [];
no_of_stars_per_site= [];
no_of_B_per_site= [];
no_of_X_per_site= [];
for k = 1:Npos
    no_of_gaps_per_site(k) = length(find(msa(:,k) == '-'));
    no_of_stars_per_site(k) = length(find(msa(:,k) == '*'));
    no_of_B_per_site(k) = length(find(msa(:,k) == 'B'));
    no_of_X_per_site(k) = length(find(msa(:,k) == 'X'));
end

if flag == 1
    figure;bar(1:Npos,no_of_gaps_per_site)
    title(sprintf('Number of gaps per site in %s',heading))
end