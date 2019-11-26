function H = ConstructHmatrix_couplings(name,ls,phi_curr,phi_cumulative)

% A = textread(name,'','delimiter','\t','emptyvalue', NaN);
A = xlsread(name); %output of John's RPROP code (opened and saved in xlsx format (important!))

% ls = size(msa,2); %Number of mutating sites (can be found in VP1_Nseq_1633_Npos_302_S0.90.rep)
%     ls = 20;
number_mutants_site = zeros(1,ls);
for kk = 1:ls
    number_mutants_site(kk) = length(find(~isnan(A(kk,:)))); %wont work for dlmread
end

%Putting the h values on the diagonal of H matrix
total_mutants = sum(number_mutants_site); %This should be equal to ls_ex
H = zeros(total_mutants,total_mutants);
mm = 0;
for kk = 1:ls
    %     if kk == 1
    %         indx = 1:number_mutants_site(kk);
    %     else
    indx = mm+1:mm+number_mutants_site(kk);
    H(indx,indx) = diag(A(kk,find(~isnan(A(kk,:)))));
    mm = mm+number_mutants_site(kk);
end

phi_cum = [0 phi_cumulative];
no_two_point_mutations = ls*(ls-1)/2;

for ss = ls+1:no_two_point_mutations+ls
    a{ss} = A(ss,find(~isnan(A(ss,:))));
end

ss = ls+1;
for kk = 1:ls
    indices = phi_cum(kk)+1;
    for mm = 2:phi_curr(kk)
        indices = [indices indices(end)+1];
    end
    indices_2 = indices(end)+1;
    for nn = kk+1:ls
        for pp = 2:phi_curr(nn)
            indices_2 = [indices_2 indices_2(end)+1];
        end
        indx1 = indices;
        indx2 = indices_2;
        ep = 1;
        for k = 1:length(indx1)
            xxx = indx1(k);
            for m = 1:length(indx2)
                yyy = indx2(m);
                H(xxx,yyy) = a{ss}(ep);
                ep = ep+1;
            end
        end
        ss=ss+1;
        indices_2 = indices_2(end)+1;
    end
    
end

