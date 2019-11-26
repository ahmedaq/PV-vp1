function save_correlations_ACE(msa_aa_ex,phi_curr,phi_cumulative,protein)

%% Finding one and two point correlations

cross_prod = ((msa_aa_ex')*msa_aa_ex)/size(msa_aa_ex,1); 
%diag = 1 pt mutation
%upper triangle = 2 pt mutations

%%gettng one point mutation from actual MSA (not binary)
% [ns,ls] = size(msa);
% one_point_mutations = zeros(ls,max(phi_curr));
% for kk = 1:length(mutant_order)    
%     for mm = 1:length(mutant_order{kk})-1 
%        one_point_mutations(kk,mm)= sum(msa(:,kk)==mutant_order{kk}(mm+1))/ns; %first one is wildtype
%     end
%     
% end

%% Saving one n two point mutations in a format compatible with John's algo

% [ns,ls] = size(msa);
ls = length(phi_curr);
[ns_ex,ls_ex] = size(msa_aa_ex);

%saving one point mutations 
index = 1;
for kk = 1:ls
    indices = index;
    for mm = 2:phi_curr(kk)        
        indices = [indices indices(end)+1];       
    end    
    dlmwrite(sprintf('correlations_ACE_%s.p',protein),diag(cross_prod(indices,indices)).','-append',....
        'delimiter','\t','precision', '%.8e', 'newline', 'pc', 'roffset', 0, 'coffset', 0)
    index = indices(end)+1;
end

%saving two point mutations
phi_cum = [0 phi_cumulative];
no_two_point_mutations = ls*(ls-1)/2;
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
        temp = cross_prod(indices_2,indices); %confirmed from John
%         temp = cross_prod(indices_2,indices).'; 
        temp = (temp(:)).';
        dlmwrite(sprintf('correlations_ACE_%s.p',protein),temp,'-append',....
        'delimiter','\t','precision', '%.8e', 'newline', 'pc', 'roffset', 0, 'coffset', 0)        
        indices_2 = indices_2(end)+1;
    end
end
