function [msa,msa_aa_ex,phi_curr,phi_cumulative,cross_prod,mutant_order,true_indices,...
    total_length, protein_length_aa, diff_amino_site_length_RawMSA, bin_matrix] ...
    = generate_msa_binary_potts_msa_input(include_ratioIn,msa)

% include_ratioIn = 0.9; % 1: full Potts, 0: binary
threshCrit2 = 0.05;
%gammaIn = [1 2 5 7 10 20 ];%[1 2 2.5 3 4 5 5.5 6 7 8 9 10];
%gamma2In = [0 0.5 1 2 2.5 4 5];%[0:0.5:5 6 7 8 9 10];
gammaIn = 5.5%6%6%[140 180]%[100 120 150]%[1.5 30 ];%[1 2 2.5 3 4 5 5.5 6 7 8 9 10];
gamma2In = 0.095%1.3%1.4%[0.5]%[0 2.5]%[0 2 5];%[0:0.5:5 6 7 8 9 10];

protInd = 1:1000;%315;%315;%315;%100;% 10 10 10 10 10 10]
numSeqIn = -1;%820;%1000;%-1;%1000; % -1 means use all available sequences as input to Ray's code

protName = 'HA';%'NP';

strain = 32;%11;%11; %H1N!: 11,   H3N2:32
if(strcmp(protName, 'NP'))
    if(strain == 11)
        strainNameStr = 'H1N1';
        nsStr = '6510';
    elseif(strain == 32)
        strainNameStr = 'H3N2';
        nsStr = '4779';
    else
        display('*******************************************************************************')
        display('Error: Chose a valid strain name and correpsonding number of seqeunces.')
        display('*******************************************************************************')
        pause
    end
elseif(strcmp(protName, 'HA'))
    if(strain == 11)
        %strainNameStr = 'H1N1';
        %nsStr = '6510';
    elseif(strain == 32)
        strainNameStr = 'H3N2';
        nsStr = '16557';
    else
        display('*******************************************************************************')
        display('Error: Chose a valid strain name and correpsonding number of seqeunces.')
        display('*******************************************************************************')
        pause
    end
end

%strDataFileName = ['msa_ns' nsStr '_ls498_' protName '_' strainNameStr '_allSeq_ic_preProc_1536_1seqPerYrPerReg.mat'];
%strDataFileName = ['msa_ns' nsStr '_ls498_' protName '_' strainNameStr '_allSeq_ic_preProc_4867_1seqPerYrPerReg.mat'];
%strDataFileName = ['msa_ns' nsStr '_ls498_' protName '_' strainNameStr '_allSeq_ic_preProc_6403_1seqPerYrPerReg.mat'];

%strDataFileName = ['msa_ns' nsStr '_ls566_' protName '_' strainNameStr '_allSeq_ic_preProc_7436_1seqPerYrPerReg.mat'];
%strDataFileName = ['msa_ns' nsStr '_ls566_' protName '_' strainNameStr '_allSeq_ic_preProc_7405_1seqPerYrPerReg_3735.mat'];
% strDataFileName = 'msa_VP1_top_sec_Nseq_1633_Npos_302.mat';
% strDataFileName = 'msa_test.mat';

%strDataFileName = 'msa_ns16557_ls566_HA_H3N2_allSeq_ic_preProc_7404_nonFix_1seqPerYrPerReg_2707.mat';
gammaInLen = length(gammaIn);
gamma2InLen = length(gamma2In);


% this is for easy reference of final sort costFunc values and directly
% referes tot he numbers that appear on the last graph
firstCol = reshape(repmat(gammaIn, gamma2InLen,1), 1, gammaInLen*gamma2InLen)';
secondCol = repmat(gamma2In', gammaInLen,1);
gammaVec = [firstCol secondCol]';


varParmSingle = zeros(gammaInLen, gamma2InLen);
varParmDouble = zeros(gammaInLen, gamma2InLen);
perSingMutWithinLimts = zeros(gammaInLen, gamma2InLen);
perDoubMutWithinLimts = zeros(gammaInLen, gamma2InLen);
numSingMutPairs = zeros(gammaInLen, gamma2InLen);
numDoubMutPairs = zeros(gammaInLen, gamma2InLen);
input_parm = cell(1,20);

J_MINFLOW_array_init = 1;
j = 1
gammaInThis = gammaIn(j);
k = 1

input_parm{1} = strainNameStr;
input_parm{2} = protName;
input_parm{3} = numSeqIn;
input_parm{4} = protInd;
input_parm{5} = gammaInThis;
input_parm{6} = gamma2In(k);
% input_parm{7} = strDataFileName;
input_parm{8} = threshCrit2;
input_parm{9} = include_ratioIn;
input_parm{10} = J_MINFLOW_array_init;




% outputParm(1) = varParmSingle;
% outputParm(2) = varParmDouble;
% outputParm(3) = perSingMutWithinLimts;
% outputParm(4) = perDoubMutWithinLimts;
% outputParm(5) = numSingMutPairs;
% outputParm(6) = numDoubMutPairs;

strainName = input_parm{1};
protName = input_parm{2};
numSeqIn = input_parm{3};
protIndIn = input_parm{4};
gammaIn = input_parm{5};
gamma2In = input_parm{6};
strDataFileName = input_parm{7};
threshCrit2 = input_parm{8};
include_ratioIn = input_parm{9};
J_MINFLOW_array_init = input_parm{10};


protLenIn = length(protIndIn);
% load(strDataFileName);

%%%%%%%%%%%%%%%
% msa = msa_sec_top;

[indices_conserved_sites,no_of_conserved_sites] = find_conserved_sites(msa);
true_indices = setdiff(1:size(msa,2),indices_conserved_sites);
msa(:,indices_conserved_sites)=[];

if ~isempty(indices_conserved_sites)
    fprintf('\n Conserved sites detected')
end

allSeqLen = size(msa, 1);
msa_unique = msa;

%%%

if(numSeqIn == -1)
    selectedSeq = 1:allSeqLen;
else
    
    
    %load('msa_gp120_prob_D.mat') % original file Ray loaded
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section randomly choses numSeqIn unique "sequence indices" from the msa mat
    % file. Unique here means no index is repeated...if 2 different indices
    % have the same sequence, they will still be selected.
    msa_uniqueAll = msa_unique;
    total_num_seqIn = size(msa_uniqueAll,1);
    
    selectedSeq = unique(randi(total_num_seqIn, 1, numSeqIn));
    diff1 = numSeqIn - length(selectedSeq);
    
    while(diff1 > 0)
        temp1 = unique(randi(total_num_seqIn, 1, diff1));
        selectedSeq = [selectedSeq temp1];
        selectedSeq = unique(selectedSeq);
        
        diff1 = numSeqIn - length(selectedSeq);
        
    end
    % for sanity check..these 2 shoud be equal
    length(selectedSeq)
    length(unique(selectedSeq))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if(protLenIn > size(msa, 2))
    protLenIn = size(msa, 2);
    protIndIn = 1:protLenIn;
end
msa_unique = msa_unique(selectedSeq, protIndIn);



% 1 ------> this unique removes copies from the raw msa
aa_seq_str_unique = msa_unique;
%[aa_seq_str_unique index_unique]= unique(msa_unique,'rows');
% head_cell_unique = head_cell_new(index_unique);


%[aa_seq_str_unique index_unique]= unique(aa_seq_str_unique,'rows'); %????????????????????????????? why another unique?
% head_cell_unique = head_cell_unique(index_unique);

protein_length_aa = size(aa_seq_str_unique,2);
num_sequences = size(aa_seq_str_unique,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2 - Convert the amino acid to a wild-type, 1st mutant, 2nd mutant
% ...


% Part 2.1: Find the entropy of the 0/1, 0/1/2/ etc.. model, and find the
% sites to keep in each one

% relative ratio of amino acids to include relative amino acid entropy to keep
% 0 = full ising model
% 1 = full potts model
include_ratio=include_ratioIn;


diff_amino_site=[];
total_mutant=30;

keep_site_array=cell(1,total_mutant); % list of sites to include in the model
exclude_site_array=cell(1,total_mutant); % list of sites to include in the model
entropy_all = zeros(protein_length_aa,total_mutant); % entropy of the model for each site

%????????????????????????????? bcb and cbc?
for bcb=1:protein_length_aa
    %     workbar(bcb/protein_length_aa);
    curr_site_array = aa_seq_str_unique(:,bcb);
    diff_amino_temp = unique(curr_site_array);
    
    % find the frequency of amino acids at each site
    diff_length=[];
    for cbc=1:length(diff_amino_temp)
        curr_amino = diff_amino_temp(cbc);
        ind = find(curr_amino==curr_site_array);
        diff_length(cbc)=length(ind);  % will contain the count of all differenr AA at the (bcb)th site
    end
    
    
    
    [ysort indsort] = sort(diff_length,'descend'); % Sorted diff_length
    freq_aminos = ysort/num_sequences; % frequency of amino acids from most frequent to least frequent
    num_amino_curr = length(freq_aminos); % number of unique amino acids
    %%%%%%%%%%
    % Entropy of all amino acids
    freq_amino_array{bcb}=freq_aminos;
    entropy_amino(bcb) = -sum(freq_aminos.*log(freq_aminos));
    
    
    %???????????????????????? from where do you get these entropy sum and
    %???????????????????????? aba=1:min(max_mutant,num_amino_curr) does
    %this mean that 21 AA and gap are not considered?
    % what happens if 'Gaps' are within top 20? do u
    
    %all equations?
    %%%%%%%%%%%%%%%%%%%%
    % construct 0-1-2....jkj=1 (0/1 model), jkj=2 (0/1/2 model) etc...
    for jkj=1:total_mutant
        
        max_mutant=jkj+1;
        
        % entropy_sum = entropy for the "bcb"th site and 0/1/..jkj model
        entropy_sum(bcb)=0;
        for aba=1:min(max_mutant,num_amino_curr)
            %??????????????????????????? discuss with Ray
            if (aba==max_mutant)
                entropy_sum(bcb) = entropy_sum(bcb)-sum(freq_aminos(aba:end))*log(sum(freq_aminos(aba:end)));
                entropy_all(bcb,jkj) = entropy_all(bcb,jkj)-sum(freq_aminos(aba:end))*log(sum(freq_aminos(aba:end)));
            else
                entropy_sum(bcb) = entropy_sum(bcb)-freq_aminos(aba)*log(freq_aminos(aba));
                entropy_all(bcb,jkj)  = entropy_all(bcb,jkj)-freq_aminos(aba)*log(freq_aminos(aba));
            end
            
        end
        
        %??????????? dont understand...ask Ray
        relative_measure(bcb) = min(entropy_sum(bcb)/entropy_amino(bcb),entropy_amino(bcb)/entropy_sum(bcb));
        
        %??????????? whats the rationale behind keep sites and exclude
        %sites?
        if ( relative_measure(bcb)>=include_ratio)
            keep_site_array{jkj}=[keep_site_array{jkj} bcb]; % list of sites to include in "0/1/..jkj"th model
        else
            exclude_site_array{jkj}=[exclude_site_array{jkj} bcb]; % list of sites to exclude in "0/1/..jkj"th model
        end
        
    end
    %%%%%%%%%%%%%%%
    
end


% find the sites to keep in the 0/1/../jkj model
keep_leftover{1}=keep_site_array{1};
for jkj=2:total_mutant
    keep_leftover{jkj} = setdiff(keep_site_array{jkj},keep_site_array{jkj-1});
    
end

% arrayline=0:0.01:max(entropy_amino);
% figure
% plot(entropy_amino,entropy_all(:,1),'rx');hold;grid;
% plot(entropy_amino,entropy_all(:,2),'bx');
% plot(entropy_amino,entropy_all(:,3),'kx');
% plot(entropy_amino,entropy_all(:,4),'gx');
% plot(arrayline,arrayline,'b')
% legend('0/1','0/1/2','0/1/2/3','0/1/2/3/4','Location','SouthEast');
% xlabel('Entropy - all amino acids')
% ylabel('Entropy')
%
% figure
% plot(entropy_amino(keep_leftover{1}),entropy_all(keep_leftover{1},1),'rx');hold;grid;
% plot(entropy_amino(keep_leftover{2}),entropy_all(keep_leftover{2},2),'bx');
% plot(entropy_amino(keep_leftover{3}),entropy_all(keep_leftover{3},3),'kx');
% plot(entropy_amino(keep_leftover{4}),entropy_all(keep_leftover{4},4),'gx');
% plot(arrayline,arrayline,'b')
% legend('0/1','0/1/2','0/1/2/3','0/1/2/3/4','Location','SouthEast');
% xlabel('Entropy - all amino acids')
% ylabel('Entropy of kept amino acids')

% Part 2.2 - Convert to binary MSA

% Part 2.2.1 Create new amino acid sequence "aa_new" based on the above
% 0/1/... model for each site. For example, if 0/1 model, replace all
% second dominant mutants with the most dominant mutant

aa_new = aa_seq_str_unique;
for bcb=1:protein_length_aa
    %     workbar(bcb/protein_length_aa);
    for aba=1:length(keep_leftover)
        ind2 = find(bcb==keep_leftover{aba});
        if (length(ind2)>0)
            mutant_value=aba;
        end
    end
    
    curr_site_array = aa_seq_str_unique(:,bcb);
    diff_amino_temp = unique(curr_site_array);
    diff_length=[];
    for cbc=1:length(diff_amino_temp)
        curr_amino = diff_amino_temp(cbc);
        ind = find(curr_amino==curr_site_array);
        diff_length(cbc)=length(ind);
    end
    
    [ysort indsort] = sort(diff_length,'descend');
    
    sort_amino = diff_amino_temp(indsort);
    
    mutants_to_change= sort_amino(mutant_value+1:end);
    
    for aba=1:length(mutants_to_change)
        ind3=find(mutants_to_change(aba)==curr_site_array);
        aa_new(ind3,bcb)=sort_amino(mutant_value+1);
    end
    mutant_order{bcb}=sort_amino(1:mutant_value+1);
    
end
aa_seq_str_unique=aa_new;


%aa_seq_str_unique = unique(aa_seq_str_unique,'rows');
num_sequences_unique = size(aa_seq_str_unique,1);
protein_length_aa =  size(aa_seq_str_unique,2);

% Part 2.2.2  Construct binary MSA - wildtype/1 mutant, 2 mutant model, and related parameteres

for aba=1:total_mutant
    temp_matrix=[];
    temp_matrix = fliplr(eye(aba));
    temp_matrix = [zeros(1,size(temp_matrix,2)) ; temp_matrix ];
    bin_matrix{aba} =temp_matrix;
    
end

testbin=1;

diff_amino_site=[];
curr_start_pos = 0;
for bcb=1:protein_length_aa
    %     workbar(bcb/protein_length_aa);
    curr_site_array = aa_seq_str_unique(:,bcb);
    diff_amino_temp = unique(curr_site_array);
    diff_length=[];
    for cbc=1:length(diff_amino_temp)
        curr_amino = diff_amino_temp(cbc);
        ind = find(curr_amino==curr_site_array);
        diff_length(cbc)=length(ind);
    end
    
    [ysort indsort] = sort(diff_length,'descend');
    diff_amino_temp=diff_amino_temp(indsort);
    
    diff_amino_temp = mutant_order{bcb};
    
    
    diff_amino_site_length{bcb} = indsort;
    diff_amino_site_length_RawMSA{bcb} = ysort;
    diff_amino_site{bcb}=diff_amino_temp;% different amino acids per site
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numcolumns = length(diff_length)-1;
    for cbc=1:num_sequences_unique
        ind = find(aa_seq_str_unique(cbc,bcb)==diff_amino_temp);
        if (ind==1)
            num_index=0;
        else
            num_index = 2^(ind-2);
        end
        
        
        if testbin==1
            curr_bin_matrix = bin_matrix{numcolumns};
            bin_value = curr_bin_matrix(ind,:);
            %             bin_value2 = fliplr(de2bi(num_index,numcolumns));
            %             diff=bin_value-bin_value2;
            %             if (diff~=0)
            %                 here=1;
            %             end
        else
            bin_value = fliplr(de2bi(num_index,numcolumns));
        end
        %         bin_value = fliplr(bitget(num_index,1:numcolumns));
        %         bin_value = bitget(num_index,1:numcolumns);
        
        msa_aa_ex(cbc,curr_start_pos+1:curr_start_pos+numcolumns)=bin_value; % binary potts matrix
    end
    curr_start_pos =   curr_start_pos +numcolumns;
    if (bcb>1) % construct start position in binary potts matrix wildtype mutnat model
        phi_cumulative(bcb) = numcolumns+ sum(phi_cumulative(bcb-1));
    else
        phi_cumulative(bcb) = numcolumns;
    end
    phi_curr(bcb) = numcolumns; % number of amino acids per site
end

cross_prod_site = (msa_aa_ex')*msa_aa_ex/num_sequences_unique;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%
% msa_aa_ex2=msa_aa_ex;
% msa_aa_ex = msa_aa_ex2(1:5,:);

for aba=1:30
    temp_matrix=[];
    temp_matrix = fliplr(eye(aba));
    temp_matrix = [zeros(1,size(temp_matrix,2)) ; temp_matrix ];
    bin_matrix{aba} =temp_matrix;
    
end

testbin=1;
total_length = size(msa_aa_ex,2);
num_sequences_unique = size(msa_aa_ex,1);
% Xtotal=sparse(zeros(total_length*num_sequences_unique,total_length));
% Dtotal=sparse(zeros(total_length*num_sequences_unique,total_length));
% Dtotal=[];
for aba=1:num_sequences_unique % go through each sequence
    curr_vector = msa_aa_ex(aba,:);
    %     workbar(aba/num_sequences_unique);
    count=1;
    count2=1;
    
    clear d;
    for bcb=1:protein_length_aa % go through each site
        length_curr = phi_curr(bcb);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        % all flip only to dominant mutant
        % find the starting position of the bcbth site
        if (bcb>1)
            curr_start = phi_cumulative(bcb-1)+1;
        else
            curr_start=1;
        end
        pos_data=curr_start:curr_start+length_curr-1;
        curr_data = curr_vector(curr_start:curr_start+length_curr-1);
        
        if testbin==1
            if (curr_data==0)
                dec_value=0;
            elseif (curr_data==1)
                dec_value=1;
            else
                dec_value = bi2de(fliplr(curr_data));
            end
            
            %                   dec_value2 = bi2de(fliplr(curr_data));
            %             diff=dec_value-dec_value2;
            %             if (diff~=0)
            %                 here=1;
            %             end
        else
            dec_value = bi2de(fliplr(curr_data));
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         count=count+size(add_bi_values,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bin_value = dec_value;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        jj=1:length_curr;
        con_values = [0 2.^(jj-1)];
        indnow = find(bin_value==con_values);
        jjcurrvalues=[0 jj];
        jjcurrvalues(indnow)=[];
        con_values2=con_values;
        con_values(indnow) = [];
        
        %         add_bi_values =   fliplr(de2bi(con_values,length_curr));
        if (testbin==1)
            curr_bin_matrix = bin_matrix{length_curr};
            add_bi_values = curr_bin_matrix(jjcurrvalues+1,:);
            
            %                           add_bi_values2 =   fliplr(de2bi(con_values,length_curr));
            %             diff=sum(sum(abs(add_bi_values-add_bi_values2)));;
            %             if (diff~=0)
            %                 here=1;
            %             end
            
        else
            add_bi_values =   fliplr(de2bi(con_values,length_curr));
        end
        
        
        
        %          add_bi_values = bitget(con_values,1:length_curr);
        indone = find(curr_data==1);
        if (sum(curr_data)>0)
            add_bi_values(:,indone)=-1;
        end
        
        d(count:count+size(add_bi_values,1)-1,curr_start:curr_start+length_curr-1)=add_bi_values;
        count=count+size(add_bi_values,1);
        %%%%%%%%%%%%%%%%%%
        
        
        
    end
    
    
    d_store1{aba}=sparse(d);
    %     d_store_temp{aba}=sparse(d);
    x_store_temp{aba} = sparse(bsxfun(@plus,curr_vector,d/2));
    
    %     Xtotal((aba-1)*total_length+1:(aba-1)*total_length+total_length,:)=sparse(bsxfun(@plus,curr_vector,d/2));
    %     Dtotal((aba-1)*total_length+1:(aba-1)*total_length+total_length,:) = sparse(d);
    %     Xtotal = [Xtotal ; sparse(bsxfun(@plus,2*curr_vector,d))];
    %     Dtotal=[Dtotal; sparse(d)];
    
end
Xtotal = sparse(cell2mat(x_store_temp'));
Dtotal = sparse(cell2mat(d_store1'));
% total_x_length =num_sequences_unique*total_length/3;

% pause
countstore=count-size(add_bi_values,1);
% dTime_construct = toc(t_samp);
xarray=[];
for aba=1:num_sequences_unique
    curr_vector = msa_aa_ex(aba,:);
    ind = find(curr_vector==1);
    xarray=[xarray ind];
    xpos(aba)=length(ind);
    if (aba==1)
        xstartpos(aba)=1;
    else
        xstartpos(aba)=sum(xpos(1:aba-1))+1;
    end
end
num_d2 = ones(1,length(d_store1))*size(d,1);


%% Finding no of mutations per sequence

no_of_mutations_per_seq = sum(msa_aa_ex,2);
x = 0:1:allSeqLen;

no_of_mutations_in_each_bin = histc(no_of_mutations_per_seq,x);

pmf_mutations = no_of_mutations_in_each_bin./sum(no_of_mutations_in_each_bin);

figure;
semilogy(x,pmf_mutations)
xlabel('Number of mutations')
ylabel('pmf')

%Saving pmf_mutations

% %in mat
% var_name_save = sprintf('%s-mutdist-data',protein); 
% save var_name_save pmf_mutations
% 
% %in xlsx format
% xlswrite(sprintf('%s-mutdist-data.xlsx',protein),pmf_mutations.')


%% Finding one and two point correlations

cross_prod = ((msa_aa_ex')*msa_aa_ex)/size(msa_aa_ex,1); 
%diag = 1 pt mutation
%upper triangle = 2 pt mutations

%%

%%gettng one point mutation from actual MSA (not binary)
% [ns,ls] = size(msa);
% one_point_mutations = zeros(ls,max(phi_curr));
% for kk = 1:length(mutant_order)    
%     for mm = 1:length(mutant_order{kk})-1 
%        one_point_mutations(kk,mm)= sum(msa(:,kk)==mutant_order{kk}(mm+1))/ns; %first one is wildtype
%     end
%     
% end

% %% Saving one n two point mutations in a format compatible with John's algo
% 
% [ns,ls] = size(msa);
% [ns_ex,ls_ex] = size(msa_aa_ex);
% 
% %saving one point mutations 
% index = 1;
% for kk = 1:ls
%     indices = index;
%     for mm = 2:phi_curr(kk)        
%         indices = [indices indices(end)+1];       
%     end    
%     dlmwrite('correlations_John_VP1.txt',diag(cross_prod(indices,indices)).','-append',....
%         'delimiter','\t','precision', '%.8e', 'newline', 'pc', 'roffset', 0, 'coffset', 0)
%     index = indices(end)+1;
% end
% 
% %saving two point mutations
% index = 1;
% no_two_point_mutations = ls*(ls-1)/2
% for kk = 1:ls
%     indices = index;
%     for mm = 2:phi_curr(kk)        
%         indices = [indices indices(end)+1];       
%     end
%     
%     indices_2 = indices(end)+1;
%     for nn = kk+1:ls
%         for pp = 2:phi_curr(nn)
%             indices_2 = [indices_2 indices_2(end)+1];
%         end
%         temp = cross_prod(indices_2,indices);
% %         temp = cross_prod(indices_2,indices).';
%         temp = (temp(:)).';
%         dlmwrite('correlations_John_VP1.txt',temp,'-append',....
%         'delimiter','\t','precision', '%.8e', 'newline', 'pc', 'roffset', 0, 'coffset', 0)        
%         indices_2 = indices_2(end)+1;
%     end
% end

