%%
% based on the conversion ot msa_aa_ex part of the
% Test990_RPROP_checkMCMC_without_k code. 

% takes an AA sequence and converts it to binary representation for the
% app. Potts model


% if input_seq contins an AA that does not occur in the mutant library 
% (mutant_order), then the AA of input_seq is replaced by the least
% occuring AA at that site (the last entry of mutant_order{x} for the xth
% site.
% hasNewAA = 1 if input_seq has AA at any site that do not appear in the
%              mutant_order library
% listOfSitesWithNewAA  is a list of sites hat have new AA not in 
%                       mutant_order library
function [out_seq_ex hasNewAA listOfSitesWithNewAA] = convertAAseq2Bin(mutant_order,diff_amino_site_length_RawMSA,bin_matrix,input_parm)

listOfSitesWithNewAA = [];
hasNewAA = false;
protein_length_aa = input_parm{1};
input_seq = input_parm{2};

testbin=1;
diff_amino_site=[];
curr_start_pos = 0;
for bcb=1:protein_length_aa
    diff_amino_temp = mutant_order{bcb};
    
    diff_length = diff_amino_site_length_RawMSA{bcb};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numcolumns = length(diff_length)-1;
    for cbc=1
        ind = find(input_seq(cbc,bcb)==diff_amino_temp);
        
        % if not present in diff_amino_temp, chose the last entry.i.e., the
        % least dominant AA at this site
        if(isempty(ind))
           ind = length(diff_amino_temp);
           hasNewAA = true;
           listOfSitesWithNewAA = [listOfSitesWithNewAA bcb];
        end
            
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
        
        out_seq_ex(cbc,curr_start_pos+1:curr_start_pos+numcolumns)=bin_value; % binary potts matrix
    end
    curr_start_pos =   curr_start_pos +numcolumns;

%     if (bcb>1) % construct start position in binary potts matrix wildtype mutnat model
%         phi_cumulative(bcb) = numcolumns+ sum(phi_cumulative(bcb-1));
%     else
%         phi_cumulative(bcb) = numcolumns;
%     end
%     phi_curr(bcb) = numcolumns; % number of amino acids per site
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
