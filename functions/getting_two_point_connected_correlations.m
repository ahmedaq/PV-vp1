function [p2_data p2_sampler] = getting_two_point_connected_correlations(input_configurations_sampler,input_data_file)

% The output file will have the extension “.p”, but it contains Potts
% configurations rather than correlations. Each line is one configuration.
% The list of numbers in each line gives the Potts state at each site, with
% the convention that 0 is the first mutant state. The “wild-type” state at
% each site i actually denoted by q_i-1, where q_i is the total number of
% states at that site. In other words, if there are three states at a site
% (and thus two 1-point correlations or fields in the .p and .j files),
% the states represent:
%
% 0 - first mutant
% 1 - second mutant
% 2 - wild-type

A = xlsread(input_configurations_sampler);

load(input_data_file)

[ns,ls] = size(msa);

%Making the modified msa according to mutant_order

msa_new = construct_msa_aa_after_entropy_compression(msa,mutant_order);


p2_data = [];
p2_sampler = [];

parfor kk = 1:ls-1
    for mm = kk+1:ls
        p22_data = [];
        p22_sampler = [];
        for kkk = 1:(length(mutant_order{kk})-1)
            for mmm = 1:(length(mutant_order{mm})-1)
                temp_data = msa_new(:,[kk mm]) == ...
                    repmat([mutant_order{kk}(kkk+1) mutant_order{mm}(mmm+1)],ns,1);
                p22_data = [p22_data sum((temp_data(:,1).*temp_data(:,2)))/ns - ...
                    ((sum(temp_data(:,1))/ns)*(sum(temp_data(:,2))/ns))];
                
                temp = A(:,[kk mm]) == repmat([kkk-1 mmm-1],size(A,1),1);
                p22_sampler = [p22_sampler sum((temp(:,1).*temp(:,2)))/size(A,1) - ...
                    ((sum(temp(:,1))/size(A,1))*(sum(temp(:,2))/size(A,1)))];
            end
        end
        p2_data = [p2_data fliplr(p22_data)];
        p2_sampler = [p2_sampler p22_sampler];
    end
end