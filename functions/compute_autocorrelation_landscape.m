
function ak = compute_autocorrelation_landscape(data_protein, H_protein, D)

load(data_protein)
load(H_protein)

[~,ls_msa] = size(msa);

%Parameters
nrs = 1e6;  %number of random start sequences
% D = 5;      %maximum number of mutations in random start sequences
no_of_steps = 50; %the k value that we need to find correlation over

[~,ak] = correlation_calculation(nrs,D,no_of_steps,...
    ls_msa,msa_aa_ex,phi_curr,phi_cumulative,mutant_order,H);