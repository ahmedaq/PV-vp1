


function [Vall samples_MCMC_double random_array_sweep]  = verifyMCMC_TC_PT_noMSA(param_verifyMCMC_TC_PT)


J_MINFLOW_mat = param_verifyMCMC_TC_PT{1};
total_length = param_verifyMCMC_TC_PT{2};
msa_aa_ex = param_verifyMCMC_TC_PT{3};
phi_curr = param_verifyMCMC_TC_PT{4};
phi_cumulative = param_verifyMCMC_TC_PT{5};
protein_length_aa = param_verifyMCMC_TC_PT{6};
strMethod = param_verifyMCMC_TC_PT{7};
strSmallorLarge = param_verifyMCMC_TC_PT{8};
temp12 = param_verifyMCMC_TC_PT{9};
temp22 = param_verifyMCMC_TC_PT{10};
seedSeq = param_verifyMCMC_TC_PT{11};

thin = temp12(1);
burnin = temp12(2);
nosim = temp12(3);

MCsweepLength = temp22(1);
numParallel = temp22(2);
betaArray = temp22(3:end);



clear param_verifyMCMC;



%% verfication

random_array_sweep = rand(1,nosim/MCsweepLength*(numParallel - 1));
    number_samples = ceil((nosim-burnin)/thin);
    
    
J_MINFLOW_mat_array = J_MINFLOW_mat(:);
t_samp = tic();

random_site_array =  randi([1 protein_length_aa],1,nosim*numParallel); % choose the random site

num_amino_array = phi_curr(random_site_array);
unique_amino = unique(num_amino_array);
for cbc=1:length(unique_amino)
    num_amino = unique_amino(cbc);
    ind = find(num_amino_array==unique_amino(cbc));
    rand_amino_array(ind) = randi([1 num_amino],1,length(ind));
end
random_array = rand(1,nosim*numParallel);
t_samp = toc(t_samp);
fprintf( 'Random generation in %f seconds \n', t_samp );


totalnosample=0;

t_samp = tic();
samples_MCMC_double = zeros(number_samples,total_length);

aba=1;

curr_vector = seedSeq;
    
[doublemutant nosample energyAll numMutAll Vall samples_MCMC ]= PT_final_chain_1(random_array,...
    random_site_array,rand_amino_array,curr_vector,J_MINFLOW_mat_array,total_length,...
    nosim,phi_cumulative,phi_curr,burnin,thin,curr_vector,number_samples, numParallel, ...
    betaArray(1:numParallel),MCsweepLength, random_array_sweep);
    
samples_MCMC_double((aba-1)*number_samples+1:(aba-1)*number_samples+number_samples,:)=...
    reshape(samples_MCMC,total_length,nosample)';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_samp = toc(t_samp);
thirdMCMC=t_samp;
fprintf( 'MCMC in %f seconds \n', t_samp );

cross_prod_site= (msa_aa_ex')*msa_aa_ex/size(msa_aa_ex,1);

double_mutant_sample = (samples_MCMC_double')*samples_MCMC_double/size(samples_MCMC_double,1);
single_mutant_sample = mean(samples_MCMC_double);

endPointThisLine = max([cross_prod_site(:) double_mutant_sample(:)]);

arrayline = [0:0.01:endPointThisLine];%+0.05;

cross_nondiag = cross_prod_site - diag(diag(cross_prod_site));
total_nondiag = double_mutant_sample - diag(diag(double_mutant_sample));

cross_diag = diag(cross_prod_site);
total_diag = diag(double_mutant_sample);
    
numpatients = size(msa_aa_ex,1);

%%%%%%%%%

figure
if(strcmpi(strSmallorLarge,'small'))
    subplot(2,2,1)
end
plot(cross_nondiag(:),total_nondiag(:),'bx');hold;
plot(arrayline,arrayline,'k'), axis([0 1 0 1]);
xlabel('Double Mutant MSA')
ylabel('Double Mutant MCMC')

if(strcmpi(strSmallorLarge,'large'))
   figure
else
  subplot(2,2,3)
end
plot(cross_diag(:),total_diag(:),'rx');hold on;
plot(arrayline,arrayline,'k'), axis([0 1 0 1]);
xlabel('Single Mutant MSA')
ylabel('Single Mutant MCMC')

[xMutMtxCell_actualMSA, allMutMtx_actualMSA] = mutCountMSA(msa_aa_ex, phi_curr, protein_length_aa);
[xMutMtxCell_MCMC, allMutMtx_MCMC] = mutCountMSA(samples_MCMC_double, phi_curr, protein_length_aa);

[freqCountNumMutPerSeq_actualMSA numSeq_actualMSA] = plotXMutFreqPMF(allMutMtx_actualMSA, protein_length_aa);
[freqCountNumMutPerSeq_MCMC numSeq_MCMC] = plotXMutFreqPMF(allMutMtx_MCMC, protein_length_aa);

pdf_actualMSA = freqCountNumMutPerSeq_actualMSA/numSeq_actualMSA./(sum(freqCountNumMutPerSeq_actualMSA/numSeq_actualMSA));
pdf_MCMC = freqCountNumMutPerSeq_MCMC/numSeq_MCMC./(sum(freqCountNumMutPerSeq_MCMC/numSeq_MCMC));

if(strcmpi(strSmallorLarge,'large'))
   figure
else
  subplot(2,2,2)
end

figure
semilogy(0:protein_length_aa, (pdf_actualMSA), 'rx-'), hold on
semilogy(0:protein_length_aa, (pdf_MCMC), '.-')
legend('MSA', 'MCMC')
title(['All mutation PMF method: ' strMethod])
xlabel('number of mutations')
ylabel('PMF')
axis([0 200 0 max(max(pdf_MCMC),max(pdf_actualMSA))+.01]);%axis([0 150 0 .25])%axis([0 20 0 .35])%axis([0 150 0 .15])

if(strcmpi(strSmallorLarge,'large'))
   figure
else
  subplot(2,2,4)
end

plot(sum(samples_MCMC_double,2))
xlabel('MCMC Run')
ylabel('Number of mutants')
