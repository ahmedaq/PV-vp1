
function energy_vs_fitness_nonpeak1(data_vp1_nonpeak1,H_vp1_nonpeak1)

run startup.m

% load data_vp1_v5_peak1
% load H_vp1_v5_peak1
load(data_vp1_nonpeak1)
load(H_vp1_nonpeak1)

Cseq = seqconsensus(msa);

%% Check to see if mutants in Mahoney strain exist in MSA!

load vp1_mahoney
vp1_mahoney_mut = vp1_mahoney(true_indices);
[~,alignment_mahoney] = nwalign(vp1_mahoney_mut,Cseq);
sites_diff_mahoney = find(alignment_mahoney(2,:)~='|');
for kk = 1:length(sites_diff_mahoney)
    seqs_with_mutants_as_in_mahoney{kk} = find(msa(sites_diff_mahoney(kk),:)==vp1_mahoney_mut(sites_diff_mahoney(kk)));
end

%% Shulman et al. 2014

site_freq_no_mutation = setdiff(1:302,true_indices);

load('seqs_Shulman2014.mat') %seqs from supplement of paper

msaToTestComplete = [vp1_mahoney;seqs_Shulman2014];

msaToTest = msaToTestComplete(:,true_indices);

numSeqMsaToTest = size(msaToTest, 1);
out_seq_ex = zeros(numSeqMsaToTest, total_length);
EnergySeq_Shulman2014 = zeros(1, numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    
    input_seq = msaToTest(i,:);

    input_parm = cell(1,2);
    input_parm{1} = protein_length_aa;
    input_parm{2} = input_seq;

    [out_seq_ex(i,:)] = convertAAseq2Bin(mutant_order,diff_amino_site_length_RawMSA,bin_matrix,input_parm);
    EnergySeq_Shulman2014(i) = calcSeqEnergy(out_seq_ex(i,:),H);
    
    indices_diff = find(vp1_mahoney~=msaToTestComplete(i,:));
    no_conserved_sites_in_indices_diff = sum(ismember(indices_diff,site_freq_no_mutation));
%     if i==2 
    EnergySeq_Shulman2014(i) = EnergySeq_Shulman2014(i) - no_conserved_sites_in_indices_diff*min(diag(H));
%     end
end

%taking inverse because it is PD50, the higher the value, the less virulence
pd50_neurovirulence = 1./[5.9 7.2 6.7 6.8];
pd50_neurovirulence_norm = (pd50_neurovirulence/pd50_neurovirulence(1));

% [rho,pval] = corr([EnergySeq_Shulman2014-EnergySeq_Shulman2014(1)]',pd50_neurovirulence_norm(:),'tail','left','type','spearman')

%% Colston et al. 1995

vp1_mahoney_P95S = vp1_mahoney;
mutant_order{rev_translation_indices(95,true_indices)};
vp1_mahoney_P95S(95) = 'S';

vp1_mahoney_P95T = vp1_mahoney;
mutant_order{rev_translation_indices(95,true_indices)};
vp1_mahoney_P95T(95) = 'T';

vp1_mahoney_V160I = vp1_mahoney;
mutant_order{rev_translation_indices(160,true_indices)};
vp1_mahoney_V160I(160) = 'I';

vp1_mahoney_P95S_V160I = vp1_mahoney;
vp1_mahoney_P95S_V160I([95 160]) = 'SI';

% msaToTestComplete = [vp1_mahoney;vp1_mahoney_T99K]
msaToTestComplete = [vp1_mahoney;vp1_mahoney_P95S;vp1_mahoney_P95T;vp1_mahoney_V160I;vp1_mahoney_P95S_V160I];

% msaToTest = [Cseq;Cseq_99K;Cseq_99T;msaToTestComplete(:,true_indices)];
msaToTest = msaToTestComplete(:,true_indices);
% msaToTest = Cseq;

numSeqMsaToTest = size(msaToTest, 1);
out_seq_ex = zeros(numSeqMsaToTest, total_length);
EnergySeq = zeros(1, numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    
    input_seq = msaToTest(i,:);

    input_parm = cell(1,2);
    input_parm{1} = protein_length_aa;
    input_parm{2} = input_seq;

    [out_seq_ex(i,:)] = convertAAseq2Bin(mutant_order,diff_amino_site_length_RawMSA,bin_matrix,input_parm);
    [EnergySeq(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end

titer = [1.6e9 1.8e9 5.5e8 7.8e8 8.3e8];
% titer_norm = titer/titer(1);

% [rho,pval] = corr([EnergySeq-EnergySeq(1)]',titer_norm','tail','left','type','spearman')


%%  LAIO et al. 1997 (deleted background)
vp1_mahoney_V160I = vp1_mahoney;
mutant_order{rev_translation_indices(160,true_indices)};
vp1_mahoney_V160I(160) = 'I';

vp1_mahoney_W170R = vp1_mahoney; 
% mutant_order{rev_translation_indices(170,true_indices)}
vp1_mahoney_W170R(170) = 'W';

vp1_mahoney_T177S = vp1_mahoney;
% mutant_order{rev_translation_indices(177,true_indices)}
vp1_mahoney_T177S(177) = 'T';

vp1_mahoney_V160I_W170R = vp1_mahoney;
vp1_mahoney_V160I_W170R([160 170]) = 'IR';

vp1_mahoney_V160I_T177S = vp1_mahoney;
vp1_mahoney_V160I_T177S([160 177]) = 'IK';

vp1_mahoney_W170R_T177S = vp1_mahoney;
vp1_mahoney_W170R_T177S([170 177]) = 'RS';

vp1_mahoney_V160I_W170R_T177S = vp1_mahoney;
vp1_mahoney_V160I_W170R_T177S([160 170 177]) = 'IRS';

% msaToTestComplete = [vp1_mahoney;vp1_mahoney_T99K]
msaToTestComplete = [vp1_mahoney;vp1_mahoney_V160I;vp1_mahoney_W170R;vp1_mahoney_T177S;vp1_mahoney_V160I_W170R;...
    vp1_mahoney_V160I_T177S;vp1_mahoney_W170R_T177S;vp1_mahoney_V160I_W170R_T177S];

% msaToTest = [Cseq;Cseq_99K;Cseq_99T;msaToTestComplete(:,true_indices)];
msaToTest = msaToTestComplete(:,true_indices);
% msaToTest(:,1:4) = '-';
% msaToTest = Cseq;
H_first4aaZero = H;

phi_cumulative(1) = phi_curr(1);
for kk = 2:length(mutant_order)
    phi_cumulative(kk) = sum(phi_curr(1:kk));
end
phi_cum = [0 phi_cumulative];

del_aa = 98:102;

H_first4aaZero(phi_cum(rev_translation_indices(del_aa(1),true_indices))+1:phi_cum(rev_translation_indices(del_aa(end)+1,true_indices)),:)=0;
H_first4aaZero(:,phi_cum(rev_translation_indices(del_aa(1),true_indices))+1:phi_cum(rev_translation_indices(del_aa(end)+1,true_indices)))=0;

numSeqMsaToTest = size(msaToTest, 1);
out_seq_ex = zeros(numSeqMsaToTest, total_length);
EnergySeq2 = zeros(1, numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    
    input_seq = msaToTest(i,:);

    input_parm = cell(1,2);
    input_parm{1} = protein_length_aa;
    input_parm{2} = input_seq;

    [out_seq_ex(i,:)] = convertAAseq2Bin(mutant_order,diff_amino_site_length_RawMSA,bin_matrix,input_parm);
    EnergySeq2(i) = calcSeqEnergy(out_seq_ex(i,:),H_first4aaZero);
    
    indices_diff = find(vp1_mahoney~=msaToTestComplete(i,:));
    no_conserved_sites_in_indices_diff = sum(ismember(indices_diff,site_freq_no_mutation));
    
    EnergySeq2(i) = EnergySeq2(i)- no_conserved_sites_in_indices_diff*min(diag(H));
    
end

titer2 = [8.5e8 2.3e8 4.3e7 8.6e7 2.6e7 2.8e7 9.3e7 2.6e7];
% titer2_norm = titer2/titer2(1);
% 
% [rho,pval] = corr([EnergySeq2-EnergySeq2(1)]',titer2_norm','tail','left','type','spearman')


%% Bouchard et al. 1995

nuc_index = [2585 2741 2749 2762 2775 2795 2879];
new_nuc = 'GGAUAAU';
for k = 1:length(nuc_index)
    [nuc_mahoney(k),aa_mahoney(k),aa_variant(k),aa_index(k)] = ...
        convert_nuc_index_to_aa_index_mahoney(nuc_index(k),new_nuc(k));
end


% aa_mahoney
% aa_variant
% aa_index
vp1_aa_index = aa_index-579;

load vp1_mahoney
vp1_mahoney_T36A = vp1_mahoney;
mutant_order{rev_translation_indices(36,true_indices)};
vp1_mahoney_T36A(36) = 'A';

vp1_mahoney_T88A = vp1_mahoney;
mutant_order{rev_translation_indices(88,true_indices)};
vp1_mahoney_T88A(88) = 'A';

vp1_mahoney_M90I = vp1_mahoney;%
mutant_order{rev_translation_indices(90,true_indices)};
vp1_mahoney_M90I(90) = 'L';

vp1_mahoney_P95S = vp1_mahoney;
mutant_order{rev_translation_indices(95,true_indices)};
vp1_mahoney_P95S(95) = 'S';

vp1_mahoney_T99K = vp1_mahoney;%
mutant_order{rev_translation_indices(99,true_indices)};
vp1_mahoney_T99K(99) = 'A';

vp1_mahoney_A106T = vp1_mahoney;
mutant_order{rev_translation_indices(106,true_indices)};
vp1_mahoney_A106T(106) = 'T';

vp1_mahoney_L134F = vp1_mahoney; %%%% problematic? %%%
% mutant_order{rev_translation_indices(134,true_indices)}
vp1_mahoney_L134F(134) = 'F';



% msaToTestComplete = [vp1_mahoney;vp1_mahoney_T99K]
msaToTestComplete = [vp1_mahoney;vp1_mahoney_T36A;vp1_mahoney_T88A;vp1_mahoney_M90I;vp1_mahoney_P95S;...
    vp1_mahoney_T99K;vp1_mahoney_A106T;vp1_mahoney_L134F];

% msaToTest = [Cseq;Cseq_99K;Cseq_99T;msaToTestComplete(:,true_indices)];
msaToTest = msaToTestComplete(:,true_indices);
% msaToTest = Cseq;

numSeqMsaToTest = size(msaToTest, 1);
out_seq_ex = zeros(numSeqMsaToTest, total_length);
EnergySeq3 = zeros(1, numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    
    input_seq = msaToTest(i,:);

    input_parm = cell(1,2);
    input_parm{1} = protein_length_aa;
    input_parm{2} = input_seq;

    [out_seq_ex(i,:)] = convertAAseq2Bin(mutant_order,diff_amino_site_length_RawMSA,bin_matrix,input_parm);
    EnergySeq3(i) = calcSeqEnergy(out_seq_ex(i,:),H);
    
    if i==8
        EnergySeq3(i)=EnergySeq3(i)-min(diag(H));
    end
end

% titer3 = [0.35 .37 0.25 0.32 .38 0.29 .29 .14]; %plaque size at 37C
titer3 = [0.35 .37 0.25 0.32 .38 0.29 .29 .14] .* [.9 .9 .5 .9 1 .8 1.2 .8]; %plaque size at 40C

% titer3_norm = titer3/titer3(1);
% [rho,pval] = corr([EnergySeq3-EnergySeq3(1)]',titer3_norm','tail','left','type','spearman')


%% Colston et al. 1994

vp1_mahoney_G225D = vp1_mahoney;
% mutant_order{rev_translation_indices(225,true_indices)}
vp1_mahoney_G225D(225) = 'D';

vp1_mahoney_D226G = vp1_mahoney;
% mutant_order{rev_translation_indices(226,true_indices)}
vp1_mahoney_D226G(226) = 'G';

vp1_mahoney_D226N = vp1_mahoney;
% mutant_order{rev_translation_indices(226,true_indices)}
vp1_mahoney_D226N(226) = 'N';

vp1_mahoney_L228F = vp1_mahoney;
% mutant_order{rev_translation_indices(228,true_indices)}
vp1_mahoney_L228F(228) = 'F';

vp1_mahoney_A231V = vp1_mahoney;
mutant_order{rev_translation_indices(231,true_indices)};
vp1_mahoney_A231V(231) = 'V';

vp1_mahoney_L234P = vp1_mahoney;
% mutant_order{rev_translation_indices(234,true_indices)}
vp1_mahoney_L234P(234) = 'P';

vp1_mahoney_D236G = vp1_mahoney;
% mutant_order{rev_translation_indices(236,true_indices)}
vp1_mahoney_D236G(236) = 'G';

vp1_mahoney_M132I = vp1_mahoney;
mutant_order{rev_translation_indices(132,true_indices)};
vp1_mahoney_M132I(132) = 'I';

vp1_mahoney_A241V = vp1_mahoney;
mutant_order{rev_translation_indices(241,true_indices)};
vp1_mahoney_A241V(241) = 'S';

vp1_mahoney_H265R = vp1_mahoney;
% mutant_order{rev_translation_indices(265,true_indices)}
vp1_mahoney_H265R(265) = 'Y';


msaToTestComplete = [vp1_mahoney;vp1_mahoney_G225D;vp1_mahoney_D226G;vp1_mahoney_D226N;...
    vp1_mahoney_L228F;vp1_mahoney_A231V;vp1_mahoney_L234P;vp1_mahoney_D236G;vp1_mahoney_M132I;...
    vp1_mahoney_A241V;vp1_mahoney_H265R];

% msaToTest = [Cseq;Cseq_99K;Cseq_99T;msaToTestComplete(:,true_indices)];
msaToTest = msaToTestComplete(:,true_indices);
% msaToTest = Cseq;

numSeqMsaToTest = size(msaToTest, 1);
out_seq_ex = zeros(numSeqMsaToTest, total_length);
EnergySeq5 = zeros(1, numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    
    input_seq = msaToTest(i,:);

    input_parm = cell(1,2);
    input_parm{1} = protein_length_aa;
    input_parm{2} = input_seq;

    [out_seq_ex(i,:)] = convertAAseq2Bin(mutant_order,diff_amino_site_length_RawMSA,bin_matrix,input_parm);
    EnergySeq5(i) = calcSeqEnergy(out_seq_ex(i,:),H);
    
    if ismember(i,[2:5 7 8 11])
        EnergySeq5(i) = EnergySeq5(i)-min(diag(H));
    end
end

EnergySeq5 = [EnergySeq5(1) mean(EnergySeq5(2:11))];

titer5 = [4.7e7 4.4e6];
% titer5_norm = titer5/titer5(1);
% [rho,pval] = corr([EnergySeq5-EnergySeq5(1)]',titer5_norm','tail','left','type','spearman')


%% Standardized data plot

Energy1_standardized = standardize_data(EnergySeq);
Energy2_standardized = standardize_data(EnergySeq2);
Energy3_standardized = standardize_data(EnergySeq3);
Energy4_standardized = standardize_data(EnergySeq_Shulman2014);
Energy5_standardized = standardize_data(EnergySeq5);


titer1_standardized = standardize_data(titer);
titer2_standardized = standardize_data(titer2);
titer3_standardized = standardize_data(titer3);
titer4_standardized = standardize_data(pd50_neurovirulence);
titer5_standardized = standardize_data(titer5);


%% Color Scheme definition
color_scheme = brewermap(9,'Set1');
color_scheme(10,:) = 115*ones(1,3)/255;


red = color_scheme(1,:);
blue = color_scheme(2,:);
green = color_scheme(3,:);
purple = color_scheme(4,:);
orange = color_scheme(5,:);
yellow = color_scheme(6,:);
brown = color_scheme(7,:);
pink = color_scheme(8,:);
gray = color_scheme(9,:);
darkgray = 115*ones(1,3)/255;
white = [1 1 1];

lightblue = [66,146,150]/255;
lightgreen = [161,217,155]/255;

%%


markersize = 7;
line_width = 0.5;

figure;
plot(Energy5_standardized,titer5_standardized,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',darkgray,'Linewidth',line_width)%,'MarkerFaceColor','r')
hold on 
plot(Energy1_standardized,titer1_standardized,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',pink,'Linewidth',line_width)
xlabel('Normalized energy')
ylabel('Normalized fitness')
% xlabel('$${\rm E}_{\rm mutant}-{\rm E}_{\rm Mahoney}$$','interpreter','latex')
% ylabel('$${\rm f}_{\rm mutant}/{\rm f}_{\rm Mahoney}$$','interpreter','latex')
% ylabel('$$\frac{f_{mutant}}{f_{Mahoney}}$$','interpreter','latex')
% plot(Energy5,titer5,'<','MarkerEdgeColor',purple,'MarkerSize',10,'Linewidth',1.5)
plot(Energy3_standardized,titer3_standardized,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',brown,'LineWidth',line_width)
plot(Energy2_standardized,titer2_standardized,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',orange,'Linewidth',line_width)%,'MarkerFaceColor','r')
plot(Energy4_standardized,titer4_standardized,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',purple,'Linewidth',line_width)%,'MarkerFaceColor','r')
% ylim([0 1.5])
% xlim([-2 2])
% xlabel('E-Eref')
% ylabel('f/fref');

% E = [Energy1_standardized Energy2_standardized Energy3_standardized Energy4_standardized Energy5_standardized]
% t = [titer1_standardized titer2_standardized titer3_standardized titer4_standardized titer5_standardized]

[rho,pval] = corr([Energy1_standardized Energy2_standardized Energy3_standardized Energy4_standardized Energy5_standardized]',...
    [titer1_standardized titer2_standardized titer3_standardized titer4_standardized titer5_standardized]','tail','left','type','spearman')
% title('$$\rho_s = -0.35$$','interpreter','latex','FontSize',12)
text(1,2.6,'$$\rho_s = -0.51$$','interpreter','latex','FontSize',12)
% % text(1.5,2.3,'$$p \;= 0.042$$','interpreter','latex','FontSize',12)
text(1,2.2,'$$p \;= 0.003$$','interpreter','latex','FontSize',12)

P = polyfit([Energy1_standardized Energy2_standardized Energy3_standardized Energy4_standardized Energy5_standardized]',...
    [titer1_standardized titer2_standardized titer3_standardized titer4_standardized titer5_standardized]',1);
% x = -10:30; %xaxis
x = -3:3; %xaxis
y = P(1)*x+P(2);
plot(x,y,'k--','LineWidth',1,'HandleVisibility','off')
h = legend('Colston1994','Colston1995','Bouchard1995','Liao1997','Shulman2014','Location','SouthWest');
rect = [0.17 0.2 0.1786 0.1393];
set(h, 'Position', rect)
legend boxoff

xlim([-3 3])
ylim([-3 3])

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'YTick'       , -3:1:3, ...
  'XTick'       , -3:1:3, ...
  'LineWidth'   , 1        );


