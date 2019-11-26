
function [msa_vp1_wt,sec1] = analysis_similarity_matrix(msa_vp1)

run startup.m
msa = msa_vp1;

%% Incorporating Sabin sequence in MSA
[header_Sabin1_aa,seq_Sabin1_aa] = fastaread('Sabin1_aa.fasta');
seq_Sabin1_aa_VP1 = seq_Sabin1_aa(580:881);

msa = [msa;seq_Sabin1_aa_VP1];

[ns,ls] = size(msa);

%% Similarity matrix and its EigenValue Decomposition (EVD)

Sim_mtrx = sim_seq(msa);

[eigvect,lambda] = fig_bar_first_eig_vectors(Sim_mtrx, 'Similarity Matrix', 3, 0);

eig_vec1 = eigvect(:,1);
eig_vec2 = eigvect(:,2);

% figure;
% scatter(eigvect(:,1),eigvect(:,2),'bo');
% % title('Based on similarity matrix')
% xlabel('Eigenvector 1')
% ylabel('Eigenvector 2')

%% 
% sec1 = find(eig_vec2>=-0.01);
sec1 = find(eig_vec2>=-0.01 & eig_vec1>=.0195);
sec2 = find(eig_vec2<-0.01);

figure;
plot(eig_vec1(sec1),eig_vec2(sec1),'o','MarkerFaceColor',green,'MarkerEdgeColor','w','MarkerSize',6);
% title('Based on similarity matrix')
hold on;
plot(eig_vec1(sec2),eig_vec2(sec2),'o','MarkerFaceColor',orange,'MarkerEdgeColor','w','MarkerSize',6);

plot(eig_vec1(ns),eig_vec2(ns),'kp','LineWidth',1,'MarkerFaceColor','k','MarkerSize',10); %sabin strain

xlabel('Principal component 1')
ylabel('Principal component 2')
axis([0.019 0.022 -0.04 0.04])
% grid on
hold on;
legend('Wild-type','VDPV','Sabin strain','Location','NorthWest')
legend boxoff
post_proc_fig

%% Removing extra sequences associated with vaccines based on journals

load seqs_vaccine_derived_vp1

sec1(seqs_vaccine_derived_vp1) = [];
msa_vp1_wt = msa_vp1(sec1,:);