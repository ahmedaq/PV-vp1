% clear all;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to obtain sequences poliovirus protein
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ahmed Abdul Quadeer
% Last revised date: 6th February, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [msa_VP1, header_VP1] = stats_seqs_PV1(inputfile)

run startup.m

%% Loading sequence data

% [header,seqs]=fastaread(inputfile);
% 
% for k = 1:size(seqs,2)
%    lengths_seqs(k) = length(seqs{k});
% end
% 
% hist(lengths_seqs,size(seqs,2))
% length_max = max(lengths_seqs)

% After using MAFFT directly on the fasta file

[header,seqs]=fastaread('PV1_aligned_MAFFT.fasta');

%% Removing sequences with gaps>3000

for kk = 1:length(header)
    msa(kk,:) = seqs{kk};
end

no_of_gaps_per_site = plot_no_of_gaps_per_site(msa,'raw MSA',0);

sites_3000gaps = find(no_of_gaps_per_site>3000);
msa(:,sites_3000gaps) = [];

plot_no_of_gaps_per_site(msa,'MSA after removal of gaps/site >3000',0);

[Nseq,Npos] = size(msa);

%% Finding Mahoney strain in msa and comparing it with the actual one
for k = 1:Nseq
    if ~isempty(strfind(header{k},'CAA24461'))
        m = k;
    end
end


[header_Mahoney,seq_Mahoney] = fastaread('Mahoney_strain_Human_poliovirus1.fasta');
[score_global_aa, globalAlignment_aa] = nwalign(msa(m,:),seq_Mahoney);
% showalignment(globalAlignment_aa)

%% Separating capsid proteins msa (based on ACCESSION P03300)

msa_capsid = msa(:,2:881);
msa_VP0 = msa(:,2:341);
msa_VP4 = msa(:,2:69);
msa_VP2 = msa(:,70:341);
msa_VP3 = msa(:,342:579);
msa_VP1 = msa(:,580:881);

msa_2A = msa(:,882:1030);
msa_2B = msa(:,1031:1128);
msa_2C = msa(:,1129:1456); %actual was starting from 1128
msa_3AB = msa(:,1457:1565);
msa_3A = msa(:,1457:1543);
msa_3B = msa(:,1544:1565); %inferred
msa_3CD = msa(:,1566:2209);
msa_3C = msa(:,1566:1747);
msa_3D =  msa(:,1748:2209); %inferred

plot_no_of_gaps_per_site(msa_capsid,'capsid',0);
plot_no_of_gaps_per_site(msa_VP0,'VP0',0);
plot_no_of_gaps_per_site(msa_VP4,'VP4',0);
plot_no_of_gaps_per_site(msa_VP2,'VP2',0);
plot_no_of_gaps_per_site(msa_VP3,'VP3',0);
plot_no_of_gaps_per_site(msa_VP1,'VP1',0);

plot_no_of_gaps_per_site(msa_2A,'2A',0);
plot_no_of_gaps_per_site(msa_2B,'2B',0);
plot_no_of_gaps_per_site(msa_2C,'2C',0);
plot_no_of_gaps_per_site(msa_3AB,'3AB',0);
plot_no_of_gaps_per_site(msa_3A,'3A',0);
plot_no_of_gaps_per_site(msa_3B,'3B',0);
plot_no_of_gaps_per_site(msa_3CD,'3CD',0);
plot_no_of_gaps_per_site(msa_3C,'3C',0);
plot_no_of_gaps_per_site(msa_3D,'3D',0);

%% Removing strains with large number of gaps 

percent_gaps_thresh_in_a_seq = 0.2; %0.2 --> Nseq Vp1 = 2261

[msa_capsid_v2,indx_seqs_removed_capsid] = removing_seqs_from_msa_with_many_gaps(msa_capsid,percent_gaps_thresh_in_a_seq,'Capsid');
[msa_VP1_v2,indx_seqs_removed_VP1] = removing_seqs_from_msa_with_many_gaps(msa_VP1,percent_gaps_thresh_in_a_seq,'VP1');
[msa_VP0_v2,indx_seqs_removed_VP0] = removing_seqs_from_msa_with_many_gaps(msa_VP0,percent_gaps_thresh_in_a_seq,'VP0');
[msa_VP4_v2,indx_seqs_removed_VP4] = removing_seqs_from_msa_with_many_gaps(msa_VP4,percent_gaps_thresh_in_a_seq,'VP4');
[msa_VP2_v2,indx_seqs_removed_VP2] = removing_seqs_from_msa_with_many_gaps(msa_VP2,percent_gaps_thresh_in_a_seq,'VP2');
[msa_VP3_v2,indx_seqs_removed_VP3] = removing_seqs_from_msa_with_many_gaps(msa_VP3,percent_gaps_thresh_in_a_seq,'VP3');

[msa_2A_v2,indx_seqs_removed_2A] = removing_seqs_from_msa_with_many_gaps(msa_2A,percent_gaps_thresh_in_a_seq,'2A');
[msa_2B_v2,indx_seqs_removed_2B] = removing_seqs_from_msa_with_many_gaps(msa_2B,percent_gaps_thresh_in_a_seq,'2B');
[msa_2C_v2,indx_seqs_removed_2C] = removing_seqs_from_msa_with_many_gaps(msa_2C,percent_gaps_thresh_in_a_seq,'2C');
[msa_3AB_v2,indx_seqs_removed_3AB] = removing_seqs_from_msa_with_many_gaps(msa_3AB,percent_gaps_thresh_in_a_seq,'3AB');
[msa_3A_v2,indx_seqs_removed_3A] = removing_seqs_from_msa_with_many_gaps(msa_3A,percent_gaps_thresh_in_a_seq,'3A');
[msa_3B_v2,indx_seqs_removed_3B] = removing_seqs_from_msa_with_many_gaps(msa_3B,percent_gaps_thresh_in_a_seq,'3B');
[msa_3CD_v2,indx_seqs_removed_3CD] = removing_seqs_from_msa_with_many_gaps(msa_3CD,percent_gaps_thresh_in_a_seq,'3CD');
[msa_3C_v2,indx_seqs_removed_3C] = removing_seqs_from_msa_with_many_gaps(msa_3C,percent_gaps_thresh_in_a_seq,'3C');
[msa_3D_v2,indx_seqs_removed_3D] = removing_seqs_from_msa_with_many_gaps(msa_3D,percent_gaps_thresh_in_a_seq,'3D');

%% Displaying MSA

indx_U = find(msa_VP1_v2=='U');
indx_J = find(msa_VP1_v2=='J');
indx_X = find(msa_VP1_v2=='X');
indx_Z = find(msa_VP1_v2=='Z');
indx_B = find(msa_VP1_v2=='B');
indx_O = find(msa_VP1_v2=='O');

msa_VP1_v3 = msa_VP1_v2;
msa_VP1_v3([indx_U;indx_J;indx_X;indx_Z;indx_O;indx_B])='B';

% seqalignviewer(msa_VP1_v3)

%% Number of sequences and positions in each protein

[Nseq_capsid,Npos_capsid] = size(msa_capsid_v2);
[Nseq_VP1,Npos_VP1] = size(msa_VP1_v2);
[Nseq_VP0,Npos_VP0] = size(msa_VP0_v2);
[Nseq_VP2,Npos_VP2] = size(msa_VP2_v2);
[Nseq_VP3,Npos_VP3] = size(msa_VP3_v2);
[Nseq_VP4,Npos_VP4] = size(msa_VP4_v2);

[Nseq_2A,Npos_2A] = size(msa_2A_v2);
[Nseq_2B,Npos_2B] = size(msa_2B_v2);
[Nseq_2C,Npos_2C] = size(msa_2C_v2);
[Nseq_3AB,Npos_3AB] = size(msa_3AB_v2);
[Nseq_3A,Npos_3A] = size(msa_3A_v2);
[Nseq_3B,Npos_3B] = size(msa_3B_v2);
[Nseq_3CD,Npos_3CD] = size(msa_3CD_v2);
[Nseq_3C,Npos_3C] = size(msa_3C_v2);
[Nseq_3D,Npos_3D] = size(msa_3D_v2);

%% Plots

figure('units','normalized','outerposition',[0 0 .7 .7])

subplot(3,1,1)
bar(1:10,[Nseq_VP4,Nseq_VP2,Nseq_VP3,Nseq_VP1,Nseq_2A,...
    Nseq_2B,Nseq_2C,Nseq_3AB,Nseq_3C,Nseq_3D],0.3,...
    'FaceColor',darkgray,'EdgeColor',darkgray)
set(gca, 'XTick',1:10, 'XTickLabel',{'vp4' 'vp2' 'vp3' 'vp1' ...
    '2A' '2B' '2C' '3AB' '3C' '3D'})
xlabel('Protein (PV1)')
ylabel('Number of sequences')
grid off
axis([0 11 0 2500])
post_proc_fig

subplot(3,1,2)
bar(1:10,[Npos_VP4,Npos_VP2,Npos_VP3,Npos_VP1,Npos_2A,...
    Npos_2B,Npos_2C,Npos_3AB,Npos_3C,Npos_3D],0.3,...
    'FaceColor',darkgray,'EdgeColor',darkgray)
set(gca, 'XTick',1:10, 'XTickLabel',{'vp4' 'vp2' 'vp3' 'vp1' ...
    '2A' '2B' '2C' '3AB' '3C' '3D'})
xlabel('Protein (PV1)')
ylabel('Number of residues')
grid off
axis([0 11 0 500])
post_proc_fig

subplot(3,1,3)
bar(1:10,[Nseq_VP4/Npos_VP4,...
    Nseq_VP2/Npos_VP2,Nseq_VP3/Npos_VP3,Nseq_VP1/Npos_VP1,Nseq_2A/Npos_2A,...
    Nseq_2B/Npos_2B,Nseq_2C/Npos_2C,Nseq_3AB/Npos_3AB,...
    Nseq_3C/Npos_3C,Nseq_3D/Npos_3D],0.3,...
    'FaceColor',darkgray,'EdgeColor',darkgray)
set(gca, 'XTick',1:10, 'XTickLabel',{'vp4' 'vp2' 'vp3' 'vp1' ...
    '2A' '2B' '2C' '3AB' '3C' '3D'})
xlabel('Protein (PV1)')
ylabel('No. of sequences / No. of residues')
grid off
axis([0 11 0 8])
post_proc_fig

%% saving msa and header in .mat format

msa_VP1 = msa_VP1_v3;
header_VP1 = header;
header_VP1(indx_seqs_removed_VP1) = [];

