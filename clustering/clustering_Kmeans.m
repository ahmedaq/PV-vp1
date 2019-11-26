function b = clustering_Kmeans(input_data_file,data_peaks_rank)

load(data_peaks_rank)
load(input_data_file)

[ns,ls] = size(msa_vp1_wt);

% Binarization

Cseq = seqconsensus(msa_vp1_wt);
Cseq_mtrx = repmat(Cseq,ns,1);
msa_bin = double(msa_vp1_wt==Cseq_mtrx); %representing mutation by 0 and no mutation by 1

% To remove 100% conserved sites
freq_single_mutation = mean(msa_bin);
site_freq_no_mutation = find(freq_single_mutation==1);
msa_bin(:,site_freq_no_mutation)=[];

A = msa_bin;
[n,p] = size(A)
data_standardization = 1;
no_of_clusters = 40;

if data_standardization == 1
    for kk = 1:p
        A_normalized(:,kk) = (A(:,kk)-mean(A(:,kk)))/sqrt(var(A(:,kk)));
    end
    A = A_normalized;
end


%%
for kk = 1:no_of_clusters
    labels(:,kk) = kmeans(A,kk);
end

%% Estimating the number of optimal clusters using Silhouette

[~, ~, ~, optimalSil] = evaClusters(A,labels);

%% Contrasting obtained clusters against peaks

% opt_clusters_true = optimalSil; %estimated
opt_clusters_true = 25;% %based on peaks
labels_optimal = labels(:,opt_clusters_true);

%
b=[];
for kk = 1:25%length(indices_seqs_that_converged_to_peak_seq)
    for mm = 1:opt_clusters_true %clusters formed
        a = find(labels_optimal==mm);
        b(kk,mm) = sum(ismember(a,indices_seqs_that_converged_to_peak_seq_vp1{kk}))...
            /length(indices_seqs_that_converged_to_peak_seq_vp1{kk})*100;
    end
end

%row i of matrix "b" shows percetage of peak i seqs in cluster 1:opt_clusters_true

[max_b,indx_b] = max(b.'); %maximum percentage of peak i in any cluster
[[1:25].' indx_b.' max_b.' ]
length(unique(indx_b));

%% Data for plotting in Python

bb = b(1:10,:);

num_peaks = 10; %top 10

peaks = 1:num_peaks;
peaks_vec = repmat(peaks.',size(bb,2),1);

clusters = 1:size(bb,2);
clusters_vec = [];
for kk = 1:length(clusters)
    clusters_vec = [clusters_vec;repmat(clusters(kk),num_peaks,1)];
end

bbb = [peaks_vec clusters_vec bb(:)];

csvwrite('clustering_kmeans_vec25.csv',bbb)