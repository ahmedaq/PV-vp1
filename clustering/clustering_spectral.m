
function b = clustering_spectral(input_data_file,data_peaks_rank)

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

% Parameter selection
data_standardization = 1;
clusters = 1:40;
normalization = 1;
simGraph_type = 3; %1 = Knn, 2 = epsilon, 3 = full

if simGraph_type == 1
    Type_simGraph = 1;  % 'Type' - Type if kNN Graph, 1 - Normal, 2 - Mutual
    sigma = 1; % 'sigma' - Parameter for Gaussian similarity function. Set
    %  this to 0 for an unweighted graph. Default is 1.
elseif simGraph_type == 2
    epsilon = 0.2;
elseif simGraph_type == 3
    sigma = 1; % 'sigma' - Parameter for Gaussian similarity function
end

Type_sc = 1;
%   'Type' - Defines the type of spectral clustering algorithm
%            that should be used. Choices are:
%      1 - Unnormalized
%      2 - Normalized according to Shi and Malik (2000)
%      3 - Normalized according to Jordan and Weiss (2002)
no_of_eigenvalues = 5;

%%

for kkk = 1:length(clusters)
    
    no_of_clusters = clusters(kkk);
    
    A = msa_bin; %Data matrix
    [n,p] = size(A);
    
    if data_standardization == 1
        for kk = 1:p
            A_normalized(:,kk) = (A(:,kk)-mean(A(:,kk)))/sqrt(var(A(:,kk)));
        end
        B = A_normalized;
    end
    
    %Normalizing data
    if normalization == 1
        A = normalizeData(A.');
    else
        A = A.';
    end
    
    
    %Constructing simgraph
    if simGraph_type == 1
        W = SimGraph_NearestNeighbors(A, 50, Type_simGraph, sigma);
    elseif simGraph_type == 2
        W = SimGraph_Epsilon(A, epsilon);
    elseif simGraph_type == 3
        W = SimGraph_Full(A, sigma);
    end
    
    %Spectral Clustering
    [C, ~, U, ~] = SpectralClustering_smallevs(W, no_of_eigenvalues, no_of_clusters, Type_sc);
    
    C = full(C);
    
    for kk = 1:size(C,2)
        labels(C(:,kk)==1,kkk) = kk;
    end
    
end

%% Estimating the number of optimal clusters using Silhouette

[~, ~, ~, optimalSil] = evaClusters(U,labels);

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