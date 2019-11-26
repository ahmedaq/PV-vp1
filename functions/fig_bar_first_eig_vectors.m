function [eigvect,lambda]=fig_bar_first_eig_vectors(C, textdata,num_eig_vectors_to_plot,flag_plot)

if nargin < 4
    flag_plot = 1;
end

[eigvect_unsorted,lambda_unsorted]=eig(C);
[lambda,lambda_order]=sort(diag(lambda_unsorted),'descend');
eigvect=eigvect_unsorted(:,lambda_order);

if flag_plot == 1
    figure;
    for kkk = 1:num_eig_vectors_to_plot
        subplot(num_eig_vectors_to_plot,1,kkk);
        bar(eigvect(:,kkk));
        xlabel(sprintf('Eigenvector %d of %s',kkk,textdata))
        ylim([-0.5 0.5])
    end
end