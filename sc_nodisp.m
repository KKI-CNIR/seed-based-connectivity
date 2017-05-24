function [cluster_labels] = sc_nodisp(A, sigma, noC)
%Spectral clustering using a sparse similarity matrix (t-nearest-neighbor).
%where
%   A - NxN sparse distance matrix (N = # of voxels in the ROI)
%   sigma - sigma value used in computing affinity; if 0, apply self-tunning technique
%   noC - number of clusters
%   cluster_labels - N-by-1 vector containing cluster labels
%
%   Based on code from Wen-Yen Chen (wychen@alumni.cs.ucsb.edu) & Chih-Jen
%   Lin (cjlin@csie.ntu.edu.tw)
%   "Parallel Spectral Clustering in Distributed Systems"
%   IEEE Transactions on Pattern Analysis and Machine Learning
%
%   Modified MB Nebel 

n = size(A, 1);

if (sigma == 0) % Selftuning spectral clustering
    % Find the # of nonzero elements for each column
    col_count = sum(A~=0, 1)';
    col_sum = sum(A, 1)';
    col_mean = col_sum ./ col_count;
    [x y val] = find(A);
    A = sparse(x, y, -val.*val./col_mean(x)./col_mean(y)./2);
    clear col_count col_sum col_mean x y val;
else % Fixed-sigma spectral clustering
    A = A.*A;
    A = -A/(2*sigma*sigma);
end

% Do exp function sequentially because of memory limitation
num = 2000;
num_iter = ceil(n/num);
S = sparse([]);
for i = 1:num_iter
    start_index = 1 + (i-1)*num;
    end_index = min(i*num, n);
    S1 = spfun(@exp, A(:,start_index:end_index)); % sparse exponential func
    S = [S S1];
    clear S1;
end
clear A;

% Calculate the Laplacian, L = D^(-1/2) * S * D^(-1/2)
DD = sum(S, 2) + eps;
DD = sqrt(1./DD); % D^(-1/2)
DD = spdiags(DD, 0, n, n);
L = DD * S * DD;

clear S;

% Do eigen decomposition, 
%   if L = D^(-1/2) * S * D(-1/2)    : set 'LM' (Largest Magnitude), or
%   if L = I - D^(-1/2) * S * D(-1/2): set 'SM' (Smallest Magnitude).
%   if L = D^(-1)*S or I - D^(-1)*S, don't have to do second normalization
%   step (random walk)

OPTS.issym = 1;
OPTS.disp = 0;

[V, val] = eigs(L, noC, 'LM', OPTS);

% v2 = D*V(:, noC-1);
% [Y, A_order] = sort(v2);

clear v2 Y DD;

% Normalize each row to be of unit length and run k-means
sq_sum = sqrt(sum(V.*V, 2)) + eps;
U = V ./ repmat(sq_sum, 1, noC);
clear sq_sum V;
cluster_labels = k_means_nodisp(U, [], noC);

% % Alternative version of k-means (completes 10 iterations and returns best)
% % from Zelnik-Manor (2005)
% [cluster_indices bestD] = kmeans2(U, noC);
% cluster_labels = zeros(n, 1);
% for ic = 1:noC
%     cluster_labels(cluster_indices{ic}) = ic;
% end

