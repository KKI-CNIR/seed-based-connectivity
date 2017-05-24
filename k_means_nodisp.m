function cluster_labels = k_means_nodisp(data, centers, nClusts)
%Euclidean k-means clustering algorithm.
%
%   Input    : data           : N-by-D data matrix, where N is the number of data,
%                               D is the number of dimensions
%              centers        : K-by-D matrix, where K is nClusts, or
%                               'random', random initialization, or
%                               [], empty matrix, orthogonal initialization
%              nClusts   : Number of clusters
%
%   Output   : cluster_labels : N-by-1 vector of cluster assignment
%
%   Author   : Dimitrios Zeimpekis, Efstratios Gallopoulos, 2006.
%              http://scgroup.hpclab.ceid.upatras.gr/scgroup/Projects/TMG/
%
%   Modified : Wen-Yen Chen (wychen@alumni.cs.ucsb.edu)
%			   Chih-Jen Lin (cjlin@csie.ntu.edu.tw)

%
% Parameter setting
%
iter = 0;
qold = inf;
threshold = 0.0001;

%
% Check if with initial centers
%
if strcmp(centers, 'random')
  disp('Random initialization...');
  centers = random_init(data, nClusts);
elseif isempty(centers)
%   disp('Orthogonal initialization...');
  centers = orth_init(data, nClusts);
end

%
% Double type is required for sparse matrix multiply
%
data = double(data);
centers = double(centers);

%
% Calculate the distance (square) between data and centers
%
n = size(data, 1);
x = sum(data.*data, 2)';
X = x(ones(nClusts, 1), :);
y = sum(centers.*centers, 2);
Y = y(:, ones(n, 1));
P = X + Y - 2*centers*data';

%
% Main program
%
while 1
  iter = iter + 1;

  % Find the closest cluster for each data point
  [val, ind] = min(P, [], 1);
  % Sum up data points within each cluster
  P = sparse(ind, 1:n, 1, nClusts, n);
  centers = P*data;
  % Size of each cluster, for cluster whose size is 0 we keep it empty
  cluster_size = P*ones(n, 1);
  % For empty clusters, initialize again
  zero_cluster = find(cluster_size==0);
  if length(zero_cluster) > 0
    disp('Zero centroid. Initialize again...');
    centers(zero_cluster, :)= random_init(data, length(zero_cluster));
    cluster_size(zero_cluster) = 1;
  end
  % Update centers
  centers = spdiags(1./cluster_size, 0, nClusts, nClusts)*centers;

  % Update distance (square) to new centers
  y = sum(centers.*centers, 2);
  Y = y(:, ones(n, 1));
  P = X + Y - 2*centers*data';

  % Calculate objective function value
  qnew = sum(sum(sparse(ind, 1:n, 1, size(P, 1), size(P, 2)).*P));
%   mesg = sprintf('Iteration %d:\n\tQold=%g\t\tQnew=%g', iter, full(qold), full(qnew));
%   disp(mesg);

  % Check if objective function value is less than/equal to threshold
  if threshold >= abs((qnew-qold)/qold)
%     mesg = sprintf('\nkmeans converged!');
%     disp(mesg);
    break;
  end
  qold = qnew;
end

cluster_labels = ind';


%-----------------------------------------------------------------------------
function init_centers = random_init(data, nClusts)
%RANDOM_INIT Initialize centroids choosing nClusts rows of data at random
%
%   Input : data         : N-by-D data matrix, where N is the number of data,
%                          D is the number of dimensions
%           nClusts : Number of clusters
%
%   Output: init_centers : K-by-D matrix, where K is nClusts
rand('twister', sum(100*clock));
init_centers = data(ceil(size(data, 1)*rand(1, nClusts)), :);

function init_centers = orth_init(data, nClusts)
%ORTH_INIT Initialize orthogonal centers for k-means clustering algorithm.
%
%   Input : data         : N-by-D data matrix, where N is the number of data,
%                          D is the number of dimensions
%           nClusts : Number of clusters
%
%   Output: init_centers : K-by-D matrix, where K is nClusts

%
% Find the nClusts centers which are orthogonal to each other
%
Uniq = unique(data, 'rows'); % Avoid duplicate centers
num = size(Uniq, 1);
first = ceil(rand(1)*num); % Randomly select the first center
init_centers = zeros(nClusts, size(data, 2)); % Storage for centers
init_centers(1, :) = Uniq(first, :);
Uniq(first, :) = [];
c = zeros(num-1, 1); % Accumalated orthogonal values to existing centers for non-centers
% Find the rest nClusts-1 centers
for j = 2:nClusts
  c = c + abs(Uniq*init_centers(j-1, :)');
  [minimum, i] = min(c); % Select the most orthogonal one as next center
  init_centers(j, :) = Uniq(i, :);
  Uniq(i, :) = [];
  c(i) = [];
end
clear c Uniq;
