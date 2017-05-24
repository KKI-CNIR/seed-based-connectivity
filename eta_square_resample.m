function [eta_sq, elpsd_time] = eta_square_resample(seeds, write_flag, outfile)
%Function to compute eta^2 matrix for specified 
%
%Usage
%   eta_square_resample(seeds, write_flag, outfile)
%   where
%       seeds - matlab variable (# voxels in ROI X # of voxels in WB)
%       write_flg - if =1, eta_sq will be saved to a .mat
%       outfile - is the path/name of the .mat to save
%       eta_sq - square matrix based on formula in Cohen et
%       al. 2008 NeuroImage 41(1):45-57.

% MBN July 21, 2011 - conn_map_all passed as input argument instead of
% being loaded; removed zsq_flag.
% MBN February 2, 2012 - added write flag; nanmean

et = tic;

%% Making the eta2 matrix after masking
[srows scols] = size(seeds);

colmeans = nanmean(seeds);

%(sum_i (a_i + bi)/2)n = 1/2*(sum_i ai/n sum_i bi/n) 
% 1/2 abar + 1/2 bbar
MEANS = (repmat(colmeans, [scols 1])' + repmat(colmeans, [scols 1]));

%ev_ss = sum a_i^2 + sum b_i^2
ev_ss = (repmat(sum(seeds.^2), [scols 1])' + repmat(sum(seeds.^2), [scols 1]));

% ev_ss = sum a_i^2 -2Msum a_i +M^2  +  sum b_i^2 - 2Msum b_i + M^2
% sum (a_i-M)^2 + (b_i -M)^2
denom = ev_ss - 1/2*srows*MEANS.^2;

aibi = seeds'*seeds;
%aibi = sum a_i b_i = 2 sum a_i bi /2

numer = ev_ss./2 - aibi;
%numer = [sum a_i^2 -  2sum a_i bi + sum b_i^2)/2 
%= sum a_i^2 -sum aibi + sum b_i ^2
% = sum (a_i-bi)^2 /2;colmeans = mean(sees_inbrain);

eta_sq = 1 - numer./denom;
% This is a full matrix, not just upper right triangular
% could use eta_sq = triu(eta_sq); if want just upper right

%spectral clustering algorithm requires non-negative values
eta_sq(isnan(eta_sq)) = 0;

if(write_flag)
save(outfile, 'eta_sq');
end

elpsd_time = toc(et);
disp(strcat(num2str(toc(et)), ' s to calculate eta^2'))

