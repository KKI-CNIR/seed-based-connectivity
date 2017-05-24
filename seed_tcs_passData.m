function [tcs, fdata] = seed_tcs_passData(fdata, tps, roi_mask, outdir)
%Function to extract timecourses from seed regions 
%for one participant
%Usage
%   [tcs, Y] = seed_tcs_passData(fdata, tps, roi, out_dir)
%   fdata - voxels x time matrix containing preprocessed brain data
%   roi  - vector; if roi is binary, a separate
%   timecourse for each voxel within the roi will be extracted; otherwise
%   average timecourse for each region within the roi image will be
%   extracted
%   outdir - directory where the .csv containing the extracted timecourses will be saved
%   tcs - regions x timepoints; saves either the time series for each voxel
%   in the roi (if the roi image is binary) or the average time series for every
%   labeled region within the roi image
%   Y - the 4D dataset saved in memory (rows are voxels; columns are time
%   points); useful if calculating whole brain connectivity maps as next
%   step
%
% MBN Jan 28, 2015 - converted seed_tcs_gz.m so that data (voxels X
% timepoints) is passed as input instead of being loaded

%% Read ROI mask files
rois = nonzeros(unique(roi_mask));
nrois = length(rois);

% if roi file is binary, extract one timecourse for every voxel with the
% roi 
if nrois<2
    %%%MAKE SURE THIS STILL WORKS
    vox.roi_indx = find(roi_mask);
    %Initialize tcs - rows are voxels and columns are timepoints
    tcs = zeros(length(vox.roi_indx), size(fdata, 2));
    tcs = fdata(squeeze(vox.roi_indx), :);

%else calculate average timecourse for every parcel within roi file
else
    roi_reshape = reshape(roi_mask, [numel(roi_mask), 1]);
    % distribution of parcel membership
    roi_seeds = double(repmat(roi_reshape,1,nrois) == repmat(rois', numel(roi_mask),1));
    %multiply data by distrib, divide by column sums to get average &
    %transpose so regions x timepoints
    tcs = (fdata'*roi_seeds)./ repmat(sum(roi_seeds, 1), size(fdata, 2), 1)';
end

% if outdir defined, save timecourses of voxels w/in roi to .csv
if(exist('outdir','var'))
    csvwrite(fullfile(outdir, strcat('roi_voxels_', num2str(size(fdata, 2)), '.csv')), tcs);
%     save(fullfile(outdir, strcat('roi_voxels_', num2str(size(fdata, 2)), '.mat')), 'vox');
end
