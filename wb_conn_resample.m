function [conn_map_all, comptime] = wb_conn_resample(WB, tsize, roi_tc)
%Function to compute whole brain correlation maps for specified ROI
%time series
%Usage
%   wb_conn_resample(WB, tsize, roi_tc)
%   where
%       WB - is the whole brain data set reshaped to 2D (timepoints by
%       voxels)
%       tsize - is the number of timepoints in WB
%       roi_tc - is the timecourse(s) for comparision (timepoints by ROIs)
%   output
%       conn_map_all - 2D connectivity maps (seeds X voxels) NOTE: can be
%       huge depending on how many "seeds" are used so this does not get
%       written out
%
% MBN Jan 12, 2011 - Modified the way that the roi_mask is loaded,
% incorporated jmuschelli's corrit.m & changed the output filename to
% reflect the # of TRs used
% MBN Feb 17, 2011 - roi_filenm_mask could contain mean tc for the roi or
% tc for each timepoint
% MBN March 3, 2011 - uses spm_read_vols instead of avw_read
% MBN July 21, 2011 incorporated jmuschelli's changes to correlation calculation
% which speeds up computation time so don't need to save conn_map_all.mat
% FC map(s) returned as output
% MBN July 26, 2011 - pass roi timecourses as an input variable instead of
% loading the file
% MBN March 13, 2012 - pass whole brain 4D data and number of timepoints as
% input variables

fc=tic;

%% Compute correlations (does NOT write files)

SDY = nanstd(WB);
means = nanmean(WB);

% allmeans = repmat(means, [tsize(4) 1]);
% allSDs = repmat(SDY, [tsize(4) 1]);
cY = (WB - repmat(means, [tsize 1]))./repmat(SDY, [tsize 1]);

clear WB

ROI_TC = roi_tc;
ROI_MEANS = repmat(mean(ROI_TC), [tsize 1]);
ROI_SD = repmat(std(ROI_TC), [tsize 1]);

ROI_TC = ((ROI_TC - ROI_MEANS)./ROI_SD/(tsize-1))';

clear ROI_MEANS ROI_SD
conn_map_all = ROI_TC * cY;

conn_map_all(isnan(conn_map_all)) = 0;
conn_map_all(isinf(conn_map_all)) = 0;

comptime = toc(fc);
disp(strcat(num2str(toc(fc)), ' s to calculate FC map'))

