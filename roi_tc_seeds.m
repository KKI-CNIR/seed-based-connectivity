function mn_roi = roi_tc_seeds(data_dir,roi_dir, ...
    data_filenm_mask, roi_filenm_mask, brain_mask_file)
%Function to create mean ROI timecourse given an ROI (mask/seed) and the
%dataset
%Usage
%   mn_roi_tc = roi_mn_tc(data_dir,roi_dir,data_filenm_mask, ...
%       roi_filename_mask,brain_mask_file)
%   data_dir - directory containing the preprocessed data files (3D nifti)
%   roi_dir  - directory containing the ROI masks
%   data_filenm_mask - string containing the common substring for the
%       preprocessed files (typically 'fswa' or 'swa')
%   roi_filenm_mask - string containing the common substring for the ROI
%       files (typically {'RSN','ACC'})
%   brain_mask_file - string with brain mask filename

% by Suresh E Joel - modified July, 2009; modified Sep 16,2009

if(nargin<4),
    error('Not enough arguements');
end;
% 
% data_dir = '/mnt/bisivity/nb/KKI/1018959';
% roi_dir = '/mnt/bisivity/nb/rois';
% data_filenm_mask = '*rest_1';
% roi_filenm_mask = '264seeds_mask_4mm';

%% Get orientation of data file (to match all masks to the same orientation)
ffiles=fullfile(data_dir,[data_filenm_mask,'.nii']);
V = spm_vol(ffiles);
% so = get_orient(V(1,:));
Y = spm_read_vols(V);
Y = reshape(Y, [numel(Y(:, :, :, 1)), numel(Y(1, 1, 1, :))]);
clear V;

%% Read brain mask is specified
% if(nargin>4)
%     %     brain_mask_file=fullfile(data_dir,[brain_mask_file,'.img']);
%     brain_mask_file=fullfile(data_dir,[brain_mask_file,'.nii']);
%     if(~strcmp(get_orient(brain_mask_file),so)),
%         change_orient(brain_mask_file,so);
%     end;
%     V=spm_vol(brain_mask_file);
%     bM=spm_read_vols(V);
% end;

%% Read the ROI mask files
hw=waitbar(0,'Reading Mask Files');
% if(ischar(roi_filenm_mask)), roi_filenm_mask{1}=roi_filenm_mask; end;
files=fullfile(roi_dir,[roi_filenm_mask,'.nii']);
%sM=zeros(V.dim(1),V.dim(2),V.dim(3),length(files));
% for i_roi=1:length(files),
% P{n}=[fullfile(roi_dir,files.name),',1']; %#ok
% if(~strcmp(get_orient(P{n}),so)),
%     change_orient(P{n},so);
% end;
[path mn_roi.name ext] = fileparts(files);
%     n=n+1;
% end;
% clear files;
% P=strvcat(P); %#ok

%warning off;%#ok
V=spm_vol(files);
sM=spm_read_vols(V);
seeds = unique(nonzeros(sM));
waitbar(1,hw);

%% Compute mean of the seed ROI region
waitbar(0,hw,'Reading & computing ROI mean timecourse');
clear files V P;
% sM(find(sM == 435)) = 265;

seedMat = false(size(Y, 1), length(seeds));
for iroi = 1:length(seeds)
   seedMat(:, iroi) = sM(:) == seeds(iroi); 
end

mn_roi.tc = zeros(length(seeds), size(Y, 2));
for iroi = 1:length(seeds)
    indices = find(seedMat(:, iroi));
    mn_roi.tc(iroi,:)= mean(Y(indices, :));
    waitbar(iroi/length(seeds),hw);
end;

disp('Saving Timecourses');
%% Save the timecourses
save([ffiles(1:end-4), '_', mn_roi.name,'_mn_tc.mat'],'mn_roi');
csvwrite([ffiles(1:end-4), '_', mn_roi.name,'_mn_tc.csv'], mn_roi.tc);
waitbar(0,hw,'Saving mean timecourses');
close(hw);