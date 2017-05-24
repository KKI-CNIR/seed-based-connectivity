%procdir: directory with preprocessed data
%tooldir: directory where processing scripts are located
%outdir: directory to write clustering files

if(ismac)
    servdir = '/Volumes/fmri/RestingState';
    procdir = '/Users/nebel/Projects/preprocessing/processed_data';
    tooldir = '/Users/nebel/Projects/preprocessing/fmri_preproc_toolbox';
    expdir = '/Users/nebel/Projects/PCG';
    shrinkdir = fullfile(expdir, 'scripts', 'shrinkIt-master');
elseif(isunix)
    addpath('/matlab/spm12/');
    procdir = '/mnt/bisivity/fmri/RestingState/Kids/rsData';
    tooldir = '/mnt/bisivity/fmri/RestingState/Jobs/sj_scripts/fmri_preproc_toolbox';
else
    procdir = '\\kki-gspnas1\LNIR_DATA$\fmri\RestingState\Kids\rsData';
    tooldir = '\\kki-gspnas1\LNIR_DATA$\fmri\RestingState\Jobs\sj_scripts\fmri_preproc_toolbox\';
end

%slist: path/name.txt containing list of IDs to be processed
%istart: line of slist to start on
%task: what to look for in each subject's raw data folder
%rename_task: e.g., task-rest_bold or task-GNG_bold; processing script will copy
%       and rename func_files to procdir/ID using the following convention -
%       sub-ID/sess-DateOfSession/func/run-??/sub-ID_sess-DateOfSession_rename_task_run-??.nii,
%       where run number is determined by the order of func_files
%specs_file: path/file.m containing processing preferences

slist = fullfile(expdir, 'scripts', 'IDsforMotion.txt');
istart = 1;
% iend = 4;
task = 'task-rest_bold';
prefix = 'fnwc50fwepia';
roifile = fullfile(expdir, 'scripts', 'wRL_precentral_gyrus.nii');
brainmask_file = fullfile(expdir, 'scripts', 'brain-PCG.nii');
k = 5;

addpath(tooldir)
addpath(shrinkdir)
%open slist and get functional IDs (SM numbers)
subjects = readtable(slist, 'Delimiter', ' ');

SMs = char(subjects.Var2);
scanDates = num2str(subjects.Var3);
rns = char(subjects.Var4);
alength = size(SMs, 1);
keepIDs = ones(alength, 1);
[~, slistname, ~] = fileparts(slist);

if(~exist('iend', 'var'))
    iend = length(SMs);
end

cwd = pwd;

%keep track of IDs in list missing data
fid2 = fopen(fullfile(cwd, 'ID_DOS_missing_parrec.txt'), 'w');

%load roifile
R = spm_vol(roifile);
roi = spm_read_vols(R);
roi = roi(:);
roi_indx = find(roi);

%load brain-PCG
B = spm_vol(brainmask_file);
bm = spm_read_vols(spm_vol(B));
bm_indx = find(bm);

%loop over subjects in slist
for isub = istart:iend
    rawSM = SMs(isub, :);
    
    %add leading zeros to make all IDs the same length
    ID = strcat('sub-', num2str(str2num(rawSM), '%04d'));
    
    DOS = strcat('ses-', scanDates(isub, :));
    fprintf('Working on: %s/%s\n', ID, DOS)
    
    %4D preprocessed data for the current subject
    ffile = fullfile(procdir, ID, DOS, 'func', strcat('run-0', rns(isub)), ...
        strcat(prefix, ID, '_', DOS, '_', task, '_run-0', rns(isub),'.nii'));
    
    if(exist(ffile, 'file'))
        subdir = fullfile(expdir, ID, DOS, strcat('run-0', rns(isub)));
        
        %if subject outdir does not exist, create it
        if(~exist(subdir, 'dir'))
            [SUCCESS, MESSAGE, MESSAGEID] = mkdir(subdir);
        end
        
        %load 4D data
        P = spm_vol(fullfile(procdir, ID, DOS, 'func', strcat('run-0', rns(isub)), ...
            strcat(prefix, ID, '_', DOS, '_', task, '_run-0', rns(isub),'.nii')));
        nframes = length(P);
        fdata = spm_read_vols(P);
        
        %reshape fdata
        vByT = double(reshape(fdata, [numel(fdata(:, :, :, 1)), numel(fdata(1, 1, 1, :))]));
        clear fdata;
        
        %extract timeseries for every voxel in ROI
        [tcs, ~] = seed_tcs_passData(vByT, nframes, roi, subdir);
        torder = 1:nframes;
        
        %apply brain-PCG mask to 4D data
        vByT = vByT(bm_indx, :);
        
        %calculate whole brain correlation for each ROI voxel
        [conn_map_all, comptime] = wb_conn_resample(vByT', nframes, ...
            tcs(:, torder)');
        
        %fisher transform correlations
        conn_map_all = fisher_r2z(conn_map_all);
        conn_map_all(conn_map_all > 5) = 5;
        conn_map_all(conn_map_all < -5) = -5;
        
        %calculate similarity matrix for voxels in ROI
        etafile = fullfile(subdir, strcat('raw_eta_sq_', prefix, '.mat'));
        [eta_sq, elpsd_time] = eta_square_resample(conn_map_all', 1, etafile);
        
        %cluster
        
        cl = sc_nodisp(eta_sq, 0 , k);
        
        % generate adjacency matrix for cluster solution
        adjacency = adjacency_matrix_resample(cl, 1, fullfile(subdir, strcat('adjacency_k', num2str(k), '_raw_', prefix, '.mat')));
        
        %write out labels
        %copy header information
        img_ca = R;
        %change file name in header
        img_ca.fname = fullfile(subdir, strcat('ca_k', num2str(k), '_raw_', prefix, '.nii'));
        img_ca.private.dat.fname = img_ca.fname;
        img_ca.dt = [16 0];
        
        Z = zeros(prod(img_ca.dim), 1);
        Z(roi_indx) = cl;
        
        Z = reshape(Z, img_ca.dim);
        spm_write_vol(img_ca, Z);
                
    else %no 4D preprocessed file; needs to copy raw data from godzilla & process
        warning(strcat('Missing ', task,' .rec from ', DOS))
        fprintf(fid2, '%s %s\n', rawSM, DOS);
        continue
    end
    
end

fclose(fid2);
