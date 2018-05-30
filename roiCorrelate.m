clear all; close all;

% subjects = {'sub-1036', 'sub-1129', 'sub-1167', 'sub-1171', 'sub-1229', 'sub-1244', 'sub-1306',...
%     'sub-1309', 'sub-1383', 'sub-1395', 'sub-1406', 'sub-1483','sub-1410', 'sub-1425', 'sub-1459',...
%     'sub-1602', 'sub-1618', 'sub-1660', 'sub-1677', ...
%     'sub-1472', 'sub-1341', 'sub-1059'};
subjects = {'ss4','ss5','ss7','ss9','SS12','SS14','SS16','SS17','SS18','SS19','SS20', 'SS21',...
      'SS26','ss30','ss32','ss33','ss34','ss35','ss36','ss37','SS38'};

% ROI = {'leftcaudate','rightcaudate','leftputamen','rightputamen','bilatcaudate','bilatputamen'};
% ROI = {'OxfordStriatal_executive', 'OxfordStriatal_sensorimotor'};
% ROI = {'OxfordStriatal_executive','left_mfg'};
ROI = {'OxfordStriatal_executive', 'OxfordStriatal_sensorimotor', 'JuelichPreMotor_L', 'JuelichPreMotor_R', 'SalletSMA', 'JuelichM1a_L', 'JuelichM1a_R'};
dirPath = '/Users/nebel/Projects/stereotypy/data/'; 
TCpath = '/Users/nebel/Projects/stereotypy/Timecourses/';

% %load brainmask
% brainmask_file = '/Users/nebel/Projects/preprocessing/fmri_preproc_toolbox/brain_mask.nii'
% B = spm_vol(brainmask_file);
% bm = spm_read_vols(spm_vol(B));
% bm_indx = find(bm);

for d = 1:length(subjects);
    ID = subjects{d};
    %     ifile = [dirPath,'fsnwc50fwepia',ID,'_ses-CH1_task-rest_bold_run-01.nii'];
%     ifile = [dirPath,'fsnwc50wepia',ID,'_ses_session1_RS_run-01.nii'];
%     disp(ifile);
    TR = 2.5;
%     disp(TR);
    for r = 1:length(ROI);
        seed = ROI{r};
        load([TCpath,'Timecourses', seed,'_',ID,'.mat']);
        %         tps = 1:d;
        %         disp(tps)
        
        tcs(:,r) = regionTC.tc;
    end %r
    
    
  temp = corrcoef(tcs);
  sscorrs(d) = temp(1, 2);
    
    
end
