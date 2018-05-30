clear all; close all;

dirPath = 'Q:\fmri\RestingState\Kids\Projects\mm_workingdirectory\allsubs\'; 
% 
roipath = 'Q:\fmri\RestingState\Kids\Projects\mm_workingdirectory\seeds\';

RSsubject_info = [];
subjects = {'ss4','ss5','ss7','ss9','SS12','SS13','SS14','SS16','SS17','SS18','SS19','SS20', 'SS21','SS22','SS23','SS25',...
      'SS26','ss27','ss28','ss30','ss31','ss32','ss33','ss34','ss35','ss36','ss37','SS38'};
outDir = [dirPath 'results'];

ROI = {'leftcaudate','rightcaudate','leftputamen','rightputamen','bilatcaudate','bilatputamen'}

for d = 1:length(subjects)
    for r = 1:length(ROI)
        seed = ROI{r}
        ID = subjects{d};
        disp(ID);
        regionTC = roi_tc_seeds(dirPath,roipath,['fsnwc50wepia',ID,'_ses_session1_RS_run-01'],seed);
        save(['Q:\fmri\RestingState\Kids\Projects\mm_workingdirectory\Timecourses\',seed,'_',ID,'.mat'],'regionTC');   
    end
end

