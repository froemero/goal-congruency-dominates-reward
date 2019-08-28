%% About:
% This script does the preprocessing of fMRI data.
% specifics for preprocessing can be found in preprocessing_BASBW_job
%% Before running this script, make sure that BASBWfMRIdata.mat is up to
% date (created first running analyzeBASB_508_copy.m till line 918 and then
% running MakeTable508.m)! Also, make sure that the fMRI data folders are
% renamed (1st level= SubID, recorded data folder within that one = Raw) 
% and that F01 has been run on all subjects. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
cd('/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/508'); % go where data is
addpath('/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts'); % include Path to scripts so we can call shit
basedir='/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/508'; % save path to data

% get subject info
load('~/Dropbox (Brown)/ShenhavLab/experiments/bas/Analysis/Compiled/BASBWfMRIdata.mat');
% get valid subject IDs
sIDs= unique(allSubData.subj_idx);
% get directories
fIDs=dir([basedir '/7*']);
%%
for nsubs = 1: length(fIDs)
    currsub= categorical(str2double(fIDs(nsubs).name));
    if ismember(currsub, sIDs)
        % go to preprocessing folder of current subject
        cd(strcat(fIDs(nsubs).name, '/Preprocessing')) % go into current subject folder
        out_dir= [basedir,'/', fIDs(nsubs).name, '/Clean/'];
        % add Clean data folder if it doesn't exist yet - though we don't
        % really use that now...
        if ~(exist( out_dir,'dir')==7)
            mkdir(out_dir);
        end
        
        EPIDIRS=dir(['EP2*']); % get the converted data per run from the EPI folders we created in the first script
        for ii = 1:length(EPIDIRS)
            
            
            % for each EPI folder, get filenames
            % save into variable to make input for batch
            P = dir([EPIDIRS(ii).name,'/BASB*.nii']); % get all the .nii files starting w BASB in the corresponding folder
            ims = spm_select('expand', fullfile(P(1).folder,P.name)); % we need to do this to get all slices in the 4D nifti
            %P2 = (strcat(P(1).folder,'/',{P.name}')); % for 3D it worked
            %like this
            P2 = cellstr(ims); % convert to cells as required format for batch (before it's char vector)
            exp=['run', num2str(ii),'data = P2;']; % save into variable "runxx" data
            eval(exp)
            
        end
        % combine data of different runs (will be 1 by 2 cell with 1 by
        % slices cells in each)
        rawdatafiles = {run1data, run2data};
        preprocessing_BASBW_job; % get the batch that does that for us.
        matlabbatch{1}.spm.spatial.realign.estwrite.data = rawdatafiles; % replace datafiles in batch
        batchname= [basedir,'/', fIDs(nsubs).name, '/Preprocessing/preprocessing_job'];
        save(batchname, 'matlabbatch'); % save the batch for later reference
        
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');

        % here the actual preprocessing happens
        spm_jobman('run', matlabbatch);
        
       cd(basedir) % go back to basedir for the next subject
          
        
    else
        
        fprintf('%s excluded from data analysis. Make sure this is on purpose.\n', fIDs(nsubs).name)
        
    end
end