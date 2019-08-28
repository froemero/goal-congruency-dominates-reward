clear all; close all;
addSPM;
addpath('/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts'); % to make sure we have access to other functions (e.g. AS_ExtractVals_auto)

% Change only if excluding some initial subjects:
firstSub = 1; 
subset = 0; % this is to get only some ROIs... adjust which in line 43
% Total trials in session (DOUBLE-CHECK):
numTrials = 144;  
isStimLocked =1;
% Spit out individual con values? (leave as 0)
printText = 1;

% Where subject folders are found:
basedir = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/508/';

%%% Directory with ROIs:
rDir = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/508/BASBW_group/ROIs/'; % in group folder

%%%% Flag to find all relevant ROI files:
rtype = 'Anat_ProbAtlas_BA*4_p50bin_1.nii'; %'BARTRA_FULL*nii';%'BAS_S1_S2*nii';'OVr_L_R.nii'
% rtype = 'vmPFC*v_*_*nii'

% Currently assuming that we'll save wherever this is run (change this):
savedir = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Export/'; 
savefile = sprintf('%s_ROIS.mat', rtype(1:11)); %%% CHANGE THIS DEPENDING ON MASK NAME


%% Compile subject and ROI info
outsubs = {'7005'};
% Compile subject IDs: 
% We're not getting them from results, but 
chdir(basedir);
allResultsB = dir(['7*']);  %%% CHANGE THIS
allResultsB = {allResultsB(:).name};

allResultsB = allResultsB(~ismember(allResultsB, outsubs));

% Compile ROIs (masks):
Rs = dir(fullfile(rDir,rtype));
Rs = {Rs(:).name}
if subset
Rs = Rs(6:7);
end
for ri = 1:length(Rs)
    allMasks{ri} = fullfile(rDir,Rs{ri});
Rs{ri}
end


%% Run extraction

% Single 3-D matrix that will contain means within each ROI for each trial within each subject:
allMs = [];  % Will need to transorm to 2-D later

% For each trial:
for ci = 1:numTrials
    allFs{1}=[];
    
    % Assigning beta file number for this trial
    tmpBetaVal = num2str(ci/10000,'%.4f');
    tmpBetaVal = tmpBetaVal(3:6);
    
    for fi = 1:length(allResultsB)
        %%%% Change these paths (and I think img should be nii for SPM12 beta files but double check:
        try 
            if isStimLocked
            allFs{1}=[allFs{1};cellstr([basedir ,allResultsB{fi},'/SPM/BAS_AllTrials_Ev_deO_rwls_allTcat/beta_',tmpBetaVal,'.nii',' '])];
            else
            allFs{1}=[allFs{1};cellstr([basedir ,allResultsB{fi},'/SPM/BAS_AllTrials_Ev_Resp_deO_rwls_allTcat/beta_',tmpBetaVal,'.nii',' '])];
            end

        catch
            if isStimLocked
            allFs{1}=[allFs{1};[basedir ,allResultsB{fi},'/SPM/BAS_AllTrials_Ev_Resp_deO_rwls_allTcat/beta_',tmpBetaVal,'.nii']];
            else
            allFs{1}=[allFs{1};[basedir ,allResultsB{fi},'/SPM/BAS_AllTrials_Ev_Resp_deO_rwls_allTcat/beta_',tmpBetaVal,'.nii']];
            end
        end
    end
    
    % Extract (m)ean, (v)ariance, (n)umber of voxels, and (r)oi names for each ROI and each subject on the current trial:
    [m,v,n,r] = AS_ExtractVals_auto(allMasks,allFs,printText);
    
    allMs(ci,:,:) = m'; % assign values for all subjects and ROIs on the current trial;
end

% ROI names:
Rs = Rs';
% This is the main output we care about (all the extracted values):
allMs=permute(allMs,[3 2 1]);

save(fullfile(savedir,savefile),'allMs','Rs','allResultsB');
