%% About:
% this script converts image files into niftis and sets up the structure of
% the relevant folders for preprocessing
% initially generated 3D niftis are converted to 4D niftis for flexibility
%% Before running this script, make sure that BASBWfMRIdata.mat is up to
% date (created running analyzeBASB_508_copy.m which also runs MakeTable508.m)! 
% Also, make sure that the fMRI data folders are
% renamed (1st level= SubID, recorded data folder within that one = Raw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
% go to data folder
cd('/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/508');
basedir='/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/508';
addpath('/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts'); % required e.g. for convert3Dto4D_job
addSPM;
load('~/Dropbox (Brown)/ShenhavLab/experiments/bas/Analysis/Compiled/BASBWfMRIdata.mat');
sIDs= unique(allSubData.subj_idx); % valid subject IDs based on table

fIDs=dir([basedir '/7*']); % all fMRI data folders

for nsubs = 31%27: length(fIDs)
    
    currsub= categorical(str2double(fIDs(nsubs).name));
    if ismember(currsub, sIDs) % exclude subject if it's not in table
        cd(strcat(fIDs(nsubs).name, '/Raw')) % jump into raw folder
        % overall output directory
        out_dir= [basedir,'/', fIDs(nsubs).name, '/Preprocessing/'];
        % add Preprocessing folder if it doesn't exist yet
        if ~(exist( out_dir,'dir')==7)
                mkdir(out_dir);     % make preprocessing folder if it doesn't exist yet
        end
        EPIDIRS=dir(['EP2*']); % find subfolders with EPI in them 
        % for now converting EPI files only
        for ii = 1:length(EPIDIRS) % for the number of EPI folders in raw
            % create directory if there is none yet (we append run number
            % so we can index that accordingly later)
            if ~(exist( [out_dir, EPIDIRS(ii).name(1:end-4), 'run', num2str(ii)],'dir')==7)
                mkdir([out_dir, EPIDIRS(ii).name(1:end-4), 'run', num2str(ii)]);
            
            % set path for output (i.e. directory we just created) 
            currout_dir = [out_dir, EPIDIRS(ii).name(1:end-4), 'run', num2str(ii)];
            % for each EPI folder, get everything that's in it, convert it and save
            % it into the preprocessing folder with run1 and run2 in the name
            P = dir([EPIDIRS(ii).name,'/*.IMA']); % get all the IMA files in the corresponding folder
            P2 = string(strcat(P(1).folder,'/',{P.name}')); % prepare for making hdr, required for conversion
            % cd(EPIDIRS(ii).name)  % full path provided, so we don't need
            % that
            hdr = spm_dicom_headers(P2);
            % this creates the 3D niftis
            spm_dicom_convert(hdr,'all','flat','nii' , currout_dir,0);
            
            %% 3D to 4D conversion
             cd(currout_dir) % go into the folder you just saved stuff in
            % get the data into format that can be batch scripted
            P = dir([currout_dir,'/*.nii']); % get all the .nii files in the corresponding folder
            P2 = (strcat(P(1).folder,'/',{P.name}')); % this is just needed to make the file names for the conversion
            exp=['data = P2;'];
            eval(exp)
            convert3Dto4D_job; % we have a job that does that, so we load it
            % then we modify the relevant datafile names
            matlabbatch{1}.spm.util.cat.vols = data;
            % and the output file names we want for that subject and run
            matlabbatch{1}.spm.util.cat.name = ['BASB', fIDs(nsubs).name, '_EPI_run', num2str(ii)];
            % cd ..
            spm('defaults', 'FMRI');
            spm_jobman('initcfg');

            % This is where the actual conversion happens by running the batch    
            spm_jobman('run', matlabbatch);
            % We then need to go back to Raw to do the same thing for the
            % other run
            cd(strcat(basedir,'/',fIDs(nsubs).name, '/Raw'))
            end
            
        end
        
        
        
        
        
        
        
        
        
      cd(basedir) % back to level 1 for next subject
    else
        
        fprintf('%s excluded from data analysis. Make sure this is on purpose.\n', fIDs(nsubs).name)
        
    end
end