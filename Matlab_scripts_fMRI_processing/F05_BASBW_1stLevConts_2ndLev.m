function F05_BASBW_1stLevConts_2ndLev

cd('/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts/');

    
    addSPM;

basepath = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/508/'; % path to fMRI data
%basepathBehav = '~/Dropbox (Brown)/ShenhavLab/experiments/bas/Results/'; % path to behavioral data
templatePATH = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts/2ndLevelStuff';
folderEnd = [];   % Something to append to end of folder  
curVersion = 0;

%%% Main things to adjust depending on what you want to run: %%%%%%%%%%%
doConts = 1;     % run 1st level (subject) contrasts
do2ndLev = 1;    % run 2nd level (group) t-tests over contrasts
do2ndLevReg = 0;  % 2nd level regression (e.g., indiv diff covariates)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Only change if excluding some subjects:
subRange = [1 nan];  % Can adjust first and last sub (use nan for last to go through end)
excludeSubs = {};   % Specific subjects to exclude

% Other settings you'll likely leave as is:
usingRWLS = 1
deorthogAnalysis = 2;  % 2 if deorthogonalizing, 0 if using default serial orthog

% Ignore these for now (these are forced to 0 for RWLS analyses)
excludeMotionControl = 0 
rmMotionOutliers = 0
if usingRWLS
    % Forcing to 0:
    excludeMotionControl = 0   
    rmMotionOutliers = 0
end

%%  GLM settings

concatAnalysis = 1; % 1 if concatenating

usingEv = 1         % 0 = full dur for For/Eng
usingRespOnset = 0  % if setting events to start relative to response onset
useSTCsl1 = 0       % if using slice-timing correction

%% GLMS and their corresponding contrasts (some examples shown)

% Can specify multiple GLMs, each should correspond to a 1st level analysis that has been run:
%allCondTypes = {'BAS_pChGoalBin','BAS_pChGoalBin_pAvBid_pChGlAvBd','BAS_cBW', 'BAS_cBW_pAvBid', 'BAS_cBW_pVD', 'BAS_cBW_pLike'}
%allCondTypes = {'BAS_pRewvGoalRT'};
allCondTypes = {'BAS_pChGoalBin','BAS_pChGoalBin_pAvBid_pChGlAvBd','BAS_cBW', 'BAS_cBW_pAvBid', 'BAS_cBW_pVD', 'BAS_cBW_pLike', 'BAS_pRewvGoalRT'};
excludeCond = 'BAS_cBW_pLike'; % we have missing data for this condition is considered further down
% Each value here corresponds to a contrast file below (make sure contrasts have been set up appropriately for your session type)
%allNumPs = [1.1, 1.3, 2.0, 2.1, 2.1, 2.1];
%allNumPs = [1.6];
allNumPs = [1.1, 1.3, 2.0, 2.1, 2.1, 2.1, 1.6];
  



%% For running 2nd level covariates (ignore for now):

% Old example for 2 individual difference variables:
% mregStruct.cname{1} = 'zBASsum';
% mregStruct.c{1} = [-0.851032727	-0.373818114	-0.61242542	0.819218419	-0.135210807	-1.089640033	1.057825726	0.580611113	1.296433033	0.342003806	-0.61242542	0.580611113	2.489469566	-1.566854647	0.580611113	0.342003806	-0.135210807	0.819218419	-0.373818114	1.535040339	0.1033965	1.296433033	-0.373818114	-1.32824734	0.342003806	-1.566854647	-0.851032727	-1.566854647	-0.373818114	-0.373818114];
% mregStruct.cname{2} = 'mean_zBIS_zNEON1';
% mregStruct.c{2} = [0.04337761	1.176766429	-0.428635099	1.06831629	0.151827749	-0.585972668	0.428290735	0.102940318	-1.96379864	1.117203721	0.319840596	-0.52640996	-1.708686209	0.162503027	0.900303443	-0.005509821	-0.537085238	-0.802872946	0.368728026	0.428290735	1.236329137	-0.683747529	-0.201059543	-2.082924056	0.477178165	0.802528582	1.816791985	-1.009097946	0.417615457	-0.439310376];
% 
% for mci=1:length(mregStruct.c)
%     mregStruct.c{mci} = AS_nanzscore(mregStruct.c{mci});
% end
% mregStruct.folderName = mregStruct.cname{1};

%% Start running subject-level contrasts and (optionally) group analysis

overwritePrevConts = 1;  % This will/won't run contrasts if they have already been run
overwritePrevSession = 1; % This will delete the previous session and run it again --> don't overwrite what we've done before

for cti = 1:length(allCondTypes)
    try
        condType = allCondTypes{cti}
        
        curNumP = allNumPs(cti)
        
         excludeSubs = [{'7005'}];
        
        if strcmp(condType, excludeCond)
            excludeSubs = [{'7004', '7005'}];
        end
        
        if curNumP==1.1
            origTemplatePath = [templatePATH,'/BAS_Generic1Stim1P_conts.mat'];
        elseif curNumP==1.2
            origTemplatePath = [templatePATH,'/BAS_Generic1Stim2P_conts.mat'];
        elseif curNumP==1.3
            origTemplatePath = [templatePATH,'/BAS_Generic1Stim3P_conts.mat'];
        elseif curNumP==1.4
            origTemplatePath = [templatePATH,'/BAS_Generic1Stim4P_conts.mat'];
        elseif curNumP==1.5
            origTemplatePath = [templatePATH,'/BAS_Generic1Stim5P_conts.mat'];
        elseif curNumP==1.6
            origTemplatePath = [templatePATH,'/BAS_Generic1Stim6P_conts.mat'];
        elseif curNumP==2.0
            origTemplatePath = [templatePATH,'/BAS_Generic2Stim0P_conts.mat'];
        elseif curNumP==2.1
            origTemplatePath = [templatePATH,'/BAS_Generic2Stim1P_conts.mat'];
        elseif curNumP==2.2
            origTemplatePath = [templatePATH,'/BAS_Generic2Stim2P_conts.mat'];
        elseif curNumP==2.3
            origTemplatePath = [templatePATH,'/BAS_Generic2Stim3P_conts.mat'];
        elseif curNumP==2.4
            origTemplatePath = [templatePATH,'/BAS_Generic2Stim4P_conts.mat'];
        elseif curNumP==4
            origTemplatePath = [templatePATH,'/BAS_Generic4Stim0P_conts.mat'];
        elseif curNumP==8
            origTemplatePath = [templatePATH,'/BAS_Generic8Stim0P_conts.mat'];
        else
            error('Bad numP!');
        end
        
        origTemplatePath
        
        if usingEv
            condType = [condType,'_Ev'];
        end
        if usingRespOnset
            condType = [condType,'_Resp'];
        end
        if deorthogAnalysis
            condType = [condType,'_deO'];
        end
        if usingRWLS
            condType = [condType,'_rwls'];
        end
        if excludeMotionControl
            condType = [condType,'_noMreg'];
        end
        
        if rmMotionOutliers
            condType = [condType,'_rmMout'];
        end
        if useSTCsl1
            condType = [condType,'_STCsl1'];
        end
        
        if concatAnalysis
            condType = [condType,'_allTcat'];
        end
        
        outFolder = [condType,folderEnd];
        
        
        chdir(basepath); 
        
        allSubs = dir(['7*']);
        if isnan(subRange(2))
            subEnd = length(allSubs);
        else
            subEnd = subRange(2);
        end
        %             allSubs = {allSubs(:).name};
        allSubs = {allSubs(subRange(1):subEnd).name};
        
       
        exSubInds = [];
        for exi = 1:length(excludeSubs)
            exSubInds = [exSubInds find(strncmp(excludeSubs{exi},allSubs,7))];
        end
        
        exSubInds
        allSubs(exSubInds) = [];
        
        % now move to fMRI path to do the other stuff...
        chdir(basepath);
        if doConts
            for subi = 1:length(allSubs)
                try
                    load(origTemplatePath);
                    
                    curSub = allSubs{subi};
                    
                    
                    tmpTconts = dir(fullfile(basepath,curSub,'SPM',condType,'spmT*nii')); % edited to match different folder naming
                    
                    if isempty(tmpTconts) || overwritePrevConts
                        
                        disp(['Running contrasts on Subject ',curSub,' for analysis ',outFolder]); % same
                        
                        chdir(fullfile(basepath,curSub));
                        
                        matlabbatch{1}.spm.stats.con.spmmat = {fullfile(basepath,curSub,'SPM',condType,'SPM.mat')};
                        
                        if curNumP==0
                            matlabbatch{1}.spm.stats.con.consess = matlabbatch{1}.spm.stats.con.consess(1);
                        elseif curNumP==1
                            matlabbatch{1}.spm.stats.con.consess = matlabbatch{1}.spm.stats.con.consess(1:2);
                        else
                            % Not changing anything
                            matlabbatch{1}.spm.stats.con.consess = matlabbatch{1}.spm.stats.con.consess;
                        end
                        matlabbatch{1}.spm.stats.con.delete = 1;
                        
                        curBatchFile = [curSub,'_spm8_1stLevConts_',condType,'.mat'];
                        save(curBatchFile,'matlabbatch');
                        
                        spm_jobman('run',  curBatchFile)
                    else
                        disp(['SKIPPING contrasts for Subject ',curSub,' for analysis ',outFolder]);
                    end
                catch me
                    keyboard
                    disp(['ERROR on Subject ',curSub,' for analysis ',outFolder]);
                end
                clear matlabbatch;
                
            end
        end
        % Main 2nd level analysis (t-tests over contrasts):
        if do2ndLev
            disp(['Running all 2nd level analyses for ',condType]);
            AS_spm8_UKB_2ndLev(condType,origTemplatePath,excludeSubs,overwritePrevSession,subRange,curVersion,curNumP,outFolder)
        end
        % 2nd level regression (optional)
        if do2ndLevReg
            if length(mregStruct.cname)>1
                mregStruct.folderName = [mregStruct.folderName,'_',num2str(length(mregStruct.cname)),'Regs'];
            end
            disp(['Running all 2nd level Multiple Regression analyses for ',condType]);
            AS_spm8_UKB_2ndLevReg(condType,origTemplatePath,excludeSubs,overwritePrevSession,subRange,curVersion,curNumP,outFolder,mregStruct)
        end
    catch me
        keyboard
        disp(['ERROR on analysis ',condType]);
    end
end




%% 2nd level analysis
function AS_spm8_UKB_2ndLev(condType,contTemplatePath,excludeSubs,overwritePrevSession,subRange,curVersion,curNumP,outFolder)

if ~exist('overwritePrevSession','var')
    overwritePrevSession = 0;
end

load(contTemplatePath);
if curNumP==0
    matlabbatch{1}.spm.stats.con.consess = matlabbatch{1}.spm.stats.con.consess(1);
elseif curNumP==1
    matlabbatch{1}.spm.stats.con.consess = matlabbatch{1}.spm.stats.con.consess(1:2);
else
    % Not changing anything
    matlabbatch{1}.spm.stats.con.consess = matlabbatch{1}.spm.stats.con.consess;
end
for cii = 1:length(matlabbatch{1}.spm.stats.con.consess)
    contNames{cii} = matlabbatch{1}.spm.stats.con.consess{cii}.tcon.name;
end
clear matlabbatch

basepath = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/508/'; % path to fMRI data
%basepathBehav = '~/Dropbox (Brown)/ShenhavLab/experiments/bas/Results/'; % path to behavioral data
if curVersion==0
    groupPath = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/508/BASBW_group/';
end
templatePATH = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts/2ndLevelStuff';
origTemplatePath = [templatePATH,'/UKB_spm8_2ndLev_template.mat'];
load(origTemplatePath);

chdir(basepath);

allSubs = dir(['7*']); % we actually get the folder names directly, so no need to use a prefix other than first digit
if isnan(subRange(2))
    subEnd = length(allSubs);
else
    subEnd = subRange(2);
end
allSubs = {allSubs(subRange(1):subEnd).name};

exSubInds = [];
for exi = 1:length(excludeSubs)
    exSubInds = [exSubInds find(strncmp(excludeSubs{exi},allSubs,8))];
end
exSubInds
allSubs(exSubInds) = [];

chdir(groupPath);
if overwritePrevSession && exist(outFolder,'dir')
    display(['DELETING PREVIOUS 2nd LEVEL CONDITION ',outFolder,'!']);
    unix(['rm -Rf ',outFolder]);
end
mkdir(outFolder);
chdir(outFolder);

try
    for subi = 1:length(allSubs)
        curSub = allSubs{subi}
        
        curPwd = pwd;
        chdir(fullfile(basepath,curSub,'SPM',condType)); % for MR folder, only use number (5:8)
        confNames = dir('con_0*nii'); % attention, new format!!! (was .img before)
        confNames = {confNames(:).name};
        chdir(curPwd);
        
        for conti = 1:length(contNames)
            % Model building:
%             keyboard
            matlabbatch{conti}.spm.stats.factorial_design.dir = {fullfile(groupPath,outFolder,contNames{conti})};
            
            
            matlabbatch{conti}.spm.stats.factorial_design.des.t1.scans{subi,1} =...
                fullfile(basepath,curSub,'SPM',condType,[confNames{conti},',1']);

            matlabbatch{conti}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{conti}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{conti}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{conti}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{conti}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{conti}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{conti}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            
            if subi==1
                mkdir(contNames{conti});
                % Model estimation:
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).sname = 'Factorial design specification: SPM.mat File';
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{conti}, '.','val', '{}',{1}, '.','val', '{}',{1});
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.method.Classical = 1;
               
                
                
                % Contrast:
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1) = cfg_dep;
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).sname = 'Factorial design specification: SPM.mat File';
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{conti}, '.','val', '{}',{1}, '.','val', '{}',{1});
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.consess{1}.tcon.name =contNames{conti} ;
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.consess{1}.tcon.convec = 1;
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.delete = 1;
                

                
            end
            
        end
    end
catch me
    keyboard
end

% keyboard


curBatchFile = ['Group_spm8_2ndLev_',outFolder,'.mat'];
save(curBatchFile,'matlabbatch');

spm_jobman('run',  curBatchFile)

clear matlabbatch;

return
%%




%% 2nd level covariates (e.g., individual differences) 
function AS_spm8_UKB_2ndLevReg(condType,contTemplatePath,excludeSubs,overwritePrevSession,subRange,curVersion,curNumP,outFolder,mregStruct)

try
    load(contTemplatePath);
    if curNumP==0
        matlabbatch{1}.spm.stats.con.consess = matlabbatch{1}.spm.stats.con.consess(1);
    elseif curNumP==1
        matlabbatch{1}.spm.stats.con.consess = matlabbatch{1}.spm.stats.con.consess(1:2);
    else
        % Not changing anything
        matlabbatch{1}.spm.stats.con.consess = matlabbatch{1}.spm.stats.con.consess;
    end
    for cii = 1:length(matlabbatch{1}.spm.stats.con.consess)
        contNames{cii} = matlabbatch{1}.spm.stats.con.consess{cii}.tcon.name;
    end
    clear matlabbatch
    
    basepath = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/508/'; % path to fMRI data
    %basepathBehav = '~/Dropbox (Brown)/ShenhavLab/experiments/bas/Results/'; % path to behavioral data
    if curVersion==0
        groupPath = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/508/BASBW_group/';
    end
    templatePATH = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts/2ndLevelStuff';
    origTemplatePath = [templatePATH, '/UKB_spm8_2ndLevReg_template.mat'];
    load(origTemplatePath);
        
    chdir(basepath);
    
    subjectNameRoot = '7*';
    allSubs = dir(subjectNameRoot);
    allSubs = {allSubs(:).name};
    
    exSubInds = [];
    for exi = 1:length(excludeSubs)
        exSubInds = [exSubInds find(strcmp(excludeSubs{exi},allSubs))];
    end
    allSubs(exSubInds) = [];
    for ii=1:length(mregStruct.c)
        mregStruct.c{ii}(exSubInds) = [];
    end

    chdir(groupPath);
    
    
    if exist('mregStruct','var')
        for mrii=1:length(mregStruct.cname)
            matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(mrii).cname = mregStruct.cname{mrii};
            matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(mrii).c = mregStruct.c{mrii};
        end
        mregFolder = mregStruct.folderName;
    else
        mregFolder = 'MultReg';
    end
    
    mkdir(outFolder);
    chdir(outFolder);
    
    if exist(mregFolder,'dir') && overwritePrevSession
        disp(['DELETING previous 2nd Level reg analysis named ',mregFolder])
        unix(['rm -Rf ',mregFolder]);
    end
    mkdir(mregFolder);
    chdir(mregFolder);
    
    for subi = 1:length(allSubs)
        curSub = allSubs{subi};
        
        curPwd = pwd;
        chdir(fullfile(basepath,curSub,'SPM',condType));
        confNames = dir('con_0*nii');
        confNames = {confNames(:).name};
        chdir(curPwd);
        
        for conti = 1:length(contNames)
            % Model building:
            matlabbatch{conti}.spm.stats.factorial_design.dir = {fullfile(groupPath,outFolder,mregFolder,contNames{conti})};
            
                        matlabbatch{conti}.spm.stats.factorial_design.des.mreg.scans{subi,1} =...
                fullfile(basepath,curSub,'SPM',condType,[confNames{conti},',1']);
            
            matlabbatch{conti}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{conti}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{conti}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{conti}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{conti}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{conti}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{conti}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            if subi==1
                mkdir(contNames{conti});
                
                for mrii=1:length(mregStruct.cname)
                    matlabbatch{conti}.spm.stats.factorial_design.des.mreg.mcov(mrii).c = matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(mrii).c;
                    matlabbatch{conti}.spm.stats.factorial_design.des.mreg.mcov(mrii).cname = matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(mrii).cname;
                    matlabbatch{conti}.spm.stats.factorial_design.des.mreg.mcov(mrii).iCC = matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).iCC;
                end
                
                matlabbatch{conti}.spm.stats.factorial_design.des.mreg.incint = matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint;
                
                % Model estimation:
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).sname = 'Factorial design specification: SPM.mat File';
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{conti}, '.','val', '{}',{1}, '.','val', '{}',{1});
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
                matlabbatch{length(contNames)+conti}.spm.stats.fmri_est.method.Classical = 1;
                
                % Contrast:
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1) = cfg_dep;
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).sname = 'Factorial design specification: SPM.mat File';
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{conti}, '.','val', '{}',{1}, '.','val', '{}',{1});
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');
                
                for mrii=1:length(mregStruct.cname)
                    matlabbatch{2*length(contNames)+conti}.spm.stats.con.consess{mrii}.tcon.name =[contNames{conti},'_',num2str(mrii)] ;
                    matlabbatch{2*length(contNames)+conti}.spm.stats.con.consess{mrii}.tcon.convec = [zeros(1,mrii-1) 1];
                    matlabbatch{2*length(contNames)+conti}.spm.stats.con.consess{mrii}.tcon.sessrep = 'none';
                end
                matlabbatch{2*length(contNames)+conti}.spm.stats.con.delete = 1;
                
            end
            
        end
    end
    
 
    
    curBatchFile = ['Group_spm8_2ndLevReg_',condType,'_',mregFolder,'.mat'];
    save(curBatchFile,'matlabbatch');
    
    spm_jobman('run',  curBatchFile)
    
    clear matlabbatch;
    
catch me
    keyboard
end

return




