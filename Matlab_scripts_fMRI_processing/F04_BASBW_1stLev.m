


%% required:
% RobustWLS toolbox: http://www.diedrichsenlab.org/imaging/robustWLS.html
% for installation (basically just put the rwls folder into the toolbox
% folder in your SPM folder.

cd('/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts/');
clear all;

addSPM; % function that adds SPM. edit Path when you need to.
rmpath('/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts/spm8_deOrthog/'); % why is this removed?

overwritePrevSession = 0  % CHECK THIS EACH TIME!

runBayesian = 0; %DEFAULT!!


concatAnalysis = 1;         % Turning off for localizer

curFileStub = 'BASB'; % changed
numBlocks = 2;

excludeMotionControl = 0;   % Forcing inclusion for localizer
rmMotionOutliers = 1;

marsExtractPrep = 0;

usingRWLS = 1;

if usingRWLS
    origTemplatePath = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts/BASBGenericLev1full_rwls.mat'; % !!!! make new template (e.g. run1stlevel) and edit; check
    excludeMotionControl = 0;   % Forcing inclusion for localizer
    rmMotionOutliers = 0;
else
    origTemplatePath = '~/matlab/UKBlev1.mat';
end

excludeSubRuns(1,:) = [9999,6];  % PLACEHOLDER

% disp('DEORTHOG TURNED OFF!!') % "Should be deleted" AS via slack.
deorthogAnalysis = 2;  % for 2, only using different spm script! % unclear what this means


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
 useSTCsl1 = 0; % That was something about slice timing correction that we don't have for now

    usingEv = 1;    % 0 = full dur 
%   usingPrePrompt = 0    % 0 = CHOICE, -1 = RESPONSE % using response
    usingRespOnset = 1; % set to 1 to use response onset
%   Onset 0 or 1; would determine which GLM is being run and also
%   goes into file name

    allCondTypes = {'BAS_AllTrials'}; % adjust SOT names? Looks like it's already somewhat overlapping
    
    % 'BAS_pChGoalBin','BAS_pChGoalBin_pAvBid_pChGlAvBd','BAS_cBW', 'BAS_cBW_pAvBid', 'BAS_cBW_pVD', 'BAS_cBW_pLike', , 'BAS_pRewvGoalRT'
% only run 
 %allCondTypes = {'BAS_AllTrials'};
 %allCondTypes = {'BAS_pRewvGoalRT'};
    % subjects we want in no analysis
    allexclude =[7005]; % rm 7026 once we have the data...
    condexclude = [7004];
    excludeCond = 'BAS_cBW_pLike';
%%
            
%   allCondTypes = {'BAS_AllTrials'}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if excludeMotionControl && rmMotionOutliers %% will most likely not do this here... let's see
    multRegType = 'MCpars_outOnly'
else
    multRegType = 'MCparsBasic'
end
% relevant to motion correction (we're not using cause we to rwls)
MexclSDthresh = 4; % 4 stds
MexclMMthreshConj = 1.00; % Increased from DTS %  in mm (in conjunction w std)
MexclMMthreshHard = 2.00; % Increased from DTS % in mm irrespective of std threshold

basepath = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/508/';
basepathBehav = '~/Dropbox (Brown)/ShenhavLab/experiments/bas/Results/';

chdir(basepath);

startTR = 5;  % KEEP IN MIND THAT the 1st 4 currently still exist in each file!

TRlength = 2.5; % secs


allSubs = dir([basepathBehav,'BASB7*']); 

% we might already have renamed the 7005 file so it doesn't show up in this
% list anyway

% To exclude bad subjects or to run subsets of subjects (e.g. for parallel
% processing)
subRange = [1 nan]; % check later whether this is useful here
if isnan(subRange(2))
    subEnd = length(allSubs);
else
    subEnd = subRange(2);
end
allSubs = {allSubs(subRange(1):subEnd).name};

if runBayesian
    disp('RUNNING AS BAYESIAN ANALYSIS!!!!');
end



%%
for cti = 1:length(allCondTypes)
    condType = allCondTypes{cti}
    
    excludeSub = allexclude;
    if strcmp(condType, excludeCond)
       excludeSub= [excludeSub , condexclude];
    end
    % all of this just just creates names based on whatever the current
    % input is and what kind of analysis we want to run.
    if usingEv
        condType = [condType,'_Ev'];
    end
     % commented out as replaced by next if...
%     if usingPrePrompt==1 %change to usingRespOnset (or something)
%         condType = [condType,'_Pre'];
%     elseif usingPrePrompt==-1
%         condType = [condType,'_Resp'];
%     end
    % % % % added 05/14 bc... cb...
    if usingRespOnset==1 %change to usingRespOnset (or something)
        condType = [condType,'_Resp'];
    end
    
    condType_noDeO = condType
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
    if runBayesian
        condType = [condType,'_BMS'];
    end
    
    % now "Budder bei die Fische" we actually do something for each
    % subject.
    for subi = 1:length(allSubs)
        tic
        clear SPM; % delete the SPM file that's potentially currently in wksp
        allOkayBlocks = [];
        
        allMparDiffs = [];
        allMexclRegs = [];
        allMexclRegs_byBlock = [];
        try
            subi
            curSub = allSubs{subi}
            curSubNum = str2num(curSub(5:8)); % adjusted for BASB
            
            if ismember(curSubNum, excludeSub)
                %exist(fullfile(basepath,curSub,'EXCLUDE.mat'),'file')
                display(['EXCLUDING ',curSub]);
            else
                load(origTemplatePath);
                
                % Forcing TR!
                % why is this specified for RT? Did they actually mean TR?
                if ~usingRWLS
                    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TRlength;
                else
                    matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.timing.RT = TRlength; % ok, that's already there
                end
                curSubBehavMat = dir(fullfile(basepathBehav,[curSub(1:8),'*mat'])); % check whether that's actually working
                curSubBehavMat = curSubBehavMat(1).name;
                load(fullfile(basepathBehav,curSubBehavMat));
                

            % get the array for all runs and all RTs from logfiles to get
            % the relevant timing
            tmpAllRunNums = results.runNum;
            tmpAllRT_eval = results.metaResults1.CHOICE.RT;
            % then go through runs
            for run =1:2

                tmpCurRunTrialStart = results.timing.trialStartSecsRelative(tmpAllRunNums==run);
                tmpCurRunPostChoiceITI = results.timing.postChoiceITIs(tmpAllRunNums==run);
                tmpCurRunRTs = tmpAllRT_eval(tmpAllRunNums==run);

                tmpRunLength = tmpCurRunTrialStart(end)+...
                expParams.CHOICE.runEndISI+...
                tmpCurRunPostChoiceITI(end)+...
                tmpCurRunRTs(end);
            
                tmpBlockEnd(run) =  results.metaResults1.timing.runStartSecs(run)+ tmpRunLength;
                runEndISI_TRadj_add = max(0,expParams.CHOICE.TRlength-mod(tmpBlockEnd(run)+expParams.CHOICE.runEndISI-results.CHOICE.timing.scanBlockStartTTL(run),expParams.CHOICE.TRlength));
                tmpBlockEnd(run) = tmpBlockEnd(run)+runEndISI_TRadj_add;

                tmpRunLength = TRlength*ceil(tmpRunLength/TRlength);
                numActualTRsPerBlock(run)= round(( tmpBlockEnd(run)-results.metaResults1.timing.runStartSecs(run))/TRlength);
            end


                % Excluding dummy scans:
% %                 numActualTRsPerBlock = round((tmpBehav.results.metaResults1.CHOICE.timing.scanBlockEnd-tmpBehav.results.metaResults1.timing.runStartSecs)/tmpBehav.expParams.CHOICE.TRlength);
                % Including dummy scans
                TRsPerBlock = numActualTRsPerBlock+(startTR-1);
                allEndTRs = TRsPerBlock;
                
                %%
                % ASSUMING 2 BLOCKS!
                startTRperBlock = [1 numActualTRsPerBlock(1)+1];
                finalTRperBlock = [numActualTRsPerBlock(1) sum(numActualTRsPerBlock)];
                % ok, this needs to be changed such that the files are in
                % different run folders.
                
                %chdir(fullfile(basepath,curSub,'NII'));%% adjust, won't work for current folder structure
                % Hence, we go into the directory that has the different
                % epi directories
                chdir([basepath,num2str(curSubNum), '/Preprocessing']);
                EPIDIRS=dir(['EP2*']);
                clear allRuns;
                % we loop through those to find the relevant files and
                % append those
                for ii = 1:length(EPIDIRS)
                    
                if ~useSTCsl1 % this is where the slice timing flag becomes active
                     allRuns(ii) = dir([EPIDIRS(ii).name,'/swr*',curFileStub,'*nii']);
                    %allRuns=dir(['swr*',curFileStub,'*nii']);
                else
                    allRuns(ii) = dir([EPIDIRS(ii).name,'swarx*',curFileStub,'*nii']);
                    %allRuns=dir(['swarx*',curFileStub,'*nii']);
                end
                allMCparsTxt=dir([EPIDIRS(ii).name,'/r*',curFileStub,'*txt']);
                end
               
                % because they are in different folders, we need to
                % actually add the folder to be able to load them later. 
                % This is what's happening here.
                allRunsMain = strcat({allRuns.folder},'/',{allRuns.name});
                allRuns = allRunsMain;
                
               % allMCparsTxt=dir(['r*',curFileStub,'*txt']); % This must
               % happen by folder as well! moved up into the loop above
                %allRunsMain = {allMCparsTxt(:).name} % same reasoning as
                %above, we also provide the location.
                allRunsMain=strcat({allMCparsTxt.folder},'/',{allMCparsTxt.name});
                allMCparsTxt = allRunsMain;
                
                curSubBlockExclusions = excludeSubRuns(find(excludeSubRuns(:,1)==curSubNum),2);
                
                block = 0;
                
                allExclTPbyBlock = {};
                
                allRunsNiiFilesCat = {};
                for runi = 1:length(allRuns)
                    if isempty(find(curSubBlockExclusions==runi))
                        block = block+1;
                        blockOkay = 1;
                    else
                        blockOkay = 0;
                        % EXCLUSION CURRENTLY ONLY WORKS IF THERE ARE
                        % NII FILES FOR ALL BLOCKS (INCL. EXCLUDED ONES)
                    end
                    
                    if blockOkay
                        allOkayBlocks = [allOkayBlocks runi];
                        curRunNiiFiles = {};
                        %curRunPath = fullfile(basepath,curSub,'NII');%
                        % maybe I don't actually need this bc filename
                        %includes folder...
                        curRunNii = allRuns{runi};
                        curNiiHdr = spm_vol_nifti(curRunNii);
                        curNumTRs = curNiiHdr.private.dat.dim(4);
                        trInd = 1;
                        for tri = startTR:allEndTRs(runi)%curNumTRs  % CHANGED!
                            %curRunNiiFiles{trInd} = [fullfile(curRunPath,curRunNii),',',num2str(tri)];
                            curRunNiiFiles{trInd} = [fullfile(curRunNii),',',num2str(tri)];
                            trInd = trInd + 1;
                        end
                    end
                    allRunsNiiFilesCat = [allRunsNiiFilesCat,curRunNiiFiles];
                    
                    if concatAnalysis
                        %curSubCondMat =  fullfile(basepath,curSub,'SPM','SOTS',[condType_noDeO,'_sots_allTcat.mat']);
                        curSubCondMat =  fullfile(basepath,num2str(curSubNum),'SPM','SOTS',[condType_noDeO,'_sots_allTcat.mat']);

                    else
                        %curSubCondMat =  fullfile(basepath,curSub,'SPM','SOTS',[condType_noDeO,'_sots_run',num2str(block),'.mat']);
                        curSubCondMat =  fullfile(basepath,num2str(curSubNum),'SPM','SOTS',[condType_noDeO,'_sots_run',num2str(block),'.mat']);

                    end
                    
                    if concatAnalysis && length(allRuns)>1 %%&& lev1Seg<10
                        %curSubRegMat =  fullfile(basepath,curSub,'SPM',[multRegType,'_AllRunsConcat.mat']);
                        curSubRegMat =  fullfile(basepath,num2str(curSubNum),'SPM',[multRegType,'_AllRunsConcat.mat']);

                    else
                        %curSubRegMat =  fullfile(basepath,curSub,'SPM',[multRegType,'_run',num2str(block),'.mat']);
                        curSubRegMat =  fullfile(basepath,num2str(curSubNum),'SPM',[multRegType,'_run',num2str(block),'.mat']);

                    end
                    if blockOkay
                        if ~concatAnalysis
                             if ~exist(curSubRegMat,'file') && rmMotionOutliers
                                %%%% SAVING EACH TIME JUST IN CASE
                                R=load(allMCparsTxt{runi});
                                % CURRENTLY adjusted under assumption that initial TRs WERE excluded from MC
                                R = R(1:(allEndTRs(runi)-(startTR-1)),:);
                                save(curSubRegMat,'R');
                             end
                             
                             
                             if ~usingRWLS
                                 % Copying 1st block config
                                 matlabbatch{1}.spm.stats.fmri_spec.sess(block) = matlabbatch{1}.spm.stats.fmri_spec.sess(1);
                                 
                                 matlabbatch{1}.spm.stats.fmri_spec.sess(block).multi = {curSubCondMat};
                                 if excludeMotionControl && ~rmMotionOutliers
                                     matlabbatch{1}.spm.stats.fmri_spec.sess(block).multi_reg = {''};
                                 else
                                     matlabbatch{1}.spm.stats.fmri_spec.sess(block).multi_reg = {curSubRegMat};
                                 end
                                 matlabbatch{1}.spm.stats.fmri_spec.sess(block).scans = curRunNiiFiles';
                                 matlabbatch{1}.spm.stats.fmri_spec.sess(block).hpf = 128;
                             else
                                 matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(block) = matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(1);
                                                                  
                                 matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(block).multi = {curSubCondMat};
                                 
                                 if excludeMotionControl && ~rmMotionOutliers
                                     matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(block).multi_reg = {''};
                                 else
                                     matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(block).multi_reg = {curSubRegMat};
                                 end                                 
                                 
                                 matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(block).scans = curRunNiiFiles';
                                 matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(block).hpf = 128;
                           end

                        end
                        if rmMotionOutliers 
                            if runi==1
                                for xruni = 1:length(allRuns)
                                    tmpMpars = load(allMCparsTxt{xruni});
                                    % x/y/z TRANSLATIONAL movement
                                    tmpMdiffs = max(abs([zeros(1,3);diff(tmpMpars(:,1:3))]),[],2)';
                                    allMparDiffs = [allMparDiffs,tmpMdiffs];
                                end
                                for xxruni = 1:length(allRuns)
                                    tmpBlockMdiffs = allMparDiffs(startTRperBlock(xxruni):finalTRperBlock(xxruni));
                                    allExclTPbyBlock{xxruni} = find((tmpBlockMdiffs>MexclMMthreshConj & tmpBlockMdiffs>=(std(allMparDiffs)*MexclSDthresh)) | tmpBlockMdiffs>MexclMMthreshHard);
                                end
                                if concatAnalysis
                                    tmpExclTPs_concat = find((allMparDiffs>MexclMMthreshConj & allMparDiffs>=(std(allMparDiffs)*MexclSDthresh)) | allMparDiffs>MexclMMthreshHard);
                                end
                            end
                            
                            if ~concatAnalysis
                                curExcludeTPs = allExclTPbyBlock{block};
                                if isempty(curExcludeTPs)
                                    if ~usingRWLS
                                        matlabbatch{1}.spm.stats.fmri_spec.sess(block).multi_reg = {''};
                                    else
                                        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(block).multi_reg = {''};
                                    end
                                else
                                    tmpOutRegs = zeros(numActualTRsPerBlock(block),length(curExcludeTPs));
                                    for cici=1:length(curExcludeTPs)
                                        tmpOutRegs(curExcludeTPs(cici),cici) = 1;
                                    end
                                    R = tmpOutRegs;
                                    save(curSubRegMat,'R');
                                    if ~usingRWLS
                                        matlabbatch{1}.spm.stats.fmri_spec.sess(block).multi_reg = {curSubRegMat};
                                    else
                                        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(block).multi_reg = {curSubRegMat};
                                    end
                                end
                            end
                        end
                    end
                end
                
                if ~usingRWLS
                    tmpSessLen = length(matlabbatch{1}.spm.stats.fmri_spec.sess);
                else
                    tmpSessLen = length(matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess);
                end

                if tmpSessLen>block && ~concatAnalysis
                    if ~usingRWLS
                        matlabbatch{1}.spm.stats.fmri_spec.sess = matlabbatch{1}.spm.stats.fmri_spec.sess(1:block);
                    else
                        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess = matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(1:block);
                    end
                end
                
                if concatAnalysis
                    if rmMotionOutliers
                        curExcludeTPs = tmpExclTPs_concat;
                        if isempty(curExcludeTPs)
                            tmpOutRegs = [];
                        else
                            tmpOutRegs = zeros(sum(numActualTRsPerBlock),length(curExcludeTPs));
                            for cici=1:length(curExcludeTPs)
                                tmpOutRegs(curExcludeTPs(cici),cici) = 1;
                            end
                        end
                        disp(['EXCLUDING ',num2str(length(curExcludeTPs)),'TRs!!!'])
                        
                        %curSubRegMat =  fullfile(basepath,curSub,'SPM',['CatCondMeansTrends_noMout.mat']);
                        curSubRegMat =  fullfile(basepath,num2str(curSubNum),'SPM',['CatCondMeansTrends_noMout.mat']);
                    else
                        tmpOutRegs = [];
                        %curSubRegMat =  fullfile(basepath,curSub,'SPM',['CatCondMeansTrends.mat']);
                        curSubRegMat =  fullfile(basepath,num2str(curSubNum),'SPM',['CatCondMeansTrends.mat']);
                    end
                    
                    %curSubCondMat =  fullfile(basepath,curSub,'SPM','SOTS',[condType_noDeO,'_sots_allTcat.mat']);
                    curSubCondMat =  fullfile(basepath,num2str(curSubNum),'SPM','SOTS',[condType_noDeO,'_sots_allTcat.mat']);

                    if ~usingRWLS
                        matlabbatch{1}.spm.stats.fmri_spec.sess = matlabbatch{1}.spm.stats.fmri_spec.sess(1);
                        matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {curSubCondMat};
                        matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = allRunsNiiFilesCat';
                        % Turning off HPF:
                        matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = inf;
                    else
                        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess = matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(1);
                        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(1).multi = {curSubCondMat};
                        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(1).scans = allRunsNiiFilesCat';
                        % Turning off HPF:
                        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(1).hpf = inf;
                    end
                    if ~exist(curSubRegMat,'file')
                        tmpAllMeans = zeros(sum(numActualTRsPerBlock),numBlocks-1);
                        tmpAllTrends = zeros(sum(numActualTRsPerBlock),numBlocks-1);
                        for runi = 1:length(allRuns)
                            if runi==1
                                curStartTR =1;
                                curEndTR = numActualTRsPerBlock(runi);
                            elseif runi==2
                                curStartTR =numActualTRsPerBlock(runi-1)+1;
                                curEndTR = sum(numActualTRsPerBlock);                                
                            else
                               error(); 
                            end
                            if runi<length(allRuns)
                                tmpAllMeans(curStartTR:curEndTR,runi) = 1;
                            end
                            tmpAllTrends(curStartTR:curEndTR,runi) = linspace(-1,1,numActualTRsPerBlock(runi));
                        end
                        R=[tmpAllMeans,tmpAllTrends,tmpOutRegs];
                        save(curSubRegMat,'R');
                    end
                    if ~usingRWLS
                        matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {curSubRegMat};
                    else
                        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(1).multi_reg = {curSubRegMat};
                    end
                end
                
                %chdir(fullfile(basepath,curSub,'SPM'));
                chdir(fullfile(basepath,num2str(curSubNum),'SPM'));
                if exist(condType,'dir') && overwritePrevSession
                    display(['DELETING PREVIOUS LEVEL 1 FILES FOR ',curSub]);
                    unix(['rm -Rf ',condType]);
                end
                if ~(exist(condType,'dir') && ~overwritePrevSession)  % Not skipping previous sessions
                    mkdir(condType); chdir(condType);
                    curCondDir = pwd;
                    
                    if ~usingRWLS
                        matlabbatch{1}.spm.stats.fmri_spec.dir = {curCondDir};
                    else
                        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.dir = {curCondDir};
                    end
                    
                    chdir(fullfile(basepath,num2str(curSubNum)));
                    
                    if deorthogAnalysis~=1
                        if runBayesian % CURRENTLY ONLY SET UP FOR DEO == 0 or 2
                            matlabbatch(2) = [];
                            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
                            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
                            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
                            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
                            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
                            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
                            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'fMRI model specification: SPM.mat File';
                            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
                            matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
                            matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.space.volume.block_type = 'Slices';
                            matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.signal = 'UGL';
                            matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.ARP = 3;
                            matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.noise.UGL = 1;
                            matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.LogEv = 'Yes';
                            matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.anova.first = 'No';
                            matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.anova.second = 'Yes';
                            matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.gcon = struct('name', {}, 'convec', {});
                        end
                        %num2str(curSubNum)
                        curBatchFile = [curSub(1:8),'_spm8_1stLev',condType,'.mat'];
                        save(curBatchFile,'matlabbatch');
                        % Specifying AND estimating SPM:
                        
                        if deorthogAnalysis==2
                            % deorth version of spm design:
                            % ok, that's a little ugly, but.... whatever...
                            addpath('/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts/spm8_deOrthog/');
                        else
                            rmpath('/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts/spm8_deOrthog/');
                        end
                        %spm('defaults', 'FMRI');
                        %spm_jobman('initcfg');
                        spm_jobman('run',  curBatchFile)
                        %spm_jobman('run',  matlabbatch)
                    else
                        allSPMs = [];
                        orig_matlabbatch = matlabbatch;
                        matlabbatch(2) = [];
                        
                        curBatchFile = [curSub(1:8),'_spm8_1stLev',condType,'.mat'];
                        save(curBatchFile,'matlabbatch');
                        % Specifying but not estimating SPM:
                        spm_jobman('run',  curBatchFile)
                        
                        tmpSPM = load(fullfile(curCondDir,'SPM.mat'));
                        allSPMs{1} = tmpSPM.SPM;
                        delete(fullfile(curCondDir,'SPM.mat'));
                                                
                        if deorthogAnalysis==1
                            lengthRegs = nan; % how many regs per block
                            for namei = 1:length(tmpSPM.SPM.xX.name)
                                if isnan(lengthRegs) && strcmp([tmpSPM.SPM.xX.name{namei}(1:5)],'Sn(2)')
                                    lengthRegs = namei-1;
                                end
                            end
                            
                            for deorthogIter = 1:(numDeorthogParams-1)
                                num2str(curSubNum)
                                for runi = 1:length(allOkayBlocks)  % RUNI SHOULD MAP ONTO BLOCKS - need to double-check!
%                                     curSubCondMatOrig =  fullfile(basepath,curSub,'SPM',[condType_noDeO,'_sots_run',num2str(runi),'.mat']);
%                                     curSubCondMatNext =  fullfile(basepath,curSub,'SPM',[condType_noDeO,'_nextDeO_sots_run',num2str(runi),'.mat']);
                                    curSubCondMatOrig =  fullfile(basepath,num2str(curSubNum),'SPM',[condType_noDeO,'_sots_run',num2str(runi),'.mat']);
                                    curSubCondMatNext =  fullfile(basepath,num2str(curSubNum),'SPM',[condType_noDeO,'_nextDeO_sots_run',num2str(runi),'.mat']);
                                    if deorthogIter==1
                                        tmpPrevCondMat = load(curSubCondMatOrig);
                                    else
                                        tmpPrevCondMat = load(curSubCondMatNext);
                                    end
                                    % onsets should remain the same:
                                    tmpPrevCondMat.pmod(1).name = circshift(tmpPrevCondMat.pmod(1).name,[0 -deorthogIter]);
                                    tmpPrevCondMat.pmod(1).poly = circshift(tmpPrevCondMat.pmod(1).poly,[0 -deorthogIter]);
                                    tmpPrevCondMat.pmod(1).param = circshift(tmpPrevCondMat.pmod(1).param,[0 -deorthogIter]);
                                    tmpPrevCondMat.pmod(2).name = circshift(tmpPrevCondMat.pmod(2).name,[0 -deorthogIter]);
                                    tmpPrevCondMat.pmod(2).poly = circshift(tmpPrevCondMat.pmod(2).poly,[0 -deorthogIter]);
                                    tmpPrevCondMat.pmod(2).param = circshift(tmpPrevCondMat.pmod(2).param,[0 -deorthogIter]);
                                    names = tmpPrevCondMat.names; onsets = tmpPrevCondMat.onsets; durations = tmpPrevCondMat.durations; pmod = tmpPrevCondMat.pmod;
                                    save(curSubCondMatNext,'onsets','pmod','durations','names');
                                    
                                    if ~usingRWLS
                                        matlabbatch{1}.spm.stats.fmri_spec.sess(runi).multi = {curSubCondMatNext};
                                    else
                                        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.sess(runi).multi = {curSubCondMatNext};
                                    end
                                end
                                
                                curBatchFile = [curSub(1:8),'_spm8_1stLev',condType,'.mat'];
                                save(curBatchFile,'matlabbatch');
                                % Specifying but not estimating SPM:
                                spm_jobman('run',  curBatchFile)
                                
                                tmpSPM = load(fullfile(curCondDir,'SPM.mat'));
                                allSPMs{deorthogIter+1} = tmpSPM.SPM;
                                delete(fullfile(curCondDir,'SPM.mat'));
                            end
                            
                            SPM = allSPMs{1}; % starting w/ original, then modifying
                            
                            %                     keyboard
                            for dori = 1:length(deorthogRegs)
                                curDOR = deorthogRegs(dori); %CPI
                                if length(deorthogRegs)==4 && dori<=2
                                    useDOR =  deorthogRegs(1)-1;
                                elseif length(deorthogRegs)==4 && dori>=3
                                    useDOR =  deorthogRegs(3)-1;
                                elseif length(deorthogRegs)==2 && dori==1
                                    useDOR = deorthogRegs(1)-1;
                                elseif length(deorthogRegs)==2 && dori==2
                                    useDOR = deorthogRegs(2)-1;  % RT
                                end
                                if length(deorthogRegs)==4 && (dori==2 || dori==4)
                                    useSPMi = 3;
                                else
                                    useSPMi = 2;
                                end
                                SPM.xX.X(:,curDOR:lengthRegs:(lengthRegs*length(allOkayBlocks))) = allSPMs{useSPMi}.xX.X(:,(useDOR):lengthRegs:(lengthRegs*length(allOkayBlocks)));
                            end
                        end
                        save(fullfile(curCondDir,'SPM.mat'),'SPM');
                        clear matlabbatch;
                        matlabbatch = {orig_matlabbatch{2}};
                        matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(curCondDir,'SPM.mat')};
                        curBatchFile = [curSub,'_spm8_1stLevEstOnly',condType,'.mat'];
                        save(curBatchFile,'matlabbatch');
                        spm_jobman('run',  curBatchFile)
                    end
                else
                    display(['SKIPPING ',curSub]);
                end
                
                clear matlabbatch; clear results; clear p; clear expParams;%clear tmpBehav;
            end
        catch me
            keyboard;
        end
        
        rmpath('/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts/spm8_deOrthog/');
        toc
    end
    
end

rmpath('/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts/spm8_deOrthog/');



