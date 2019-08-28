clear all;

basepath = '~/Dropbox (Brown)/ShenhavLab/experiments/bas/';
basepathMR = '/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/508/';
basepathBehav = '~/Dropbox (Brown)/ShenhavLab/experiments/bas/Results/';
addpath('/Volumes/Research/CLPS_Shenhav_Lab/ScanningData/BASB/Scripts'); % to access e.g. AS_nanzscore and later potentially other useful functions

numBlocks = 2;
TRlength = 2.5; 

%% PLACEHOLDERS:
exSubs = [7005]; % rm 7026 here once we have that data
excludeSubRuns(1,:) = [9999,1];
numActualBlockTRs = nan; % Excluding initial dummy TRs (DEAL WITH BELOW)
%% Start group loop

chdir(basepath);

allResultsB = dir(['Results/BASB7*mat']); 
allResultsB = {allResultsB(:).name};

for subi=1:length(allResultsB)
    curSub = allResultsB{subi}
    curSubNum = str2num(curSub(5:8)); 
    
    chdir(basepathMR);
    curSubMR = curSub(5:8); 

    load(fullfile(basepathBehav,allResultsB{subi}));
    
    tmpSubNum = str2num(p.subID(1:4)); 
    
    curFinalExclude = ~isempty(find(exSubs==tmpSubNum));
    
    if ~curFinalExclude
        % moved so it doesn't crash when excluding subjects
            chdir(fullfile(curSubMR));
            if ~exist('SPM','dir')
                mkdir('SPM');
            end
            chdir('SPM');
        %% Start subject processing
        tmpAllRunNums = results.runNum; 
        tmpAllOns_eval = results.metaResults1.timing.progStartSecsRelative; 
        tmpAllRT_eval = results.metaResults1.CHOICE.RT; 
        
        % 
        tmpChooseGoal= ones(1,length(tmpAllRT_eval));
        if expParams.isBestFirst
        tmpChooseGoal(tmpAllRunNums==2) = tmpChooseGoal(tmpAllRunNums==1)*-1; %%%%%% FILL IN +1 for CHOOSE-BEST and -1 for CHOOSE-WORST
        else
        tmpChooseGoal(tmpAllRunNums==1) = tmpChooseGoal(tmpAllRunNums==2)*-1;    
        end
        
        % we need catches here, because of saving issues for some
        % participants
        try
        tmpAllRate_anx = results.SUBEVAL(1).Rating; 
        catch
         tmpAllRate_anx = nan(1,length(tmpChooseGoal));   
        end
        try
        tmpAllRate_conf = results.SUBEVAL(2).Rating; 
        catch
        tmpAllRate_conf = nan(1,length(tmpChooseGoal));    
        end
        try
        tmpAllRate_like = results.SUBEVAL(3).Rating; 
        catch
        tmpAllRate_like = nan(1,length(tmpChooseGoal));     
        end
        
        allTrProdChosen = results.metaResults1.CHOICE.prodNumChosen_wKeys; 
        
        %% Product set
        for tri=1:length(expParams.CHOICE.product1Bids)
       
            tmpTrBids = [expParams.CHOICE.product1Bids{tri} expParams.CHOICE.product2Bids{tri}];
            tmpTrBidsSorted = sort(tmpTrBids);
            allTrAvBid(1,tri) = mean(tmpTrBids);
            allTrMinBid(1,tri) = min(tmpTrBids);
            allTrMaxBid(1,tri) = max(tmpTrBids);
            allTrMaxvMinBid(1,tri) = max(tmpTrBids)-min(tmpTrBids);
            allTrMaxvNextBid(1,tri) = tmpTrBidsSorted(4)-tmpTrBidsSorted(3);
            allTrMaxvAvRemBid(1,tri) = tmpTrBidsSorted(4)-mean(tmpTrBidsSorted(1:3));
            allTrMinvNextBid(1,tri) = tmpTrBidsSorted(1)-tmpTrBidsSorted(2);
            allTrMinvAvRemBid(1,tri) = tmpTrBidsSorted(1)-mean(tmpTrBidsSorted(2:4));
            
            tmpChosenProdNum = results.CHOICE.prodNumChosen_wKeys(tri);

            allTrChosBidFIXED(1,tri) = results.BID(1).finalBid(find(expParams.BID.productOrder==tmpChosenProdNum));
            allTrChosvMinBidFIXED(1,tri) = allTrChosBidFIXED(tri)-min(tmpTrBids);
            allTrSumUnchBidFIXED(1,tri) = sum(tmpTrBids)-allTrChosBidFIXED(tri);
            allTrChosvAvRemBidFIXED(1,tri) = allTrChosBidFIXED(tri)-(allTrSumUnchBidFIXED(tri)/(expParams.CHOICE.numProductsPerChoice(1)-1));
        end
        
        tmpAllRate = results.metaResults1.CHOICE.Rating;  % check for -1; There can't be (no DL)
        tmpMissedTrials = (tmpAllRate==-1); 
        
        
        %% Missed trial nan-ing:
        tmpAllRate_anx(tmpMissedTrials) = nan;
        tmpAllRate_conf(tmpMissedTrials) = nan;
        tmpAllRate_like(tmpMissedTrials) = nan;
        %%%%% etc., etc.
        
        
        %% make goal vs av Rem from max/min vs avRem and choose best
        
            allTrgoalvAvRemBid = zeros(1, length(allTrMaxvAvRemBid));
            allTrgoalvAvRemBid(tmpChooseGoal==1) = allTrMaxvAvRemBid(tmpChooseGoal==1);
            allTrgoalvAvRemBid(tmpChooseGoal==-1) = allTrMinvAvRemBid(tmpChooseGoal==-1);
            
            %% compute chosen vs unchosen in goal space + avgoal value
            
            allTrAvGoalV = allTrAvBid;
            allTrAvGoalV(tmpChooseGoal==-1) = 10-allTrAvBid(tmpChooseGoal==-1); % convert choose worst vals to goal space
            allTrGoalChosenvAvRem = allTrChosvAvRemBidFIXED;
            allTrGoalChosenvAvRem(tmpChooseGoal==-1)= allTrChosvAvRemBidFIXED(tmpChooseGoal==-1)*-1; % flip values for CW
            
            
        %% Demeaning/z-scoring:
            z_AvBid = AS_nanzscore(allTrAvBid)'; % OVr
            z_AvGoalV =  AS_nanzscore(allTrAvGoalV)';
            z_MinBid = AS_nanzscore(allTrMinBid)';
            z_MaxBid = AS_nanzscore(allTrMaxBid)';
            z_MaxvMinBid = AS_nanzscore(allTrMaxvMinBid)';
            z_MaxvNextBid = AS_nanzscore(allTrMaxvNextBid)';
            z_MaxvAvRemBid = AS_nanzscore(allTrMaxvAvRemBid)';
            z_GoalvAvRemBid = AS_nanzscore(allTrgoalvAvRemBid)'; 
            z_MinvNextBid = AS_nanzscore(allTrMinvNextBid)';
            z_MinvAvRemBid = AS_nanzscore(allTrMinvAvRemBid)';
            z_ChosBid = AS_nanzscore(allTrChosBidFIXED)';
            z_ChosvMinBid = AS_nanzscore(allTrChosvMinBidFIXED)';
            z_SumUnchBid = AS_nanzscore(allTrSumUnchBidFIXED)';
            z_ChosvAvRemBid= AS_nanzscore(allTrChosvAvRemBidFIXED)'; %VDr
            z_GoalChosenvAvRem= AS_nanzscore(allTrGoalChosenvAvRem)'; %VDg
            z_anx = AS_nanzscore(tmpAllRate_anx)';
            z_conf = AS_nanzscore(tmpAllRate_conf)';
            z_like = AS_nanzscore(tmpAllRate_like)';
            z_RT = AS_nanzscore(tmpAllRT_eval)';
        %% Timing:
        clear tmpBlockEnd;
        for run =1:2
            
            tmpCurRunTrialStart = results.timing.trialStartSecsRelative(tmpAllRunNums==run);
            tmpCurRunPostChoiceITI = results.timing.postChoiceITIs(tmpAllRunNums==run);
            tmpCurRunRTs = tmpAllRT_eval(tmpAllRunNums==run);
            
            
            % gets the last trial onset, adds the last RT and ITIs
            tmpRunLength = tmpCurRunTrialStart(end)+...
                expParams.CHOICE.runEndISI+...
                tmpCurRunPostChoiceITI(end)+...
                tmpCurRunRTs(end);
            
            
            tmpBlockEnd(run) =  results.metaResults1.timing.runStartSecs(run)+ tmpRunLength;
            runEndISI_TRadj_add = max(0,expParams.CHOICE.TRlength-mod(tmpBlockEnd(run)+expParams.CHOICE.runEndISI-results.CHOICE.timing.scanBlockStartTTL(run),expParams.CHOICE.TRlength));
            tmpBlockEnd(run) = tmpBlockEnd(run)+runEndISI_TRadj_add;
            tmpRunLength = TRlength*ceil(tmpRunLength/TRlength);

            numActualBlockTRs(run) = round(( tmpBlockEnd(run)-results.metaResults1.timing.runStartSecs(run))/TRlength);

        end
        
        %%%% CHECK THESE: Ok, so what that does is get the time to be added
        %%%% to the second block trials based on the number of TRs and TR
        %%%% length in block 1; IMHO that should be fine.
        tmpAllOnsetAdder = [zeros(1,length(find(tmpAllRunNums==1))) TRlength*numActualBlockTRs(1)*ones(1,length(find(tmpAllRunNums==2)))];
        tmpAllOns_eval = tmpAllOnsetAdder + results.metaResults1.timing.progStartSecsRelative;
        
        
        % Basic version modeling all trials as single condition
        allSOTSconds_basic = []; % setting up this basic model
        allSOTSconds_basic(1).name = 'AllTrials'; % which is called 'Alltrials'
        allSOTSconds_basic(1).ons = tmpAllOns_eval(~tmpMissedTrials); % It includes all onsets (for run 1 and 2)
        allSOTSconds_basic(1).dur_Ev = 0;
        allSOTSconds_basic(1).dur_Full = tmpAllRT_eval(~tmpMissedTrials); % the full duration of the trial is the RT
        allSOTSconds_basic(1).ons_Resp = allSOTSconds_basic(1).ons+allSOTSconds_basic(1).dur_Full; % here we model the onset of the responses relative to the beginning of the task (cause that's were they are in time space)
        allSOTSconds_basic(1).P(1).name = 'none'; % we don't have regressors here...
        if any(tmpMissedTrials) % Need to model for eval and exec separately --> this can be omitted in practice, cause no miss trials
            allSOTSconds_basic(2).name = 'MissedTrials';
            allSOTSconds_basic(2).ons = tmpAllOns_eval(tmpMissedTrials);
            allSOTSconds_basic(2).dur_Ev = 0;
            allSOTSconds_basic(2).dur_Full = tmpAllRT_eval(tmpMissedTrials);
            allSOTSconds_basic(2).ons_Resp = allSOTSconds_basic(3).ons+allSOTSconds_basic(3).dur_Full;
            allSOTSconds_basic(2).P(1).name = 'none';
        end
        
        % Parametric setting for choose goal 
        allSOTSconds_pChGoalBin = allSOTSconds_basic; % we copy the previous model and give it a new name
        allSOTSconds_pChGoalBin(1).P(1).name = 'isChooseBest';  % 1 vs. -1 contrast for task
        allSOTSconds_pChGoalBin(1).P(1).P = tmpChooseGoal(~tmpMissedTrials);
        allSOTSconds_pChGoalBin(1).P(1).h = 1; 
        
        allSOTSconds_pRewvGoal = allSOTSconds_basic;
        allSOTSconds_pRewvGoal(1).P(1).name = 'isChooseBest'; 
        allSOTSconds_pRewvGoal(1).P(1).P = tmpChooseGoal(~tmpMissedTrials);
        allSOTSconds_pRewvGoal(1).P(1).h = 1;
        allSOTSconds_pRewvGoal(1).P(2).name = 'OVr';    
        allSOTSconds_pRewvGoal(1).P(2).P = z_AvBid(~tmpMissedTrials);
        allSOTSconds_pRewvGoal(1).P(2).h = 1;
        allSOTSconds_pRewvGoal(1).P(3).name = 'OVg'; 
        allSOTSconds_pRewvGoal(1).P(3).P = z_AvGoalV(~tmpMissedTrials);
        allSOTSconds_pRewvGoal(1).P(3).h = 1;
        allSOTSconds_pRewvGoal(1).P(4).name = 'RVr'; 
        allSOTSconds_pRewvGoal(1).P(4).P = z_ChosvAvRemBid(~tmpMissedTrials);
        allSOTSconds_pRewvGoal(1).P(4).h = 1;
        allSOTSconds_pRewvGoal(1).P(5).name = 'RVg'; 
        allSOTSconds_pRewvGoal(1).P(5).P = z_GoalChosenvAvRem(~tmpMissedTrials);
        allSOTSconds_pRewvGoal(1).P(5).h = 1;

        allSOTSconds_pRewvGoalRT = allSOTSconds_pRewvGoal;
        allSOTSconds_pRewvGoalRT(1).P(6).name = 'RT'; 
        allSOTSconds_pRewvGoalRT(1).P(6).P = z_RT(~tmpMissedTrials);
        allSOTSconds_pRewvGoalRT(1).P(6).h = 1;
        
        
        
        % Parametric setting for choose goal * avBid
        allSOTSconds_pChGoalBin_pAvBid_pChGlAvBd = allSOTSconds_pChGoalBin;
        allSOTSconds_pChGoalBin_pAvBid_pChGlAvBd(1).P(2).name = 'zAvBid';
        allSOTSconds_pChGoalBin_pAvBid_pChGlAvBd(1).P(2).P = AS_nanzscore(allTrAvBid(~tmpMissedTrials))'; % this is actually reduntant. We could just use the pre- z-scored variables...
        allSOTSconds_pChGoalBin_pAvBid_pChGlAvBd(1).P(2).h = 1;
        allSOTSconds_pChGoalBin_pAvBid_pChGlAvBd(1).P(3).name = 'ChGlAvBd_inter';
        allSOTSconds_pChGoalBin_pAvBid_pChGlAvBd(1).P(3).P = allSOTSconds_pChGoalBin_pAvBid_pChGlAvBd(1).P(2).P .* allSOTSconds_pChGoalBin_pAvBid_pChGlAvBd(1).P(2).P;
        allSOTSconds_pChGoalBin_pAvBid_pChGlAvBd(1).P(3).h = 1;
                
        % Alternate version w/ separate best/worst conditions
        allSOTSconds_cBW = [];
        allSOTSconds_cBW(1).name = 'ChooseBest';
        allSOTSconds_cBW(1).ons = tmpAllOns_eval(~tmpMissedTrials & tmpChooseGoal==1);
        allSOTSconds_cBW(1).dur_Ev = 0;
        allSOTSconds_cBW(1).dur_Full = tmpAllRT_eval(~tmpMissedTrials & tmpChooseGoal==1);
        allSOTSconds_cBW(1).ons_Resp = allSOTSconds_cBW(1).ons+allSOTSconds_cBW(1).dur_Full;
        allSOTSconds_cBW(1).P(1).name = 'none';
        allSOTSconds_cBW(2).name = 'ChooseWorst';
        allSOTSconds_cBW(2).ons = tmpAllOns_eval(~tmpMissedTrials & tmpChooseGoal==-1);
        allSOTSconds_cBW(2).dur_Ev = 0;
        allSOTSconds_cBW(2).dur_Full = tmpAllRT_eval(~tmpMissedTrials & tmpChooseGoal==-1);
        allSOTSconds_cBW(2).ons_Resp = allSOTSconds_cBW(2).ons+allSOTSconds_cBW(2).dur_Full;
        allSOTSconds_cBW(2).P(1).name = 'none';
        if any(tmpMissedTrials) 
            allSOTSconds_cBW(3) = allSOTSconds_basic(2);
        end
        
        % make AvBid within B/W model
        allSOTSconds_cBW_pAvBid = allSOTSconds_cBW;
        allSOTSconds_cBW_pAvBid(1).P(1).name = 'zAvBid'; % add 
        allSOTSconds_cBW_pAvBid(1).P(1).P    = z_AvBid(~tmpMissedTrials & tmpChooseGoal==1);
        allSOTSconds_cBW_pAvBid(1).P(1).h = 1;
        allSOTSconds_cBW_pAvBid(2).P(1).name = 'zAvBid';
        allSOTSconds_cBW_pAvBid(2).P(1).P    = z_AvBid(~tmpMissedTrials & tmpChooseGoal==-1);
        allSOTSconds_cBW_pAvBid(2).P(1).h = 1;
        
        
        
        % make VD within B/W model
        allSOTSconds_cBW_pVD = allSOTSconds_cBW;
        allSOTSconds_cBW_pVD(1).P(1).name = 'zGoalvAvRemBid'; % add 
        allSOTSconds_cBW_pVD(1).P(1).P    = z_GoalvAvRemBid(~tmpMissedTrials & tmpChooseGoal==1);
        allSOTSconds_cBW_pVD(1).P(1).h = 1;
        allSOTSconds_cBW_pVD(2).P(1).name = 'zGoalvAvRemBid';
        allSOTSconds_cBW_pVD(2).P(1).P    = z_GoalvAvRemBid(~tmpMissedTrials & tmpChooseGoal==-1);
        allSOTSconds_cBW_pVD(2).P(1).h = 1;
        
        
        %% make liking model but do try catch bc e.g. sub 7004 is missing eval data
        % actually... that's more of an issue in later steps... cause it
        % will just put nans everywhere in this one... I think....
        
        try
        allSOTSconds_cBW_pLike = allSOTSconds_cBW;
        allSOTSconds_cBW_pLike(1).P(1).name = 'zlike'; % add 
        allSOTSconds_cBW_pLike(1).P(1).P    = z_like(~tmpMissedTrials & tmpChooseGoal==1);
        allSOTSconds_cBW_pLike(1).P(1).h = 1;
        allSOTSconds_cBW_pLike(2).P(1).name = 'zlike';
        allSOTSconds_cBW_pLike(2).P(1).P    = z_like(~tmpMissedTrials & tmpChooseGoal==-1);
        allSOTSconds_cBW_pLike(2).P(1).h = 1;
        catch
            fprintf('%d is missing this data', curSubNum);
        end
        
        %% for single trial extraction:
               allSOTSconds_AllTrials = [];
                for trialind = 1:length(tmpAllOns_eval)
                    allSOTSconds_AllTrials(trialind).name = ['Trial',num2str(trialind)];
                    allSOTSconds_AllTrials(trialind).ons = tmpAllOns_eval(trialind);
                    allSOTSconds_AllTrials(trialind).dur_Ev = 0;
                    allSOTSconds_AllTrials(trialind).dur_EvFull = 3.0;
                    allSOTSconds_AllTrials(trialind).dur_Full = tmpAllRT_eval(trialind);
                    allSOTSconds_AllTrials(trialind).ons_Resp = tmpAllOns_eval(trialind) + tmpAllRT_eval(trialind);
                    allSOTSconds_AllTrials(trialind).P(1).name = 'none';
                end
        
        
        
        
        
        %% These are the SOTS structs I've created above:
        tmpVars=whos; % gets all variables in workspace
        tmpVars = {tmpVars(:).name}; % keeps only their names
        tmpVars = tmpVars(find(strncmp(tmpVars,'allSOTSconds',12))); % selects the ones with allSOTSconds in the beginning
        % MEMO: important to have a unique common beginning for this to
        % work
        
        for tvi=1:length(tmpVars) % ok, so now we loop through these
            curSOTSname = tmpVars{tvi};
            curSOTS = eval(curSOTSname); % gotta love evil eval
            
            names = {curSOTS(:).name};
            onsets = {curSOTS(:).ons};
            durations = {curSOTS(:).dur_Full}; 
            pmod = struct('name',{''},'param',{},'poly',{});
            for pmii = 1:length(curSOTS)
                if isfield(curSOTS(pmii).P,'P')
                    pmod(pmii).name = {curSOTS(pmii).P(:).name};
                    for pmsubi = 1:length(curSOTS(pmii).P)
                        curSOTS(pmii).P(pmsubi).P = (curSOTS(pmii).P(pmsubi).P);
                    end
                    pmod(pmii).param = {curSOTS(pmii).P(:).P};
                    pmod(pmii).poly = {curSOTS(pmii).P(:).h};
                end
            end
            
            subAllSOTSnames{tvi} = curSOTSname;
            
            subAllNames{tvi,1} = names;
            subAllOnsets{tvi,1} = onsets;
            subAllDursFull{tvi,1} = durations;
            subAllPmod{tvi,1} = pmod;
            
            durations = {curSOTS(:).dur_Ev};
            subAllDursEv{tvi,1} = durations;
            
            onsets = {curSOTS(:).ons_Resp};
            subAllOnsetsResp{tvi,1} = onsets;
            
            clear('names','onsets','durations','pmod');
        end
        
        %% Saving out all of the sots files
        for tviii = 1:size(subAllPmod,1)            
            curSOTSname = subAllSOTSnames{tviii};
            
            modSubAllPmods = subAllPmod(tviii,:);
            
            names = subAllNames{tviii,1};
            onsets = subAllOnsets{tviii,1};
            durations = subAllDursFull{tviii,1};
            pmod = modSubAllPmods{1};
            
            for oi =1:length(onsets)
                if length(find(onsets{oi}<0))>0
                    display(['CRITICAL!!!!!! BAD ONSET FOR SUBJECT: ',p.subID,', ANALYSIS: ',curSOTSname(13:end)]);
                end
                if length((onsets{oi}))==0
                    display(['MISSING EVENTS FOR SUBJECT: ',p.subID,', ANALYSIS: ',curSOTSname(13:end)]);
                end
            end
            
            for p1 = 1:length(pmod)
                for p2 = 1:length(pmod(p1).param)
                    if length(unique(pmod(p1).param{p2}))<=1
                        display(['NO VARIATION IN PMOD ',num2str(p2),...
                            ' FOR COND ',num2str(p1),', SUBJECT: ',p.subID,', ANALYSIS: ',curSOTSname(13:end)]);
                    end
                end
            end
            
            if ~exist('SOTS','dir')
                unix(['rm BAS*sots*mat']);
                mkdir('SOTS');
            end
            if tviii==1
                unix(['rm -f SOTS/BAS*sots*mat']);
            end
            
            save(['SOTS/BAS',curSOTSname(13:end),'_sots_allTcat.mat'],'names','onsets','durations','pmod');
            
            durations = subAllDursEv{tviii,1};
            save(['SOTS/BAS',curSOTSname(13:end),'_Ev_sots_allTcat.mat'],'names','onsets','durations','pmod');
            
            onsets = subAllOnsetsResp{tviii,1};
            durations = subAllDursEv{tviii,1};
            save(['SOTS/BAS',curSOTSname(13:end),'_Ev_Resp_sots_allTcat.mat'],'names','onsets','durations','pmod');
            
            clear('names','onsets','durations','pmod');
            
            clear modSubAllPmods;
        end
        clear('subAllNames','subAllOnsets','subAllDurs','subAllPmod');
        
        
    else
        %% Excluded subject
        display(['EXCLUDING ',p.subID]);
        try
            blah = [];
            save(fullfile(basepathMR,curSubMR,'EXCLUDE.mat'),'blah');
        end
    end
    
    %% Precautionary cleanup:
    clear tmp*; clear allSOTScond*;    
    clear expParams; clear results;
end

