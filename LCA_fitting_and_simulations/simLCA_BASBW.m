%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Goal: simulate choices and RTs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) Simulate LCA just for data as is (reward model)
%% 2) Simulate data with flipping (goal value model)


clear; clc; close all;
PATH = '/Users/Romy/Dropbox (Brown)/ShenhavLab/experiments/bas/Analysis/HDDM/';
% load behavioral DatResearch
allSubDataTable=readtable(sprintf('%sBASB_108Data_041118.csv', PATH)); 

before = allSubDataTable;


%% flags for type of model
% Model 1: reward only model
% Model 2: value flip model

for models = 1:2
    
    
if models ==1
 mType = 1;   
else
mType = 2;
end
% flip model
if mType == 2
    isflip = 1;
    
else
    isflip=0;
    
end


%% Product values are sorted by value... 

VALS = [allSubDataTable.allTrProdBid1, allSubDataTable.allTrProdBid2, allSubDataTable.allTrProdBid3, allSubDataTable.allTrProdBid4];


%% get choice in order
AACTCHOICE=[];
for iii = 1:length(VALS)
    tmp=Shuffle([find(VALS(iii,:)==allSubDataTable.ChosenV_FX(iii))]);
   ACTCHOICE(iii) = tmp(1); 
      
end

ACTRT = allSubDataTable.rt;

%% now flip values in goal value model
if isflip
    fprintf('I just flipped!\n')
VALS(allSubDataTable.isChooseBestTrial==0,:) = -VALS(allSubDataTable.isChooseBestTrial==0,:)+unique(max(VALS));

end

%% Scale input values to max 1
VALS = (VALS)/unique(max(VALS));
VALSBW = [VALS, allSubDataTable.isChooseBestTrial]; % worst --> highest
%value first; best --> highest value last

%% value distributions
figure()
histogram(mean(VALSBW(VALSBW(:,5)==0,1:4),2));
hold on
histogram(mean(VALSBW(VALSBW(:,5)==1,1:4),2));

%%
figure()
ksdensity(std(VALSBW(VALSBW(:,5)==0,:)'), 'Bandwidth', 0.1)
hold on
ksdensity(std(VALSBW(VALSBW(:,5)==1,:)'), 'Bandwidth', 0.1)

%%

[VALS, IA, IC] = unique(VALSBW(:,1:4),'rows');

VALS = VALS(:, 1:4);


%%

doQuickPlot = 1;
nOPT = 4; % number of options to choose from
nTr = length(VALS); % might be usefull to have

SIMCHOICE = nan(length(allSubDataTable.allTrProdBid1),2);
SIMRT     = nan(length(allSubDataTable.allTrProdBid1),1);
SIMERR     = zeros(length(allSubDataTable.allTrProdBid1),1);
%% 

% Simulation parameters:
realizations = 1000; % how many paths to simulate 
stepsize = 0.05; % temporal resolution of process (e.g., proceed in 50ms steps)
% 
% z = threshold
% z = 2.3679; %1.8; % originally 1 from 1.5 % previously 1.2/ 1.3/ 1.4
% s = noise coefficient
% s = 0.6017; %0.5; %0.6; 0.4
% x0 = starting point of drift (between -z and +z; 0 = unbiased)

% t0 = non-decision time
% t0 = 0.6863; %1;  % 350ms
% k = decay rate of integrators
% k = 0.1; % first try was .1
% w = weight of mutual inhibition
% w = 0.5; % first try was .8? 0.4
% set lower-bound of accumulation at 0:
accumBound = 0;
z= 2.243; s = 0.587; t0=0.849; k= 0.153; w = 0.671;

%%
tic;
initprogress = 0;
fprintf('>>Simulation progress: 0 >>');
for trial= 1 : length(VALS(:,1))

    I = [VALS(trial,:)]; 
    Choice=nan(1,realizations);
    Errors=zeros(1,realizations);
    RT=nan(1,realizations);
    Xs=nan(nOPT,10000,realizations); %  to save the traces for all options
    %the way this is set up, the max choice duration is 50 s (1 sample every 5 ms)
    
    % Laurence Hunt's approach (mutinh.m) was different: he has the accumulators
    % separately in different variables rather than in different lines in
    % one variable. That way he updates all trials in those two variables
    % at the same time. That only works for two choices then, but might be
    % translated to something with more options as well... might be ugly
    % though...
    
    
    % this little call is making sure that we use the same random values accross all simulations (at least as long as the length of the accumulation process is the same) 
    rng(1);
    
    %% Goal: vectorize these 
    % run Monte Carlo simulation:
    
        t=0;      % time start
        l=1;      % nth timestepsize
        %x(:,l)=x0;  % path start
        Xs(:,l,:)=0;
        hasConverged = zeros(realizations,1);
        convInts =[];
        stop=0;
        while stop==0 % stop when all iterations have converged or when preset vector length is reached
            if length(I)==2
                % CORE LINE OF CODE:
              %  x(:,l+1) = x(:,l) + stepsize * (I + [-k,-w;-w,-k] * x(:,l)) + sqrt(stepsize) * s * randn(length(I),1); % advances the accumulation process
            else
                
                
                for opti=1:length(I)

                    Xs(opti,l+1,~hasConverged) = Xs(opti,l,~hasConverged) + stepsize * (I(opti) -k * Xs(opti,l,~hasConverged) - sum(w * Xs(setdiff(1:length(I),opti),l,~hasConverged)))...
                        + sqrt(stepsize) * s * randn(size(Xs(opti,l,~hasConverged))); % advances the accumulation process
                end
                
                
                
            end
            
            if accumBound
               Xs(Xs<0)=0;
            end
            
            t=t+stepsize;  % advance time
            
            l=l+1;     % advance timestepsize
            % check 
            convInts = find(any(Xs(:, l,:)>=z,1));
            if ~isempty(convInts)% check if any options hit threshold
                hasConverged(convInts)=1;
                RT(convInts)=t;            % record the RT (i.e., first threshold passage time)
                isMaxAct =   squeeze(Xs(:,l,convInts))==max(squeeze(Xs(:,l,convInts)));
                [choiceNb,~,~]=find(isMaxAct(:,sum(isMaxAct)==1));   % record which option hit threshold
                Choice(convInts(sum(isMaxAct)==1))=choiceNb;
                if any(sum(isMaxAct)~=1)
                Choice(convInts(sum(isMaxAct)~=1))=  randi([1 4],1,sum(sum(isMaxAct)~=1)); % determine choice randomly if there's no winnner
                end
            end
            if (sum(hasConverged) == realizations || l== size(Xs,2))
                stop=1;
                
            end
        end
    
    DT(trial,:) = RT;      % DT = decision process time
    RT = RT + t0; % fixed offset
    % Full and conditional RT distributions
    simDistRTall = RT;
    % when we flip values: worst --> highest value first; best --> highest value last
    % when we don't flip highest value is always last --> this is max
    
    % chosen in reward space...
     if I(1) < I(2)
        Errors(Choice~=4)= 1;
     elseif  I(2) < I(1)
        Errors(Choice~=1)= 1;
     end
    
    %% get most common choice and median RT
    % when we flip values: worst --> highest value first; best --> highest value last
    % when we don't flip highest value is always last --> this is max
    % chosen in reward space...
    % get P(highest chosen and P(lowest chosen) Assign later
    % based on choice type.
     SIMCHOICE(trial,1 )   = length(Choice(Choice==4))/realizations;
     SIMCHOICE(trial,2 )   = length(Choice(Choice==1))/realizations;  
   
    SIMRT(trial)     = nanmedian(RT); 
    SIMERR(trial)    = mean(Errors);

    progress = round(trial/nTr);
    if progress>initprogress
    fprintf('\b\b\b\b%3d%%', progress);
    end
    initprogress=progress;
end

fprintf('\n');
toc
%% replicate the actual numbers of observation per value combination

OVALS= [allSubDataTable.allTrProdBid1, allSubDataTable.allTrProdBid2, allSubDataTable.allTrProdBid3, allSubDataTable.allTrProdBid4];
if isflip
OVALS(allSubDataTable.isChooseBestTrial==0,:) = -OVALS(allSubDataTable.isChooseBestTrial==0,:)+10;
end

OVALS = OVALS/10;


%% now find rows
SIMRTO = SIMRT;
SIMCHOICEO = SIMCHOICE;
SIMERRO = SIMERR;
%%
% go through simulated Values
for nchoices= 1:length(VALS)
    
    SIMRT(ismember(OVALS, VALS(nchoices,:),'rows')) = SIMRTO(nchoices);
    
    SIMCHOICE(ismember(OVALS, VALS(nchoices,:),'rows'),1) = SIMCHOICEO(nchoices,1);
    SIMCHOICE(ismember(OVALS, VALS(nchoices,:),'rows'),2) = SIMCHOICEO(nchoices,2);
    
    SIMERR(ismember(OVALS, VALS(nchoices,:),'rows')) = SIMERRO(nchoices);
end
%%
CPSIMCHOICE= SIMCHOICE; % copy in case I fucked this up.

SIMCHOICE = CPSIMCHOICE(:,1); % P(highest option chosen)

SIMCHOICE(allSubDataTable.isChooseBestTrial==0) = CPSIMCHOICE(allSubDataTable.isChooseBestTrial==0,2); % P(lowest option chosen in choose worst)

%% subset the data table (original values are stored in "before")

SIMallSubDataTable= allSubDataTable;


%%

meanRTMIX = nanmean(SIMRT(SIMallSubDataTable.CondType==1 & allSubDataTable.rt<15));
meanRTLL = nanmean(SIMRT(SIMallSubDataTable.CondType==2 & allSubDataTable.rt<15));
meanRTMM = nanmean(SIMRT(SIMallSubDataTable.CondType==3 & allSubDataTable.rt<15));
meanRTHH = nanmean(SIMRT(SIMallSubDataTable.CondType==4 & allSubDataTable.rt<15));

simRTCond = [meanRTMIX, meanRTLL, meanRTMM, meanRTHH];

omeanRTMIX = nanmean(allSubDataTable.rt(allSubDataTable.CondType==1 & allSubDataTable.rt<15));
omeanRTLL = nanmean(allSubDataTable.rt(allSubDataTable.CondType==2 & allSubDataTable.rt<15));
omeanRTMM = nanmean(allSubDataTable.rt(allSubDataTable.CondType==3 & allSubDataTable.rt<15));
omeanRTHH = nanmean(allSubDataTable.rt(allSubDataTable.CondType==4 & allSubDataTable.rt<15));
RTCond = [omeanRTMIX, omeanRTLL, omeanRTMM, omeanRTHH];

figure();
set(gcf,'Color', [1 1 1])
hold on
subplot(2,1,1)
bar(simRTCond)
title('simulated')
set(gca,'XTickLabel',{'Mixed','Low','Mid','High'});
ylim([0 10])
ylabel('RT')
subplot(2,1,2)
bar(RTCond)
title('empirical')
set(gca,'XTickLabel',{'Mixed','Low','Mid','High'});
ylim([0 10])
ylabel('RT')


%%

meanRTMIXB = nanmean(SIMRT(SIMallSubDataTable.CondType==1 & SIMallSubDataTable.isChooseBestTrial==1 & allSubDataTable.rt<15));
meanRTLLB = nanmean(SIMRT(SIMallSubDataTable.CondType==2 & SIMallSubDataTable.isChooseBestTrial==1 & allSubDataTable.rt<15));
meanRTMMB = nanmean(SIMRT(SIMallSubDataTable.CondType==3 & SIMallSubDataTable.isChooseBestTrial==1 & allSubDataTable.rt<15));
meanRTHHB = nanmean(SIMRT(SIMallSubDataTable.CondType==4 & SIMallSubDataTable.isChooseBestTrial==1 & allSubDataTable.rt<15));

simRTCondB = [meanRTMIXB, meanRTLLB, meanRTMMB, meanRTHHB];

meanRTMIXW = nanmean(SIMRT(SIMallSubDataTable.CondType==1 & SIMallSubDataTable.isChooseBestTrial==0 & allSubDataTable.rt<15));
meanRTLLW = nanmean(SIMRT(SIMallSubDataTable.CondType==2 & SIMallSubDataTable.isChooseBestTrial==0 & allSubDataTable.rt<15));
meanRTMMW = nanmean(SIMRT(SIMallSubDataTable.CondType==3 & SIMallSubDataTable.isChooseBestTrial==0 & allSubDataTable.rt<15));
meanRTHHW = nanmean(SIMRT(SIMallSubDataTable.CondType==4 & SIMallSubDataTable.isChooseBestTrial==0 & allSubDataTable.rt<15));

simRTCondW = [meanRTMIXW, meanRTLLW, meanRTMMW, meanRTHHW];

omeanRTMIXB = nanmean(allSubDataTable.rt(allSubDataTable.CondType==1 & allSubDataTable.isChooseBestTrial==1 & allSubDataTable.rt<15));
omeanRTLLB = nanmean(allSubDataTable.rt(allSubDataTable.CondType==2 & allSubDataTable.isChooseBestTrial==1 & allSubDataTable.rt<15));
omeanRTMMB = nanmean(allSubDataTable.rt(allSubDataTable.CondType==3 & allSubDataTable.isChooseBestTrial==1 & allSubDataTable.rt<15));
omeanRTHHB = nanmean(allSubDataTable.rt(allSubDataTable.CondType==4 & allSubDataTable.isChooseBestTrial==1 & allSubDataTable.rt<15));
RTCondB = [omeanRTMIXB, omeanRTLLB, omeanRTMMB, omeanRTHHB];

omeanRTMIXW = nanmean(allSubDataTable.rt(allSubDataTable.CondType==1 & allSubDataTable.isChooseBestTrial==0 & allSubDataTable.rt<15));
omeanRTLLW = nanmean(allSubDataTable.rt(allSubDataTable.CondType==2 & allSubDataTable.isChooseBestTrial==0 & allSubDataTable.rt<15));
omeanRTMMW = nanmean(allSubDataTable.rt(allSubDataTable.CondType==3 & allSubDataTable.isChooseBestTrial==0 & allSubDataTable.rt<15));
omeanRTHHW = nanmean(allSubDataTable.rt(allSubDataTable.CondType==4 & allSubDataTable.isChooseBestTrial==0 & allSubDataTable.rt<15));
RTCondW = [omeanRTMIXW, omeanRTLLW, omeanRTMMW, omeanRTHHW];

figure();
set(gcf,'Color', [1 1 1])
hold on
subplot(2,2,1)
bar(simRTCondB)
title('simulated Best')
set(gca,'XTickLabel',{'Mixed','Low','Mid','High'});
ylim([0 10])
ylabel('RT')
subplot(2,2,2)
bar(RTCondB)
title('empirical Best')
set(gca,'XTickLabel',{'Mixed','Low','Mid','High'});
ylim([0 10])
ylabel('RT')
subplot(2,2,3)
bar(simRTCondW)
title('simulated Worst')
set(gca,'XTickLabel',{'Mixed','Low','Mid','High'});
ylim([0 10])
ylabel('RT')
subplot(2,2,4)
bar(RTCondW)
title('empirical Worst')
set(gca,'XTickLabel',{'Mixed','Low','Mid','High'});
ylim([0 10])
ylabel('RT')


%%

meanERRMIX = nanmean(SIMCHOICE(SIMallSubDataTable.CondType==1));
meanERRLL = nanmean(SIMCHOICE(SIMallSubDataTable.CondType==2));
meanERRMM = nanmean(SIMCHOICE(SIMallSubDataTable.CondType==3));
meanERRHH = nanmean(SIMCHOICE(SIMallSubDataTable.CondType==4));

simERRCond = [meanERRMIX, meanERRLL, meanERRMM, meanERRHH];

omeanERRMIX = nanmean(allSubDataTable.response(allSubDataTable.CondType==1));
omeanERRLL = nanmean(allSubDataTable.response(allSubDataTable.CondType==2));
omeanERRMM = nanmean(allSubDataTable.response(allSubDataTable.CondType==3));
omeanERRHH = nanmean(allSubDataTable.response(allSubDataTable.CondType==4));
ERRCond = [omeanERRMIX, omeanERRLL, omeanERRMM, omeanERRHH];

figure();
set(gcf,'Color', [1 1 1])
hold on
subplot(2,1,1)
bar(simERRCond)
ylim([0 1])
title('simulated')
set(gca,'XTickLabel',{'Mixed','Low','Mid','High'});
ylabel('P(goal chosen)')
subplot(2,1,2)
bar(ERRCond)
ylim([0 1])
title('empirical')
set(gca,'XTickLabel',{'Mixed','Low','Mid','High'});
ylabel('P(goal chosen)')


%%

meanERRMIXB = nanmean(SIMCHOICE(SIMallSubDataTable.CondType==1& SIMallSubDataTable.isChooseBestTrial==1));
meanERRLLB = nanmean(SIMCHOICE(SIMallSubDataTable.CondType==2& SIMallSubDataTable.isChooseBestTrial==1));
meanERRMMB = nanmean(SIMCHOICE(SIMallSubDataTable.CondType==3& SIMallSubDataTable.isChooseBestTrial==1));
meanERRHHB = nanmean(SIMCHOICE(SIMallSubDataTable.CondType==4& SIMallSubDataTable.isChooseBestTrial==1));

simERRCondB = [meanERRMIXB, meanERRLLB, meanERRMMB, meanERRHHB];

omeanERRMIXB = nanmean(allSubDataTable.response(allSubDataTable.CondType==1& allSubDataTable.isChooseBestTrial==1));
omeanERRLLB = nanmean(allSubDataTable.response(allSubDataTable.CondType==2& allSubDataTable.isChooseBestTrial==1));
omeanERRMMB = nanmean(allSubDataTable.response(allSubDataTable.CondType==3& allSubDataTable.isChooseBestTrial==1));
omeanERRHHB = nanmean(allSubDataTable.response(allSubDataTable.CondType==4& allSubDataTable.isChooseBestTrial==1));
ERRCondB = [omeanERRMIXB, omeanERRLLB, omeanERRMMB, omeanERRHHB];

meanERRMIXW = nanmean(SIMCHOICE(SIMallSubDataTable.CondType==1& SIMallSubDataTable.isChooseBestTrial==0));
meanERRLLW = nanmean(SIMCHOICE(SIMallSubDataTable.CondType==2& SIMallSubDataTable.isChooseBestTrial==0));
meanERRMMW = nanmean(SIMCHOICE(SIMallSubDataTable.CondType==3& SIMallSubDataTable.isChooseBestTrial==0));
meanERRHHW = nanmean(SIMCHOICE(SIMallSubDataTable.CondType==4& SIMallSubDataTable.isChooseBestTrial==0));

simERRCondW = [meanERRMIXW, meanERRLLW, meanERRMMW, meanERRHHW];

omeanERRMIXW = nanmean(allSubDataTable.response(allSubDataTable.CondType==1& allSubDataTable.isChooseBestTrial==0));
omeanERRLLW = nanmean(allSubDataTable.response(allSubDataTable.CondType==2& allSubDataTable.isChooseBestTrial==0));
omeanERRMMW = nanmean(allSubDataTable.response(allSubDataTable.CondType==3& allSubDataTable.isChooseBestTrial==0));
omeanERRHHW = nanmean(allSubDataTable.response(allSubDataTable.CondType==4& allSubDataTable.isChooseBestTrial==0));
ERRCondW = [omeanERRMIXW, omeanERRLLW, omeanERRMMW, omeanERRHHW];

figure();
set(gcf,'Color', [1 1 1])
hold on
subplot(2,2,1)
bar(simERRCondB)
ylim([0 1])
title('simulated Best')
set(gca,'XTickLabel',{'Mixed','Low','Mid','High'});
ylabel('P(goal chosen)')
subplot(2,2,2)
bar(ERRCondB)
ylim([0 1])
title('empirical Best')
set(gca,'XTickLabel',{'Mixed','Low','Mid','High'});
ylabel('P(goal chosen)')
subplot(2,2,3)
bar(simERRCondW)
ylim([0 1])
title('simulated Worst')
set(gca,'XTickLabel',{'Mixed','Low','Mid','High'});
ylabel('P(goal chosen)')
subplot(2,2,4)
bar(ERRCondW)
ylim([0 1])
title('empirical Worst')
set(gca,'XTickLabel',{'Mixed','Low','Mid','High'});
ylabel('P(goal chosen)')
%%
allSubDataTable.SIMRT= SIMRT;
allSubDataTable.SIMERR = SIMERR;
allSubDataTable.SIMCHOICE = SIMCHOICE;
%%
if mType == 2
writetable(allSubDataTable, sprintf('%sallSubDataTableSIMBASBWflip_01152019.xls', PATH));

elseif mType == 1
writetable(allSubDataTable, sprintf('%sallSubDataTableSIMBASBWreward_01152019.xls', PATH));
 
    
else
    writetable(allSubDataTable, sprintf('%sallSubDataTableSIMBASBWthresh_01152019.xls', PATH));

end
end