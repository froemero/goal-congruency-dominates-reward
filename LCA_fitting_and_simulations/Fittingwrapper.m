%% MINIMIZATION

% gets input from MYCOSTFUNCTION

clear; clc; close all;
PATH = '/Users/Romy/Dropbox (Brown)/ShenhavLab/experiments/bas/Analysis/HDDM/';
% load behavioral DatResearch
allSubDataTable=readtable(sprintf('%sBASB_108Data_041118.csv', PATH)); 

before = allSubDataTable;


%% flags for type of model
% Model 1: reward only model
% Model 2: value flip model
% Model 3: threshold flip model
mType = 2;
% flip model
if mType == 2
    isflip = 1;
    
else
    isflip=0;
    
end


%% scale input values max 1
% Product values are sorted by value... 

VALS = [allSubDataTable.allTrProdBid1, allSubDataTable.allTrProdBid2, allSubDataTable.allTrProdBid3, allSubDataTable.allTrProdBid4];


%% get choice in order
ACTCHOICE=[];
for iii = 1:length(VALS)
    tmp=Shuffle([find(VALS(iii,:)==allSubDataTable.ChosenV_FX(iii))]);
   ACTCHOICE(iii) = tmp(1); 
      
end

ACTRT = allSubDataTable.rt;

%% now 
if isflip
    fprintf('I just flipped!\n')
VALS(allSubDataTable.isChooseBestTrial==0,:) = -VALS(allSubDataTable.isChooseBestTrial==0,:)+unique(max(VALS));

end

%%
VALS = (VALS)/unique(max(VALS));
 % worst --> highest
%value first; best --> highest value last

%%

startPars = [1.8 0.5 1 0.1 0.5];
lowerL = [0.5 0.1 0.350 0 0];
upperL = [3 1 2 1 1];

DATA = [VALS, ACTCHOICE', ACTRT];
step=0.05;
realizations=1000;
accumBound=0;
mintype='MLE'; 


 [fitpars,objfunval,exitflag] = ...
        fminsearchbnd(@(mypars) MYCOSTFUNCTION(mypars, DATA, realizations, step, accumBound, mintype),... 
        startPars,lowerL,upperL,optimset('MaxFunEvals',1e6,'MaxIter',1e6));   
