function [SIMERR,SIMRT, SIMCHOICE, ALLSIMRT, ALLSIMChoice, didquit ] = sim_LCA_bound_par_pacc(VALS, realizations, stepsize, z, s, t0, k, w, accumBound)



%% required input:
% realizations % how many paths to simulate
% % stepsize 
% z 
% % s = noise coefficient
%

% t0 = round(t0, 3);
% z  = round(z,3);
% s  = round(s,3);
% w  = round(w,3);
% k  = round(k,3);
% % % round all parameter values!

didquit =0;

nOPT = 4; % number of options to choose from

nTr = length(VALS); % might be usefull to have

ALLSIMRT=nan(length(VALS),realizations); % saves all RTs for all trials and sims...
ALLSIMChoice=nan(length(VALS),realizations);


SIMCHOICE = nan(length(VALS),1);
SIMRT     = nan(length(VALS),1);
SIMERR     = zeros(length(VALS),1);


parfor trial= 1 : length(VALS)
    
    % LCA parameters:

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
    
    %% run Monte Carlo simulation:
    t=0;      % time start
    l=1;      % nth timestepsize
    Xs(:,l,:)=0;
    hasConverged = zeros(realizations,1);
    convInts =[];
    stop=0;
    while stop==0 % stop when all iterations have converged or when preset vector length is reached
        if length(I)==2
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
        convInts = find(any(Xs(:, l,:)>=z,1));% check if any options hit threshold
        if ~isempty(convInts)
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
    

    RT = RT + t0; % fixed offset
    % Full and conditional RT distributions
    simDistRTall = RT;
    
    if I(1) < I(2)
        Errors(Choice~=4)= 1;
    elseif  I(2) < I(1)
        Errors(Choice~=1)= 1;
    end
    
    
    %% for each trial get most frequent choice with mode
    
    %% get most common choice and median RT
    % when we flip values: worst --> highest value first; best --> highest value last
    % when we don't flip highest value is always last --> this is max
    % chosen in reward space...
    if I(1) < I(2)
        SIMCHOICE(trial) = length(Choice(Choice==4))/realizations;
    elseif  I(2) < I(1)
        SIMCHOICE(trial) = length(Choice(Choice==1))/realizations;
    elseif length(find(I==max(I)))==4
        SIMCHOICE(trial) = length(Choice(Choice==1))/realizations;
    end
    
    SIMRT(trial)     = nanmedian(RT); %% test what median gives us instead
    SIMERR(trial)    = mean(Errors);
    
    
    %% save choice and RT distributions
    
    ALLSIMRT(trial,:)=RT; % saves all RTs for all trials and sims...
    ALLSIMChoice(trial,:)=Choice;
    
    
end

fprintf('\n');


end

