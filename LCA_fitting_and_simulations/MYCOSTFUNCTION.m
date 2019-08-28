function [MYCOST] = MYCOSTFUNCTION(MYFITPARS, DATA, realizations, step, accumBound, mintype, varargin)
% runs simulation of LCA through sim_LCA and computes cost function

%  INPUT:
%  1) MYFITPARS: vector of parameter values ordered -->
%             z (threshold), s(noise),t0 (non-decision time), k (leakage),
%             w (mutual inhibition weight) --> not all parameters must be
%             there, but those to be fitted must be ordered
%  2) DATA: matrix with 6 columns: Values ordered by reward value, actual choice,
%             actual RT
%  3) realizations: number of runs
%  4) step: sampling interval in s
%  5) accumbound: set negative values to zero? (1 = yes)
%  6) optional: flag method
%  7) fixed parameters (must be entered with their respective names)

%  OUTPUT:
%  MYCOST

try
for i=1:length(varargin)
    
    pnames{i,1}=inputname(6+i);
    
end
        
for parname = 1: length(pnames)
    
    expression = [pnames{parname} '=' num2str(varargin{parname}) ';'] ;
    eval(expression);
    
end
paramnames= {'z'; 's'; 't0'; 'k'; 'w'};

paramstofit=find(~ismember(paramnames, pnames));
catch
   paramnames= {'z'; 's'; 't0'; 'k'; 'w'}; 
   paramstofit= [1:5]; 
    
end
for fitparams= 1: length(paramstofit)
    expression=[paramnames{paramstofit(fitparams)} '=' num2str(MYFITPARS(fitparams)) ';'];
    eval(expression);
    
end


fprintf('\n   >> current parameters: z= %.3f, s = %.3f, t0=%.3f k= %.3f, w = %.3f\n', z, s, t0, k, w)


VALS= unique(DATA(:,1:4),'rows');

%%


tic; 
[SIMERR,SIMRT, SIMCHOICE, ALLSIMRT, ALLSIMChoice] = sim_LCA_bound_par_pacc(VALS, realizations, step, z, s, k, w, t0, accumBound);
toc


OVALS= DATA(:,1:4);




ACTCHOICE = DATA(:,5);
% make responses 0 and 1 for comparison with likelihood
if max(ACTCHOICE) >1
    ACTCHOICEp = ACTCHOICE - 1;
else
    ACTCHOICEp = ACTCHOICE;
end

ACTRT = DATA(:,6);


%% COST FUNCTION
switch mintype
    
    case 'LSE'
        
        % Least squares (LSE)
        
        MYCOST = sum(((ACTRT-SIMRT)).^2);

  
        fprintf('current cost: %d, current mean observed RT: %d, current mean simulated RT: %d \n >> \n', MYCOST, nanmean(ACTRT), nanmean(SIMRT))
  
    case 'MLE'

        
fprintf('starting likelihood estimation...\n')
        for n=1: size(ACTRT,1) % I need to go through the data not the simulations
            
            % now I need to find the relevant simulation:
            [tf, index]=ismember(OVALS(n,:), VALS,'rows');
            
            val= ACTRT(n);
            
            % get choice likelihood by choice 
            AllTSimChoice = ALLSIMChoice(index, :);
            Probr(1)= length(AllTSimChoice(AllTSimChoice==1))/length(AllTSimChoice);
            Probr(2)= length(AllTSimChoice(AllTSimChoice==2))/length(AllTSimChoice);
            Probr(3)= length(AllTSimChoice(AllTSimChoice==3))/length(AllTSimChoice);
            Probr(4)= length(AllTSimChoice(AllTSimChoice==4))/length(AllTSimChoice);

            %%
           
                
                if Probr(ACTCHOICE(n))>0
                    try
                 tmpALLSIMRTl= ALLSIMRT(n,ALLSIMChoice(index,:)==ACTCHOICE(n));
                  [likl,centersl] = ksdensity(tmpALLSIMRTl, 'Support','positive');

                liklnew = likl*(Probr(ACTCHOICE(n)));

                xi=[min(centersl):step:(max(centersl))]';
                yi= interp1q( centersl', liklnew', xi);
                
                [ ~, ix ] = min( abs( xi-val ) );
                ProbR(n)= yi(ix)+0.01;
                catch
                    ProbR(n) = 0.01;
                    end
            
                else
                fprintf('outlier trial detected\n')    
                ProbR(n) = 0.01;   
                    
                end
            if ProbR(n)==0
                
               ProbR(n) = 1e-16; 
                
            end
            
        end
        
        MYCOST= -sum(log(ProbR)); 
        
        fprintf('current cost: %.3f, current mean observed RT: %.3f, current mean simulated RT: %.3f <<\n\n', MYCOST, nanmean(ACTRT), nanmean(SIMRT))

        
    case 'chi2'
        % chi2                 % test on counts per quantile (6 quantiles, see Ratcliff)

end


end

