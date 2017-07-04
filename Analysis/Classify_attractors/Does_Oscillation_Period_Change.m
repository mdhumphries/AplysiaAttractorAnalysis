% cyclical attractor dynamics implies that the magnitude of the oscillation
% changes, but the period does not
% Therefore, the visual changes in burst size are *not* due to slowing of
% the oscillation, but instead due to reduction in magnitude of the
% oscillation, as realised by a spiking element!
%
% Here: test this. Check that the oscillation period does not decrease,
% once the initial transient is accounted for
%
% Updates:
% June 2016: added Pearson's correlation to check magnitude of changes
% 5/7/16: checking mean recurrence times (of all recurrent points)
% 6/7/16: check mean recurrence times (in main Period); add exclusion of
%           windows outside definite attractor
% 6/7/16: added WeightedCorrelation measure
% 11/7/16: fixed Inf SD bug; added weights combining SD of delays and density of recurrence points 
% Mark Humphries 6/7/16

% get eigenvalues of recurrence points
clear all; close all
addpath ../../Functions/

% pick which threshold to look at
oscpars.ixT = 3; % threshold of 10%
oscpars.pDensity = 0.9; 

fname_Pre = '';
load da01_RecurrencePointsData_FilterPts RcrPtData
fnames = {'da01','da02','da03'};

nfiles = numel(RcrPtData);

load([fname_Pre 'JumpAnalysis2']);
load([fname_Pre 'RecurrenceStats_FilterPts']);

oscpars.blnCoalesce = 1;  % only look at change after coalescence
oscpars.blnOff = 1; % only look at change before final fall from attractor
% oscpars.blnStable = 1; % only look at change for all stable windows [in which case, weighted correlation is not needed]
blnSpearman = 1;  % use Spearman's rank in summary

alpha = 0.01;
bins = -1:0.1:1;

%% analyse
for iPrep = 1:nfiles
   iPrep
   for iProg = 1:numel(fnames) 
       % load 
       load([fname_Pre fnames{iProg} '_RecurrencePointsData_FilterPts'],'RcrPtData','pars');
       load(['../' fname_Pre fnames{iProg} '_DataProperties_FunctionAndWindowSize'],'Data','FileTable')
       
       nPoints = [RcrPtData(iPrep).Orbit(oscpars.ixT).Period(:).nPoints];
       ixPrd = find(nPoints == max(nPoints));  
       
%        % raw eigenvalues
%        usedPts = RcrPtData(iPrep).Orbit(statpars.ixT).Period(ixPrd).Npts >= pars.minNeighbours;  % min neighbours from the original analysis file
%        Egs = imag(RcrPtData(iPrep).Orbit(statpars.ixT).Period(ixPrd).Emax(usedPts));   % imaginary part of max eigenvalue
%        Rots = RcrPtData(iPrep).Orbit(statpars.ixT).Period(ixPrd).RotPeriod(usedPts);   % corresponding local estimate of rotation period
%        Ts = RcrPtData(iPrep).Orbit(statpars.ixT).Period(ixPrd).IDs(usedPts);    % indices of recurrence points in this period
%        actualTime = Data(iPrep).stimbins(RcrPtData(iPrep).ix0 + Ts);   % conversion to time
%     
%        figure
%        plot(actualTime,Egs,'ko')
       
       PeriodChange(iPrep,iProg).blnJump = any(Summary.JumpList(:,1) == iPrep & Summary.JumpList(:,2) == iProg);
       
       % get recurrence density
       tsDensity = RcrPtData(iPrep).Orbit(oscpars.ixT).AllRecurrentPoint.pRecur;

       % get time windows
       ixW = RcrPtData(iPrep).winmids;
       tsW = Data(iPrep).stimbins(ixW);
       
       % find admissible windows
       ixwins = 1:numel(ixW);
       ixWinStrt = 1; ixWinEnd = numel(ixwins);
       if oscpars.blnCoalesce  % then only check after coalescence
            ixC = find(Data(iPrep).stimbins == RcrStats(iPrep,iProg).firstW);
            ixWinStrt = find(ixW == ixC);
       end
       if oscpars.blnOff
            ixWinEnd = find(tsDensity > oscpars.pDensity,1,'last');  % final window with P% density of recurrence
       end
       ixwins = ixWinStrt:ixWinEnd; 
       PeriodChange(iPrep,iProg).ixwins = ixwins;
       PeriodChange(iPrep,iProg).Nwins = numel(ixwins);
        
       %% moving window average of eigenvalues estimates from main Period
       RotW = RcrPtData(iPrep).Orbit(oscpars.ixT).Period(ixPrd).mRotPeriod;
       % Spearman's rank to measure for monotonic change
       [PeriodChange(iPrep,iProg).rho,PeriodChange(iPrep,iProg).pval] = corr(tsW(ixwins)',RotW(ixwins)','type','Spearman','rows','complete');
       
       % Pearson to get magnitude
       [PeriodChange(iPrep,iProg).R,PeriodChange(iPrep,iProg).pvalR] = corr(tsW(ixwins)',RotW(ixwins)','type','Pearson','rows','complete');
       
       %% recurrence times of main Period
        DelayW = RcrPtData(iPrep).Orbit(oscpars.ixT).Period(ixPrd).mDelay;
         
         
       % Spearman's rank to measure for monotonic change
       [PeriodChange(iPrep,iProg).rhoRecur,PeriodChange(iPrep,iProg).pvalRecur] = corr(tsW(ixwins)',DelayW(ixwins)','type','Spearman','rows','complete');

       % Pearson to get magnitude
       [PeriodChange(iPrep,iProg).RRecur,PeriodChange(iPrep,iProg).pvalRRecur] = corr(tsW(ixwins)',DelayW(ixwins)','type','Pearson','rows','complete');
      
       %% all recurrence delays
       DelayWall = RcrPtData(iPrep).Orbit(oscpars.ixT).AllRecurrentPoint.mDelay;
       tsDelaySD = RcrPtData(iPrep).Orbit(oscpars.ixT).AllRecurrentPoint.sdDelay;
       WDelaySD = 1./tsDelaySD;
       WDelaySD(isnan(WDelaySD)) = 0; % set nans to weights of 0: not included
       WDelaySD(isinf(WDelaySD)) = 0; % set INFs to weights of 0: perfect = no noise = suspicious...
       
       
       % combined weights?
       WDelaySD = WDelaySD .* tsDensity;
       
       % Spearman's rank to measure for monotonic change
       [PeriodChange(iPrep,iProg).rhoRecurAll,PeriodChange(iPrep,iProg).pvalRecurAll] = corr(tsW(ixwins)',DelayWall(ixwins)','type','Spearman','rows','complete');
       
       % weighted by density
       % [PeriodChange(iPrep,iProg).WrhoRecurAll,PeriodChange(iPrep,iProg).WpvalRecurAll] = WeightedCorrelation(tsW(ixwins)',DelayWall(ixwins)',tsDensity(ixwins)','Spearman');
       
       % weighted by variation
       [PeriodChange(iPrep,iProg).WrhoRecurAll,PeriodChange(iPrep,iProg).WpvalRecurAll] = WeightedCorrelation(tsW(ixwins)',DelayWall(ixwins)',WDelaySD(ixwins)','Spearman');
        
       if isnan(PeriodChange(iPrep,iProg).WrhoRecurAll)
            keyboard
       end
       % Pearson to get magnitude
       [PeriodChange(iPrep,iProg).RRecurAll,PeriodChange(iPrep,iProg).pvalRRecurAll] = corr(tsW(ixwins)',DelayWall(ixwins)','type','Pearson','rows','complete');

       % weighted by density
       % [PeriodChange(iPrep,iProg).WRRecurAll,PeriodChange(iPrep,iProg).WpvalRRecurAll] = WeightedCorrelation(tsW(ixwins)',DelayWall(ixwins)',tsDensity(ixwins)','Pearson');
       % weighted by variation
       [PeriodChange(iPrep,iProg).WRRecurAll,PeriodChange(iPrep,iProg).WpvalRRecurAll] = WeightedCorrelation(tsW(ixwins)',DelayWall(ixwins)',WDelaySD(ixwins)','Pearson');

%        figure
%        plot(tsW,RotW,'ko')
       

       % fit models
%        tbl = table(tsW',RotW','VariableNames',{'Time','Period'});
%        mdlLinear = fitlm(tbl,'linear');
%        
%        % nonlinear models - none of these work!!!
%        modelfun = @(b,x)b(1) + b(2)*(1-exp(-b(3)*x));  
%        beta0 = [1 1 1];
%        mdlExpMax = fitnlm(tbl,modelfun,beta0);
%        
%        modelfun = @(b,x)b(1) + b(2)*exp(-x/b(3));  % exponential decay
%        beta0 = [1 1 1];
%        mdlExpDecay = fitnlm(tbl,modelfun,beta0);
% 
%        modelfun = @(b,x)b(1) + b(2)*exp(-x/b(3)).^b(4);  % stretched exponential decay
%        beta0 = [1 1 1 1];
%        mdlExpSDecay = fitnlm(tbl,modelfun,beta0);
%       
%        % plot
%        yLin = predict(mdlLinear,tsW');
%        yExpMax = predict(mdlExpMax,tsW');
%        yExpDecay = predict(mdlExpDecay,tsW');
%        yExpSDecay = predict(mdlExpSDecay,tsW');
%        
       figure; 
       subplot(211),plot(tsW(ixwins),RotW(ixwins),'ko'); hold on
       plot(tsW,DelayW,'o','Color',[0.3 0.3 0.7]);
%        plot(tsW,yLin,'k')
%        plot(tsW,yExpMax,'r')
%        plot(tsW,yExpDecay,'b')
%        plot(tsW,yExpSDecay,'b-')
       title(['Prep ' num2str(iPrep) '; Program ' num2str(iProg)])
       xlabel('Time (s)')
       ylabel('Mean rotation period (s)')
       
       subplot(212),plot(tsW(ixwins),DelayWall(ixwins),'ko'); hold on
%        plot(tsW,yLin,'k')
%        plot(tsW,yExpMax,'r')
%        plot(tsW,yExpDecay,'b')
%        plot(tsW,yExpSDecay,'b-')
       title(['Prep ' num2str(iPrep) '; Program ' num2str(iProg)])
       xlabel('Time (s)')
       ylabel('Mean recurrence time (All) (s)')
   end
end

%% plot stuff
% h = histc([PeriodChange.rho],bins);
% figure
% bar(bins,h,'histc');

%% rotation period
if blnSpearman
    % Spearman's
    r = [PeriodChange.rho]; p =[PeriodChange.pval];
else
    % Pearons's
    r = [PeriodChange.R]; p = [PeriodChange.pvalR];
end

figure
plot(r(p > alpha), p(p>alpha),'o','Color',[0.6 0.6 0.6],'MarkerSize',5,'MarkerFaceColor',[0.6 0.6 0.6]); hold on
plot(r(p <= alpha), p(p<=alpha),'ro','MarkerSize',5,'MarkerFaceColor','r')
title('Rotation period p-values');
xlabel('Correlation'); ylabel('P value')

hNon = histc(r(p > alpha),bins);
hYes = histc(r(p <= alpha),bins);
figure
hno = bar(bins,hNon,'histc'); hold on
hyes = bar(bins,hYes,'histc'); hold on
%shading flat
set(hno,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[1 1 1]);
set(hyes,'FaceColor',[0.8 0.2 0.2],'EdgeColor',[1 1 1]);
line([-1 1],[0 0],'Color',[0 0 0])
axis([-1 1 0 6])
xlabel('Correlation \rho')
ylabel('No. programs')

% joint scatter
matR = reshape(r,nfiles,numel(fnames));
figure
line([0.5 3.5],[0 0],'Color',[0 0 0],'LineWidth',0.5); hold on
plot([1 2 3],matR','.-','MarkerSize',5); 
set(gca,'XLim',[0.5 3.5]);
set(gca,'XTick',[1 2 3])
xlabel('Stimulations')
ylabel('Correlation \rho')
title('Rotation period program effect');

%% Main Period recurrence time 
if blnSpearman
    % Spearman's
    r = [PeriodChange.rhoRecur]; p =[PeriodChange.pvalRecur];
else
    % Pearons's
    r = [PeriodChange.RRecur]; p = [PeriodChange.pvalRRecur];
end

figure
plot(r(p > alpha), p(p>alpha),'o','Color',[0.6 0.6 0.6],'MarkerSize',5,'MarkerFaceColor',[0.6 0.6 0.6]); hold on
plot(r(p <= alpha), p(p<=alpha),'ro','MarkerSize',5,'MarkerFaceColor','r')
title('Recurrence time (Period) p-values');
xlabel('Correlation'); ylabel('P value')

hNon = histc(r(p > alpha),bins);
hYes = histc(r(p <= alpha),bins);
figure
hno = bar(bins,hNon,'histc'); hold on
hyes = bar(bins,hYes,'histc'); hold on
%shading flat
set(hno,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[1 1 1]);
set(hyes,'FaceColor',[0.8 0.2 0.2],'EdgeColor',[1 1 1]);
line([-1 1],[0 0],'Color',[0 0 0])
axis([-1 1 0 6])
xlabel('Correlation \rho')
ylabel('No. programs')
title('Recurrence time (Period) distribution');

% joint scatter
matR = reshape(r,nfiles,numel(fnames));
figure
line([0.5 3.5],[0 0],'Color',[0 0 0],'LineWidth',0.5); hold on
plot([1 2 3],matR','.-','MarkerSize',5)
set(gca,'XLim',[0.5 3.5]);
set(gca,'XTick',[1 2 3])
xlabel('Stimulations')
ylabel('Correlation \rho')
title('Recurrence time (Period) program effect');

%% Main Period recurrence time, weighted 

if blnSpearman
    % Spearman's
    % r = [PeriodChange.rhoRecurAll]; p =[PeriodChange.pvalRecurAll];
    r = [PeriodChange.WrhoRecurAll]; p =[PeriodChange.WpvalRecurAll];
    % look only at programs without Divergent periods ("jumps")
%     blnJumps = [PeriodChange.blnJump];
%     r(blnJumps) = nan; p(blnJumps) = nan;
else
    % Pearons's
    r = [PeriodChange.RRecurAll]; p = [PeriodChange.pvalRRecurAll];
    % r = [PeriodChange.WRRecurAll]; p =[PeriodChange.WpvalRRecurAll];

end


figure
plot(r(p > alpha), p(p>alpha),'o','Color',[0.6 0.6 0.6],'MarkerSize',5,'MarkerFaceColor',[0.6 0.6 0.6]); hold on
plot(r(p <= alpha), p(p<=alpha),'ro','MarkerSize',5,'MarkerFaceColor','r')
title('Recurrence time (All) p-values');
xlabel('Correlation'); ylabel('P value')

hNon = histc(r(p > alpha),bins);
hYes = histc(r(p <= alpha),bins);
figure
hno = bar(bins,hNon,'histc'); hold on
hyes = bar(bins,hYes,'histc'); hold on
%shading flat
set(hno,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[1 1 1]);
set(hyes,'FaceColor',[0.8 0.2 0.2],'EdgeColor',[1 1 1]);
line([-1 1],[0 0],'Color',[0 0 0])
axis([-1 1 0 6])
xlabel('Correlation \rho')
ylabel('No. programs')
title('Recurrence time (All) distribution');

% joint scatter
matR = reshape(r,nfiles,numel(fnames));
figure
line([0.5 3.5],[0 0],'Color',[0 0 0],'LineWidth',0.5); hold on
plot([1 2 3],matR','.-','MarkerSize',5)
set(gca,'XLim',[0.5 3.5]);
set(gca,'XTick',[1 2 3])
xlabel('Stimulations')
ylabel('Correlation \rho')
title('Recurrence time (All) program effect');

% EFFECTS here?
matP = reshape(p,nfiles,numel(fnames));

Rsize = mean(abs(matR));
Peffect = sum(matP < 0.01);


%% save
save([fname_Pre 'PeriodChange'],'PeriodChange','oscpars')
