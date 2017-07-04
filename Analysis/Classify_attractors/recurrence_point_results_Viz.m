%%% look at recurrence point time-series; and get statistics over all
%%% programs
% statspars.X:
%   struct of analysis parameters           
%               ixT: index of recurrence threshold results to use throughout
%          pDensity: min recurrent point density to count as "stably on attractor" 
%          Tperturb: threshold of point density that indicates non-recurrent period
%
% RcrStats 
%   for each preparation (iR), for each program (iP), for the Period with the maximum number of recurrence points:
%          nUsedPts: number of recurrence points 
%              mRot:
%             sdRot:
%             mEigR:
%            sdEigR:
%             mEigI: 
%            sdEigI:
%        prctStable:
%        densPeriod:
%     densAllPeriod:
%       densNoRecur:
%           nPeriod:
%            firstW:  middle of first window that meets stability criterion (sec. since start of recording)
%         tsPerturb:  time of each potential jump (sec. since start of recording): once fallen below perturbation criterion threshold, the middle of the window with the lowest density of recurrence points
%             tsOff:  final perturbation: [threshold crossing; minimum recurrence] (sec. since start of the recording): middle of the windows with these properties, after which stability is *not* obtained  
%         tsAllRots:

clear all; close all
addpath ../../Functions/

% pick which threshold to look at
statpars.ixT = 3; % threshold of 10%
statpars.pDensity = 90; % min point density to count as "coalescence time"
statpars.Tperturb = 50; % falling below threshold of point density indicates non-recurrent period
% nPerturb = 2;  % min number of consecutive windows below threshold

load da01_RecurrencePointsData_FilterPts RcrPtData
fname_Pre = '';
fnames = {'da01','da02','da03'};

nfiles = numel(RcrPtData);

PerColor = [0.4 0.4 0.4]; % dominant period
AllColor = [1 0 0];  % all recurrence points 
NoPerColor = [0.6 0.6 0.9];  % no period at all, but recurrent
%% visualise...
for iR = 1:nfiles
   for iF = 1:numel(fnames) 
       % load 
       load([fname_Pre fnames{iF} '_RecurrencePointsData_FilterPts'],'RcrPtData','pars');
       load(['../' fname_Pre  fnames{iF} '_DataProperties_FunctionAndWindowSize'],'Data','FileTable')

       %% get summary stats for dominant period
       nPoints = [RcrPtData(iR).Orbit(statpars.ixT).Period(:).nPoints];
       ixPrd = find(nPoints == max(nPoints));  
       
       usedPts = RcrPtData(iR).Orbit(statpars.ixT).Period(ixPrd).Npts >= pars.minNeighbours;
       RcrStats(iR,iF).nUsedPts = numel(usedPts); % points used to compute each mean etc!
       RcrStats(iR,iF).mRot =  mean(RcrPtData(iR).Orbit(statpars.ixT).Period(ixPrd).RotPeriod(usedPts));
       RcrStats(iR,iF).sdRot =  std(RcrPtData(iR).Orbit(statpars.ixT).Period(ixPrd).RotPeriod(usedPts));
       RcrStats(iR,iF).mEigR = mean(real(RcrPtData(iR).Orbit(statpars.ixT).Period(ixPrd).Emax(usedPts)));
       RcrStats(iR,iF).sdEigR = std(real(RcrPtData(iR).Orbit(statpars.ixT).Period(ixPrd).Emax(usedPts)));
       RcrStats(iR,iF).mEigI = mean(imag(RcrPtData(iR).Orbit(statpars.ixT).Period(ixPrd).Emax(usedPts)));
       RcrStats(iR,iF).sdEigI = std(imag(RcrPtData(iR).Orbit(statpars.ixT).Period(ixPrd).Emax(usedPts)));
       
       RcrStats(iR,iF).prctStable = sum(imag(RcrPtData(iR).Orbit(statpars.ixT).Period(ixPrd).Emax) < 1e-6) ./ numel(RcrPtData(iR).Tdelay(:,1));   
      
       % change over time
%        figure
%        plot([RcrPtData(iR).Orbit(statpars.ixT).Period(ixPrd).IDs(usedPts)],[RcrPtData(iR).Orbit(statpars.ixT).Period(:).RotPeriod(usedPts)],'k.')
%        plot([RcrPtData(iR).Orbit(statpars.ixT).Period(ixPrd).IDs(usedPts)],[RcrPtData(iR).Orbit(statpars.ixT).Period(:).Emax(usedPts)],'k.')
      
        % stats of density of period and other period-related properties
       RcrStats(iR,iF).densPeriod = nPoints ./numel(RcrPtData(iR).Tdelay(:,1));  % each set of points > tmin
       RcrStats(iR,iF).densAllPeriod = sum(nPoints)/numel(RcrPtData(iR).Tdelay(:,1));  % all points > tmin
       RcrStats(iR,iF).densNoRecur = RcrPtData(iR).Orbit(statpars.ixT).NoRecurrence.nPoints ./numel(RcrPtData(iR).Tdelay(:,1));
       RcrStats(iR,iF).nPeriod = RcrPtData(iR).Orbit(statpars.ixT).nPeriods;
       RcrStats(iR,iF).propInMaxPeriod = nPoints(ixPrd) ./ (sum(nPoints) + RcrPtData(iR).Orbit(statpars.ixT).NoPeriod.nPoints);  % proportion of all recurrent points in main Period
       RcrStats(iR,iF).propComparePeriods = nPoints ./ sum(nPoints);  % proportion of periodic points in main Period
       
       % get first window with density *any* recurrent points > X 
       PrctWnone = 100*RcrPtData(iR).Orbit(statpars.ixT).NoRecurrence.nWinPoints/RcrPtData(iR).winDT;  % 90% of these
       PrctWPeriod = 100*RcrPtData(iR).Orbit(statpars.ixT).Period(ixPrd).nWinPoints/RcrPtData(iR).winDT;
       
       % keyboard
       % PrctWanyperiod = (RcrPtData(iR).Orbit(statpars.ixT).NoRecurrence.nWinPoints + RcrPtData(iR).Orbit(statpars.ixT).NoPeriod.nWinPoints)/RcrPtData(iR).winDT; % 75% of these
       ixFirst = find(PrctWnone < 100-statpars.pDensity,1,'first');
       if isempty(ixFirst) % never reaches attractor
            tsFirst = nan; 
            tsPerturb = nan;
       else
            % store coalscence point
            tsFirst = Data(iR).stimbins(RcrPtData(iR).winmids(ixFirst));
            % find possible perturbations and transitions: when recurrence
            % density falls and rises again...
            
            % for each possible perturbation...
            % check it is sustained (> N windows)
            % check for local minima in entire period between first fall
            % and when it recover to program (90% density)... (key to
            % finding perturbation, not stopping!)
            
            allixs = 1:numel(PrctWnone(ixFirst:end));
            ixs = ixFirst-1 + find(PrctWnone(ixFirst:end) > 100-statpars.Tperturb);  % points below perturbation threshold
            
            % search each possible perturbation for return to attractor
            tsPerturb = []; tsOff = nan;
            if ~isempty(ixs)  % there are potential perturbations
                allstable = find(PrctWnone < 100-statpars.pDensity);  % points that meet stability criteria
                blnCont = 1; 
                while blnCont
                    nextStable = allstable(find(allstable > ixs(1),1));
                    if ~isempty(nextStable)  % then returns to attractor state (not e.g. end of recording)
                        temp = find(PrctWnone(ixs(1):nextStable) == max(PrctWnone(ixs(1):nextStable))); % within this set, the window with least recurrence
                        ixPerturb = ixs(1) + temp(1)-1;  % get window index
                        tsPerturb = [tsPerturb; Data(iR).stimbins(RcrPtData(iR).winmids(ixPerturb))];
                        ixs(ixs < nextStable) = []; % remove all indices from this section
                    else
                        % then doesn't return to attractor state: so
                        % store...
                        % time it falls below density threshold (i.e. after last "stable" window)...
                        % tsOff = Data(iR).stimbins(RcrPtData(iR).winmids(allstable(end)+1));
                        % time of final minimum density before end of
                        % recording...
                        temp = find(PrctWnone(ixs(1):end) == max(PrctWnone(ixs(1):end))); % within this set, the window with least recurrence
                        ixPerturb = ixs(1) + temp(1)-1;  % get window index

                        tsOff = [Data(iR).stimbins(RcrPtData(iR).winmids(ixs(1))) Data(iR).stimbins(RcrPtData(iR).winmids(ixPerturb))];  % last time it falls below perturbation threshold...
                        % keyboard
                    end

                    if isempty(ixs) || isempty(nextStable)
                        % if checked all windows of low recurrence, or doesn't
                        % return to attractor, then stop
                        blnCont = 0;
                    end
                end
            end
            if isempty(tsPerturb) tsPerturb = nan; end
            
       end
       RcrStats(iR,iF).firstW = tsFirst;
       RcrStats(iR,iF).tsPerturb = tsPerturb;
       RcrStats(iR,iF).tsOff = tsOff;
%    end
% end

       
   %% do time-series plots over windows including all periods
       cmap = flipud(cbrewer('seq','Reds',RcrPtData(iR).Orbit(statpars.ixT).nPeriods)); % most dense first
       winL = RcrPtData(iR).winDT * Data(iR).GaussQt;
       figure; 
       ylim = [0 0]; 
       ylimRot = 0; ylimT = 0; ylimW = 0;
       
       % keyboard

       % plot proportion of recurrence
       subplot(411), hold on
       PrctWall = 100*(1-RcrPtData(iR).Orbit(statpars.ixT).NoRecurrence.nWinPoints/RcrPtData(iR).winDT);
       PrctWper = 100*RcrPtData(iR).Orbit(statpars.ixT).Period(ixPrd).nWinPoints/RcrPtData(iR).winDT;
       ylimW = 100; %max(PrctW);
       plot(Data(iR).stimbins(RcrPtData(iR).winmids),PrctWall,'o','Color',AllColor)
       plot(Data(iR).stimbins(RcrPtData(iR).winmids),PrctWper,'o','Color',PerColor)
%        PrctWnoper = 100*RcrPtData(iR).Orbit(statpars.ixT).NoPeriod.nWinPoints/RcrPtData(iR).winDT;
%        plot(Data(iR).stimbins(RcrPtData(iR).winmids),PrctWnoper,'o','Color',NoPerColor)
       
       ylabel('#recurrence (% window)')     
       title(['Recording ' num2str(iR) ' ' fnames{iF}])
       
       % keyboard
       
       RcrStats(iR,iF).tsAllRots = zeros(numel(RcrPtData(iR).Tdelay(:,statpars.ixT)),1);
       for iP = 1:RcrPtData(iR).Orbit(statpars.ixT).nPeriods  
           RcrStats(iR,iF).tsAllRots(RcrPtData(iR).Orbit(statpars.ixT).Period(iP).IDs) = RcrPtData(iR).Orbit(statpars.ixT).Period(iP).RotPeriod;
           
           % plot
           subplot(412), hold on
           plot(Data(iR).stimbins(RcrPtData(iR).winmids),RcrPtData(iR).Orbit(statpars.ixT).Period(iP).mRotPeriod,'o','Color',cmap(iP,:))
           ylimRot = max(ylimRot,max(RcrPtData(iR).Orbit(statpars.ixT).Period(iP).mRotPeriod));
           ylabel('Rotation period (s)')
           subplot(413), hold on
           PrctT = 100*RcrPtData(iR).Orbit(statpars.ixT).Period(iP).Tstable/winL;
           plot(Data(iR).stimbins(RcrPtData(iR).winmids),PrctT,'o','Color',cmap(iP,:))
           ylimT = max(ylimT,max(PrctT));
           ylabel('Time stable (% window)')
           
%            subplot(413), hold on
%            PrctW = 100*RcrPtData(iR).Orbit(iT).Period(iP).nWinPoints/RcrPtData(iR).winDT;
%            plot(Data(iR).stimbins(RcrPtData(iR).winmids),PrctW,'o','Color',cmap(iP,:))
%            ylimW = max(ylimW,max(PrctW));
%            ylabel('#recurrence (% window)')

           subplot(414), hold on
           plot(Data(iR).stimbins(RcrPtData(iR).winmids),RcrPtData(iR).Orbit(statpars.ixT).Period(iP).mEigR,'o','Color',cmap(iP,:))
           ylabel('Mean Real Eig')
           xlabel('Time (s)')
           ylim = [min(ylim(1),min(RcrPtData(iR).Orbit(statpars.ixT).Period(iP).mEigR)) max(ylim(2),max(RcrPtData(iR).Orbit(statpars.ixT).Period(iP).mEigR))];
           line([30 Data(iR).stimbins(RcrPtData(iR).ixEnd)],[0 0],'Color',[0.5 0.5 0.5])
 
%            figure; hold on
%            title(['Recording ' num2str(iR) ' ' fnames{iF}])
%            xlabel('Time (s)'); ylabel('Rotation period (s)')

%            % get first window with density points > X 
%            tsFirst = Data(iR).stimbins(RcrPtData(iR).winmids(find(PrctW > pars.pDensity,1,'first')));
%            RcrStats(iR,iF).allfirstW(iP) = ;
%            if iP == ixPrd
%                 RcrStats(iR,iF).firstW = RcrStats(iR,iF).allfirstW(ixPrd);
%            end
       end
       
       subplot(411);
       set(gca,'XLim',[30 max(Data(iR).stimbins)])
       if ylimW > 0 
           set(gca,'YLim',[0 ylimW]); 
           line([Data(iR).stimbins(RcrPtData(iR).ix0) Data(iR).stimbins(RcrPtData(iR).ix0)],[0 ylimW],'Color',[0.5 0.5 0.5]);
           line([Data(iR).stimbins(RcrPtData(iR).ixEnd) Data(iR).stimbins(RcrPtData(iR).ixEnd)],[0 ylimW],'Color',[0.5 0.5 0.5]);
           % add line for coalescence 
           if ~isnan(RcrStats(iR,iF).firstW)
                line([RcrStats(iR,iF).firstW RcrStats(iR,iF).firstW],[0 ylimW],'Color',[0.3 0.5 1]);
           end       
           
           % add line for IDed perturbation time
           if ~isnan(RcrStats(iR,iF).tsPerturb)
               for i = 1:numel(RcrStats(iR,iF).tsPerturb)
                    line([RcrStats(iR,iF).tsPerturb(i) RcrStats(iR,iF).tsPerturb(i)],[0 ylimW],'Color',[0.3 0.8 0.2]);
               end
           end
       end      

       % now add no-period and no-recurrence points
       RcrStats(iR,iF).tsAllRots(RcrPtData(iR).Orbit(statpars.ixT).NoPeriod.IDs) = RcrPtData(iR).Orbit(statpars.ixT).NoPeriod.RotPeriod;

       subplot(412), hold on
       plot(Data(iR).stimbins(RcrPtData(iR).winmids),RcrPtData(iR).Orbit(statpars.ixT).NoPeriod.mRotPeriod,'o','Color',NoPerColor)
       
       subplot(413), hold on
       PrctT = 100*RcrPtData(iR).Orbit(statpars.ixT).NoPeriod.Tstable/winL;
       plot(Data(iR).stimbins(RcrPtData(iR).winmids),PrctT,'o','Color',NoPerColor)
       
       subplot(414), hold on
       plot(Data(iR).stimbins(RcrPtData(iR).winmids),RcrPtData(iR).Orbit(statpars.ixT).NoPeriod.mEigR,'o','Color',NoPerColor)
       ylim = [min(ylim(1),min(RcrPtData(iR).Orbit(statpars.ixT).NoPeriod.mEigR)) max(ylim(2),max(RcrPtData(iR).Orbit(statpars.ixT).NoPeriod.mEigR))];

       
       % add limits and lines over all periods
       subplot(412); % lines for start and end of recurrence counting
       set(gca,'XLim',[30 max(Data(iR).stimbins)])
       set(gca,'YLim',[0 ylimRot])
       line([Data(iR).stimbins(RcrPtData(iR).ix0) Data(iR).stimbins(RcrPtData(iR).ix0)],[0 ylimRot],'Color',[0.5 0.5 0.5]);
       line([Data(iR).stimbins(RcrPtData(iR).ixEnd) Data(iR).stimbins(RcrPtData(iR).ixEnd)],[0 ylimRot],'Color',[0.5 0.5 0.5]);
       
       subplot(413);
       set(gca,'XLim',[30 max(Data(iR).stimbins)])
       if ylimT > 0 
           set(gca,'YLim',[0 ylimT]); 
           line([Data(iR).stimbins(RcrPtData(iR).ix0) Data(iR).stimbins(RcrPtData(iR).ix0)],[0 ylimT],'Color',[0.5 0.5 0.5]);
           line([Data(iR).stimbins(RcrPtData(iR).ixEnd) Data(iR).stimbins(RcrPtData(iR).ixEnd)],[0 ylimT],'Color',[0.5 0.5 0.5]);
       end
       
       
       subplot(414); 
       set(gca,'XLim',[30 max(Data(iR).stimbins)])
       % set(gca,'YLim',ylim)
       %line([Data(iR).stimbins(RcrPtData(iR).ix0) Data(iR).stimbins(RcrPtData(iR).ix0)],ylim,'Color',[0.5 0.5 0.5]);
       %line([Data(iR).stimbins(RcrPtData(iR).ixEnd) Data(iR).stimbins(RcrPtData(iR).ixEnd)],ylim,'Color',[0.5 0.5 0.5]);
    
        
%         % raw data
%         figure; hold on
%         plot(Data(iR).stimbins(RcrPtData(iR).Orbit(iT).Period(1).IDs+RcrPtData(iR).ix0),real(RcrPtData(iR).Orbit(iT).Period(1).Emax),'k.') 
%         plot(Data(iR).stimbins(RcrPtData(iR).Orbit(iT).Period(2).IDs+RcrPtData(iR).ix0),real(RcrPtData(iR).Orbit(iT).Period(2).Emax),'g.') 
%         plot(Data(iR).stimbins(RcrPtData(iR).Orbit(iT).NoPeriod.IDs+RcrPtData(iR).ix0),real(RcrPtData(iR).Orbit(iT).NoPeriod.Emax),'r.')        
%         if iR == 10 && iF == 3
%             subplot(411),title('');
%             exportPPTfig(gcf,['Recur_time_series_Prep' num2str(iR) '_da' num2str(iF)],[15 5 9 13])
% 
%         end

%           exportPPTfig(gcf,['Recur_time_series_Prep' num2str(iR) '_da' num2str(iF)],[15 5 9 13])
        
          
        %% plot the relationship between rotation period and recurrence 
        patchinfo = struct('Vertices',[],'Faces',[],'FaceColor',[],'EdgeColor',[]);
        barmap = cbrewer('seq','PuBu',RcrPtData(iR).Orbit(statpars.ixT).nPeriods);
        
        figure
        subplot(211),bar(pars.bins,100*RcrPtData(iR).hstRecurrence(statpars.ixT,:),'FaceColor',[1 1 1]);  
        patchinfo.Vertices = [0 0; ... 
                           pars.tmin 0; ...
                           pars.tmin 25; ...
                           0 25];
        patchinfo.Faces = [1 2 3 4];
        patchinfo.FaceColor = [0.8 0.8 0.8];
        patchinfo.EdgeColor = 'none';
        patch(patchinfo)
        
%         ixFilled = find(RcrPtData(iR).hstRecurrence(statpars.ixT,:) > 0);
%         ixGap = diff(ixFilled);
        for iP = 1:RcrPtData(iR).Orbit(statpars.ixT).nPeriods 
            x = RcrPtData(iR).Orbit(statpars.ixT).Period(iP).mTime;
            line([x x],[0 25],'Color',[1 0 0])
            
            % add bar plot separate bar groups per period - tricky...
        end
        xlabel('Recurrence time (s)'); ylabel('Proportion of recurrence points (%)')
        
        subplot(212),plot(RcrPtData(iR).Tdelay(:,statpars.ixT),RcrStats(iR,iF).tsAllRots,'k.')
        xlabel('Recurrence time (s)'); ylabel('Rotation period (s)')
        title('Comparing recurrence time and linear systems estimates')
   end
   % pause
end

save([fname_Pre 'RecurrenceStats_FilterPts'],'RcrStats','statpars');

close all

%% plot summary stats
P10preps = [2,4,5]; % preps with P10

cmap = flipud(cbrewer('seq','Blues',nfiles)); 
cmap = cbrewer('qual','Paired',nfiles); 

% summary plot: 
% (1) color code by prep to show consistency within prep
% (2) link within prep to draw eye to groups
% (3) markersize according to stim sequence: no patterns
% (4) markersize according to "sensitised"

% Msize = [5 7 9];  % stim order
Msize = [3 9 5];   % for 5HT & washout: pick out the 5HT prep
% Msize = [5 15 10];  % does sensitised prep show clear effect?

hm = figure; hold on
hs = figure; hold on
hstab = figure; hold on
hF = figure; hold on

sym = 'o';
fLabel = 2; % which stim to label?

for iR = 1:nfiles
    % ID P10 preps using squares
%     if any(iR == P10preps')
%         sym = 's'; 
%     else
%         sym = 'o';
%     end
    
    figure(hm)
    plot([RcrStats(iR,:).mRot],[RcrStats(iR,:).mEigR],'o-',...
        'Color',cmap(iR,:),'MarkerSize',2,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
        
    figure(hs)
    CVrot = abs([RcrStats(iR,:).sdRot] ./ [RcrStats(iR,:).mRot]);
    CVeigR = abs([RcrStats(iR,:).sdEigR] ./ [RcrStats(iR,:).mEigR]);
%     plot([RcrStats(iR,:).sdRot],[RcrStats(iR,:).sdEigR],'o-',...
%       'Color',cmap(iR,:),'MarkerSize',2,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
    plot(CVrot,CVeigR,'o-',...
        'Color',cmap(iR,:),'MarkerSize',2,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);

    figure(hstab)
    % plot(1-[RcrStats(iR,:).densNoRecur],[RcrStats(iR,:).prctStable],'o-',...      
    plot([RcrStats(iR,:).densAllPeriod],[RcrStats(iR,:).prctStable],'o-',... 
        'Color',cmap(iR,:),'MarkerSize',2,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
    
   figure(hF)
   % plot(1-[RcrStats(iR,:).densNoRecur],[RcrStats(iR,:).firstW],'o-',...
   plot([RcrStats(iR,:).densAllPeriod],[RcrStats(iR,:).firstW],'o-',...
        'Color',cmap(iR,:),'MarkerSize',2,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
   
   
    % increasing MarkerSize for sequence
    for iF = 1:numel(fnames)
        figure(hm)
        plot(RcrStats(iR,iF).mRot,RcrStats(iR,iF).mEigR,sym,...
            'Color',cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',cmap(iR,:));
        if iF == fLabel
            text([RcrStats(iR,iF).mRot],[RcrStats(iR,iF).mEigR],num2str(iR),'Fontsize',8,'Color',[1 1 1])
        end
        
        
        figure(hs)
%         plot([RcrStats(iR,iF).sdRot],[RcrStats(iR,iF).sdEigR],sym,...
%             'Color',cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',cmap(iR,:));
        plot(CVrot(iF),CVeigR(iF),sym,...
            'Color',cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',cmap(iR,:));
        if iF == fLabel
            text(CVrot(iF),CVeigR(iF),num2str(iR),'Fontsize',8,'Color',[1 1 1])
        end

        figure(hstab)
        plot(RcrStats(iR,iF).densAllPeriod,RcrStats(iR,iF).prctStable,sym,...
            'Color',cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',cmap(iR,:));
        if iF == fLabel
            text(RcrStats(iR,iF).densAllPeriod,RcrStats(iR,iF).prctStable,num2str(iR),'Fontsize',8,'Color',[1 1 1])
        end
        
        figure(hF)
        plot(RcrStats(iR,iF).densAllPeriod,RcrStats(iR,iF).firstW,sym,...
            'Color',cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',cmap(iR,:));
        if iF == fLabel
            text(RcrStats(iR,iF).densAllPeriod,RcrStats(iR,iF).firstW,num2str(iR),'Fontsize',8,'Color',[1 1 1])
        end

    end
end
figure(hm)
axis tight
xlim = get(gca,'XLim');
set(gca,'XLim',[0 xlim(2)]);
line([0 xlim(2)],[0 0],'Color',[0.5 0.5 0.5])
xlabel('Orbit rotation period (s)')
ylabel('Rate of contraction (real eigenvalue)')
exportPPTfig(gcf,[fname_Pre 'rotation_vs_contraction'],[10 15 7 7])

figure(hs)
xlabel('CV of Orbit rotation period')
ylabel('CV of Rate of contraction (real eigenvalue)')

figure(hstab)
xlabel('Density of recurrence points')
ylabel('Proportion of time stable')
%exportPPTfig(gcf,'density_vs_time_stable',[10 15 7 7])

figure(hF)
xlabel('Density of recurrence points')
ylabel('Time to coalscence (s)')
%title(['Coalscence as first window with ' num2str(pars.pDensity) '% recurrence'])
exportPPTfig(gcf,[fname_Pre 'density_vs_coalescence'],[10 15 7 7])

