%% check out jumps: perturbations or transitions
%% 2nd edition: checking properly 

clear all; close all

fname_Pre = '';
fnames = {'da01','da02','da03'};

% load recurrence stats
% load RecurrenceStats 
load([fname_Pre 'RecurrenceStats_FilterPts'])  % use recurrence point statistics based on using more robust model fits - min 1s and contiguous points

npreps = size(RcrStats,1);

% parameters
jumppars.theta = 0.2;  % threshold for a different attractor
jumppars.cycles = 2;  % how many cycles before end...

% for each prep, for each program
ctr = 1;
Summary.allJumps = 0;
Summary.JumpList = [];  % list of all jumps (Prep,Program,Jump Number)
Summary.blnReturned = []; % returned to attractor?
Summary.blnSame = [];   % same or different attractor?
Summary.pAfter = [];
Summary.pAfterEnd = [];

PropPeriodLocal = zeros(npreps,numel(fnames));
PropPeriodHist = zeros(npreps,numel(fnames));

for iPrep = 1:npreps
    iPrep
   for iProg = 1:numel(fnames)
       % if exists, load recurrence data
       load([fname_Pre fnames{iProg} '_RecurrencePointsData_FilterPts'],'RcrPtData','pars');
       % load([fnames{iF} '_RecurrencePointsData_Kmeans']);
       load(['../' fname_Pre fnames{iProg} '_DataProperties_FunctionAndWindowSize'],'Data','FileTable')           

        % time-series of recurrence density per window
        PrctWrecur = 100* RcrPtData(iPrep).Orbit(statpars.ixT).AllRecurrentPoint.pRecur;
        
        nTs = RcrPtData(iPrep).ixEnd - RcrPtData(iPrep).ix0;  % time-steps of checked recurrence points
        allIDs = 1:nTs;  % indexes of checked recurrence points
        
        % for all jumps returning to attractor: check if they reach the
        % same or a different attractor
       if ~isnan(RcrStats(iPrep,iProg).tsPerturb)

           for iT = 1:numel(RcrStats(iPrep,iProg).tsPerturb)
               Summary.allJumps = Summary.allJumps+1;
               Summary.blnReturned = [Summary.blnReturned; 1]; 
               Summary.JumpList = [Summary.JumpList; iPrep iProg iT];

               % get index of time(s) in winmids...
               ixP = find(Data(iPrep).stimbins == RcrStats(iPrep,iProg).tsPerturb(iT)); % index of jump time
               ixW = find(RcrPtData(iPrep).winmids == ixP);   % index of window of jump time
               Jump(iPrep,iProg).Transition(iT).tsPerturbMin = RcrStats(iPrep,iProg).tsPerturb(iT);
               
               % find Jump start, and window just prior to it                             
               allposs = find(PrctWrecur(1:ixW-1) > statpars.pDensity); % all possible windows before the threshold crossing (ixW)
               Jump(iPrep,iProg).Transition(iT).ixWinJumpStart = allposs(end)+1;
               Jump(iPrep,iProg).Transition(iT).tsJumpStart = Data(iPrep).stimbins(RcrPtData(iPrep).winstrt(allposs(end)+1)); % start of jump period
               Jump(iPrep,iProg).Transition(iT).ixWinBefore = allposs(end); % point just before jump
               Jump(iPrep,iProg).Transition(iT).tsWinBefore = Data(iPrep).stimbins(RcrPtData(iPrep).winstrt(allposs(end))); % start of jump period
              
               % find Jump end
               allposs = find(PrctWrecur(ixW+1:end) > statpars.pDensity) + ixW;
               Jump(iPrep,iProg).Transition(iT).ixWinJumpEnd = allposs(1);
               Jump(iPrep,iProg).Transition(iT).tsJumpEnd = Data(iPrep).stimbins(RcrPtData(iPrep).winstrt(allposs(1))); % end of jump period
               
               
               % Same or different attractor?: find what proportion of delays in that window are after end of jump
               % period, or within it
               startwin = RcrPtData(iPrep).winstrt(Jump(iPrep,iProg).Transition(iT).ixWinBefore);
               endwin = startwin + RcrPtData(iPrep).winDT;
               ixPts = find(allIDs > startwin-RcrPtData(iPrep).ix0 & allIDs <= endwin-RcrPtData(iPrep).ix0);
               Delays = RcrPtData(iPrep).Tdelay(ixPts,statpars.ixT); % set of all delays in this window
               % time of every trajectory point (absilute time since start of recording
               tsThesePoints = (ixPts + RcrPtData(iPrep).ix0) * Data(iPrep).GaussQt + pars.stimstart;
               % next recurrence time of each point on the trajectory
               NextRecur = tsThesePoints' + Delays;
               NextRecur = NextRecur(~isnan(NextRecur));  % remove non-recurrent....
               Summary.pAfter = [Summary.pAfter; sum(NextRecur > RcrStats(iPrep,iProg).tsPerturb(iT)) ./ numel(NextRecur)];
               Summary.pAfterEnd = [Summary.pAfterEnd; sum(NextRecur > Jump(iPrep,iProg).Transition(iT).tsJumpEnd) ./ numel(NextRecur)];
                            
               % after = return to same attractor
               % within = different attractor
               Summary.blnSame = [Summary.blnSame; Summary.pAfter(end) > jumppars.theta];
              
%                allSummary(ctr,:) = [Jump(iPrep,iProg).Transition(iT).medDmin, Jump(iPrep,iProg).Transition(iT).iqrDmin, iPrep, iProg, iT, ...
%                    Jump(iPrep,iProg).Transition(iT).nBefore Jump(iPrep,iProg).Transition(iT).nAfter Jump(iPrep,iProg).Transition(iT).densityDmin, ...
%                    Jump(iPrep,iProg).Transition(iT).pBeforeRot, Jump(iPrep,iProg).Transition(iT).pAfterRot]; 
%                ctr = ctr+1;
           end
       end
       
               
        % is there a final "off" trajectory perturbation? i.e. permanent
        % divergence from the attractor?
        if ~isnan(RcrStats(iPrep,iProg).tsOff) 
            % orbit period as mean of all main Period recurrence points
            nPoints = [RcrPtData(iPrep).Orbit(statpars.ixT).Period(:).nPoints];
            ixPrd = nPoints == max(nPoints);  % dominant period
            mPeriod(iPrep,iProg) = RcrPtData(iPrep).Orbit(statpars.ixT).Period(ixPrd).mTime; % mean of periods in main peak
            % local orbit period as mean recurrence time in the final stable window...
            allposs = find(PrctWrecur > statpars.pDensity,1,'last');
            localcycle(iPrep,iProg) = RcrPtData(iPrep).Orbit(statpars.ixT).AllRecurrentPoint.mDelay(allposs(end));  % mean delay in the final stable window window
            
            
            % there is a final perturb threshold crossing that is not followed by recurrence
            PropPeriodHist(iPrep,iProg,1) = (Data(iPrep).bins(end)-RcrStats(iPrep,iProg).tsOff(1)) ./ mPeriod(iPrep,iProg);
            % minimum recurrence window
            PropPeriodHist(iPrep,iProg,2) = (Data(iPrep).bins(end)-RcrStats(iPrep,iProg).tsOff(2)) ./ mPeriod(iPrep,iProg);

            PropPeriodLocal(iPrep,iProg,1) = (Data(iPrep).bins(end)-RcrStats(iPrep,iProg).tsOff(1)) ./ localcycle(iPrep,iProg);
            PropPeriodLocal(iPrep,iProg,2) = (Data(iPrep).bins(end)-RcrStats(iPrep,iProg).tsOff(2)) ./ localcycle(iPrep,iProg);

        end
        
        % if so, does it fall clearly before the finite time effect of the end of the recording? 
        % if PropPeriodLocal(iPrep,iProg,1) > jumppars.cycles % criterion: use perturb crossing or final min point?
        if PropPeriodHist(iPrep,iProg,1) > jumppars.cycles % criterion: use perturb crossing or final min point?

            Summary.JumpList = [Summary.JumpList; iPrep iProg 0];
            Summary.blnReturned = [Summary.blnReturned; 0]; 
            Summary.allJumps = Summary.allJumps + 1;  % if so: add to all jumps
            Summary.pAfter = [Summary.pAfter; nan];
            Summary.pAfterEnd = [Summary.pAfterEnd; nan];
            Summary.blnSame = [Summary.blnSame; nan];
        end
        

   end
end

%% summary plot


save([fname_Pre 'JumpAnalysis2'],'Jump','Summary','jumppars')

