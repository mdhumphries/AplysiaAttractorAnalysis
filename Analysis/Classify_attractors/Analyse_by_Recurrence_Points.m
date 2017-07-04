% script to classify projections
% does periodic orbit recurrence ~ Lathrop & Kostelich (1989)
% ADD osc and burst to this: 
% 
% NB: for a given point in PC projection, computes its recurrence time for
% a given threshold; first discards all time-points less then threshold
% that are contiguous in future time with the given point - i.e. removes
% the local tangent
% 
% pars.X
%   struct of analysis parameters:
%         stimstart: time of stimulation onset (s)
%                T0: initial time of recurrence checking (s)
%         maxWindow: time before end of time-series to stop checking for recurrence (s)
%      VarExplained: target proportion of variance explained by chosen set of PC dimensions
%         prctTheta: vector (Nx1): percentages of recurrence distance
%                           distribution to use as recurrence threshold
%             bstep: step-size for recurrence point histogram (s)
%              bins: vector of bins for recurrence point histogram (s)
%              tmin: smallest recurrence time to count as a histogram peak (s)
%            minbin: index of bin which corresponds to tmin
%              Nmin: min. number of points in a peak to include it in analysis
%              Kmax: maximum number of k-means clusters to use in consensus
%                       clustering of recurrence points
%        NeighbourD: mulitiplier of current recurrence point's time to get
%                           time-series neighbourhood
%     minNeighbours: min. number of found neighbours to fit linear model
%              winL: length of sliding window (s)
%           winstep: step of sliding window (s)
%
% RcrPtData 
%   for each recording (iR), stores:
%              maxD: maximum dimension of state-space (number of PCs),
%                       according to pars.VarExplained
%           currPCs: set of PCs (1...D) since stimulation (normalised, from PRINCOMP)
%             Theta: set of recurrence thresholds (defined by set of
%                       percentiles in pars.prctTheta
%               ix0: index in currPCs of pars.T0  [NB wrt to pars.stimstart]
%             ixEnd: index in currPCs of last checked point for recurrence
%            Tdelay: NxM matrix of delays between current point and
%                       recurrence point for given threshold (N thresholds; M checked points).
%     hstRecurrence: NxB matrix of recurrence point histograms (N thresholds; B bins)
%             Orbit: struct of recurrence point analysis - see below
%             winDT: number of time-steps per window
%           winstrt: start times of each sliding window (s)
%           winmids: mid-point times of each sliding window (s)
%
% Orbit
%   for each recording (iR), for each threshold (iH) (theta), stores:
%          nPeriods: the number of different orbit periods (peaks in histogram)
%            Period: struct of dynamical results for all recurrence points in an orbit period (see below)
%          NoPeriod: struct of dynamical results for the remaining recurrence points not in any orbit period
%      NoRecurrence: struct of dynamical results for all non-recurrent points
%           frstRcr: over all recurrence points in orbit periods, the
%                       earliest recurrence point
%
% AllRecurrentPoint
%       for each recording (iR), for each threshold (iH) (theta), stores:
%         pRecur: proportion of trajectory points in this window that recur
%         mDelay: mean recurrence time in each sliding window
%        sdDelay: SD of recurrence time in each sliding window
%
% Period
%   for each recording (iR), for each threshold (iH) (theta), for each orbit period (iP), stores:
%        nPoints: number of points in this orbit period
%          mTime: mean recurrence time of points in this period
%        frstRcr: earliest recurrence point of this period
%            IDs: vector of all recurrence points in this orbit period (IDs into Tdelay matrix)
%           Npts: number of neighbourhood points used to esimate local linear model (per recurrence point)
%           Emax: maximum eigenvalue of linearised dynamics in
%                       neighbourhood of recurrent point (per recurrence point)
%         allEgs: matrix of all eigenvalues of linearised dynamics (M dims
%                       x N points) (per recurrence point)
%           Cond: condition number of A matrix (per recurrence point)
%             R2: variance-explained by fit of model to all dimensions (per recurrence point)  
%      RotPeriod: estimated rotation period, from imaginary part of max eigenvalue
%     nWinPoints: number of recurrence points from this orbit in each sliding window
% nWinUsedPoints: number of recurrence points with enough neighbourhood points to estimate local linear model
%     mRotPeriod: mean rotation period per window (from used points)
%        Tstable: duration of window with no rotation period (from used points)
%          mEigR: mean of real part of eigenvalue per window (from used points)
%          mEigI: mean of imaginary part of eigenvalue per window (from used points)
%       mEigRrot: mean of real part of eigenvalue per window from all points with rotation (from used points)
%       mEigIrot: mean of imaginary part of eigenvalue per window from all points with rotation (from used points)
%        mEigAbs: mean of abs(eigenvalue) (i.e. modulus) per window (from used points)
%
% NoPeriod, No Recurrence
%     as for struct Period, without "mTime" and "firstRcr" fields; and
%     without mEigRot, mEigIrot
%
% NOTES:
%   (1) points need not be recurrent! So all points in a given
%   window need not be in the combined Periods and the NoPeriod sets
%   (2) Updated 16/10/15 to: (i) only use contiguous periods for fitting
%   local linear model and (ii) using min. number of points to reject local
%   fits of less than 1 second of data. Resulting files now saved with
%   "Filter_Pts" on data-file
%   (3) Updated 5/11/15 to return every eigenvalue, so we can track all of
%   these in time if need be (and compute dominance)
%   (4) Updated 16/11/15 to return condition number of every A matrix (from
%       x(t+1) = Ax(t); and R^2 of fit of linear model around recurrence
%       point [R^2 computed from forward and backward solution, and correlated in all dimensions: to be double-checked]
%
%   (5) Updated 5/7/2016 to compute sliding windows of recurrence delay for
%   all delay points, and for each Period
% Mark Humphries 5/7/16

clear all; close all

addpath ../../Functions/

fname = 'da03';
prefix = ''; % '5HT';  % ''; for main data

% get PCA projections - main data
load(['../' prefix fname '_DataProperties_FunctionAndWindowSize'],'PCAdata','Data','SDF')


nfiles = numel(PCAdata);

% % get cell types
% load(['../Ensembles/' prefix fname '_Analyses_Neurons_and_Groups'],'groupdata')
% load(['../Ensembles/' prefix fname '_Ensemble_Types'],'AcorrTypes')

% parameters
pars.stimstart = 30; % (s)
pars.T0 = 35;  % start counting after stim (s)
pars.maxWindow = 10;  % stop before end of time-series (s)
pars.VarExplained = 0.80;    % how many PCs to keep

pars.prctTheta = [2 5 10]; %[1,3,5]; %  % of distance distribution to use as threshold; max distance to check (spikes/s)
% maxTheta = 10; % max distance to check (spikes/s)
% stepTheta = 0.5;
% Theta = stepTheta:stepTheta:maxTheta;

pars.bstep = 1;  % s
pars.bins = 0:pars.bstep:50;  % recurrence time pars.bins (s);

pars.tmin = 5;  % s; minimum recurrence time to count as a peak
pars.Nmin = 100;      % minimum number of points to use that peak
pars.Kmax = 5;   % maximum number of k-means to try

pars.NeighbourD = 2.5;  % multiplier of threshold to use as neighbourhood around (m,e) point
pars.minNeighbours = 100; % minimum number of points to regression = 1 second of time around (m,e) point

pars.winL = 5;  %sliding window in which to compute local rotation time and average eigenvalues
pars.winstep = 1;  % step of window

pars.minbin = pars.bins(find(pars.bins <= pars.tmin,1,'last')); % lowest retained pars.bins in histogram

%% find attractor points
for iR = 1:nfiles
    iR
%     %% do oscillator PCA and projection
%     % ensemble IDs of this recording in groupdata structure
%     ixGrps = find([groupdata(:).Recording] == iR);
% 
%     % ensemble-types of this subset of ensembles: 
%     osctypes = zeros(numel(ixGrps),1); 
%     for iG = 1:numel(ixGrps)
%         osctypes(iG) = AcorrTypes(ixGrps(iG));
%     end
%     
%     
%     % select subset of groups for PCA: specify here what to use, then save
%     % with appropriate names.... NOT CURRENTLY USED
%     %PCAsubset(iR).ixGrps = ixGrps(find(osctypes == 1)); % tonic...
%     PCAsubset(iR).ixGrps = ixGrps(find(osctypes == 2)); % use only groups with significant peak AND trough
%     % PCAsubset(iR).ixGrps = ixGrps(find(osctypes == 3)); % bursters
%     
%     % now look up those neuron IDs
%     PCAsubset(iR).neuronIDs = [];
%     for iP = 1:numel(PCAsubset(iR).ixGrps)
%         PCAsubset(iR).neuronIDs = [PCAsubset(iR).neuronIDs; groupdata(PCAsubset(iR).ixGrps(iP)).IDs];
%     end
    [r,N] = size(SDF(iR).spkfcn);
    Tend = min(r,numel(Data(iR).bins)); % urrgh - fix this in DataProperties...
    ixSpon = find(Data(iR).bins < pars.stimstart);
    ixStim = find(Data(iR).bins(1:Tend) >= pars.stimstart);
    
%     % get PCA of stim period
%     [PCAsubset(iR).coeffs, PCAsubset(iR).scores, PCAsubset(iR).eigvalues] = princomp(SDF(iR).spkfcn(ixStim,PCAsubset(iR).neuronIDs));
%      PCAsubset(iR).prctVar = cumsum(PCAsubset(iR).eigvalues) / sum(PCAsubset(iR).eigvalues);
%      PCAsubset(iR).maxD = sum(PCAsubset(iR).prctVar <= pars.VarExplained);
%     

    % pick the subset of dimensions for whole thing 
    RcrPtData(iR).maxD = sum(PCAdata(iR).prctVar <= pars.VarExplained);
    RcrPtData(iR).currPCs = PCAdata(iR).PCproj(:,1:RcrPtData(iR).maxD);  % from stimulation
    
    
    % get thresholds
    dPs = pdist(RcrPtData(iR).currPCs,'euclidean');
    RcrPtData(iR).Theta = prctile(dPs,pars.prctTheta);
    
    % for every time-point, find nearest time-point within threshold
    RcrPtData(iR).ix0 = find(Data(iR).stimbins == pars.T0); % initial time-point
    RcrPtData(iR).ixEnd = find(Data(iR).stimbins == Data(iR).stimbins(end)-pars.maxWindow);
    nTs = RcrPtData(iR).ixEnd - RcrPtData(iR).ix0;  % loop from start to before end
   
    for iN = 1:nTs
        pt0 = RcrPtData(iR).currPCs(RcrPtData(iR).ix0+iN,:);  % current x_i
        currixs = RcrPtData(iR).ix0+iN+1:numel(Data(iR).stimbins);
        D = sqrt(sum(bsxfun(@minus,RcrPtData(iR).currPCs(currixs,:),pt0).^2,2)); % Euclidean distance to all points ahead in time
        for iT = 1:numel(RcrPtData(iR).Theta)
            allpts = find(D <= RcrPtData(iR).Theta(iT));  % every point within this threshold
        
            % remove those immediately after T that are within theta -
            % short time-scale noise...
            
%             % Need the below for looking backward and forward in time
%             ix0pt = find(allpts == ix0);  % pars.T0 in this subset
%             dpts = diff(allpts);
%             ulim = find(dpts(ix0pt+1:end) > 1,1,'first'); if isempty(ulim) ulim = numel(allpts); end
%             llim = find(dpts(1:ix0pt-1) > 1,1,'last'); if isempty(llim) llim = 0; end
%             allpts(llim+1:ulim) = [];
            
            dpts = diff(allpts);
            ulim = find(dpts > 1,1,'first'); if isempty(ulim) ulim = numel(allpts); end
            allpts(1:ulim) = [];
            
            dly = min(Data(iR).stimbins(currixs(allpts))) - Data(iR).stimbins(RcrPtData(iR).ix0+iN); % keep shortest delay
            if dly                 
                % store delay: any point with delay is a (m,e) point by
                % definition
                RcrPtData(iR).Tdelay(iN,iT) = dly;
            else
                % nothing
                RcrPtData(iR).Tdelay(iN,iT) = nan;
            end
            
        end
        
    end
    
    
    
    %% get (m,theta) points: find period 1, 2 etc of recurrence points
    tic
    RcrPtData(iR).hstRecurrence = zeros(numel(RcrPtData(iR).Theta),numel(pars.bins));
    for iT = 1:numel(RcrPtData(iR).Theta)
         hst = histc(RcrPtData(iR).Tdelay(:,iT),pars.bins);
         RcrPtData(iR).hstRecurrence(iT,:) = hst / sum(hst); % approx PDF
         % find histogram peaks to get points in each period
         % (1) histogram edge finding
         % keyboard
  
        filledbins = [find(RcrPtData(iR).hstRecurrence(iT,pars.minbin:end) > 0)];
        % keyboard
        if ~isempty(filledbins)
            length = diff(filledbins);
            hstedges = [[filledbins(1); filledbins(find(length > 1)+1)'] [filledbins(length > 1)'; filledbins(end)]] + pars.minbin-1;
            npeaks = size(hstedges,1);
            ctr = 0;
            for i=1:npeaks
                % get all points within that histogram peak
                ixts = find(RcrPtData(iR).Tdelay(:,iT) >= pars.bins(hstedges(i,1)) & RcrPtData(iR).Tdelay(:,iT) <= pars.bins(hstedges(i,2))+1);

                if numel(ixts) >= pars.Nmin
                    % if enough points, keep as candidate Period N set
                    ts = RcrPtData(iR).Tdelay(ixts,iT);
                    ctr = ctr+1;
                    RcrPtData(iR).Orbit(iT).Period(ctr).IDs = ixts; % get IDs of these time-points in Tdelay array (need adjusting by ix0 to get time)
                    RcrPtData(iR).Orbit(iT).Period(ctr).mTime = median(ts);
                    RcrPtData(iR).Orbit(iT).Period(ctr).nPoints = numel(ts); 
                end

            end
            RcrPtData(iR).Orbit(iT).nPeriods = ctr;
            
        else
             RcrPtData(iR).Orbit(iT).nPeriods = 0;
        end
        
        % now assign all time-stamp IDs not in Periods to either:
        % (a) "NoPeriod": recurrent, but not strongly periodic
        % (b) "NoRecurrence": not a recurrent point.... (perturbation?)
        assignedIDs = []; allIDs = 1:nTs;
        for iP = 1:RcrPtData(iR).Orbit(iT).nPeriods
            assignedIDs = [assignedIDs; RcrPtData(iR).Orbit(iT).Period(iP).IDs];
        end
        RcrPtData(iR).Orbit(iT).NoRecurrence.IDs = find(isnan(RcrPtData(iR).Tdelay(:,iT)));
        RcrPtData(iR).Orbit(iT).NoRecurrence.nPoints = numel(RcrPtData(iR).Orbit(iT).NoRecurrence.IDs);
        RcrPtData(iR).Orbit(iT).NoPeriod.IDs = setdiff(allIDs,[assignedIDs; RcrPtData(iR).Orbit(iT).NoRecurrence.IDs]);
        RcrPtData(iR).Orbit(iT).NoPeriod.nPoints = numel(RcrPtData(iR).Orbit(iT).NoPeriod.IDs);
        
        
%          %% (2) kmeans - wroks very well, but MUCH slower, and more
%          % complex to explain
% 
%          X = RcrPtData(iR).Tdelay(~isnan(RcrPtData(iR).Tdelay(:,iT)),iT);  % all recurrence times 
%          X(X < pars.tmin) = []; % remove all that likely tangenital motion 
%          if ~isempty(X)
%              Idx = zeros(numel(X),pars.Kmax);
%              for k= 1:pars.Kmax
%                 [Idx(:,k),C,SumD]= kmeans(X,k,'Start','uniform');  % k-means at each k
%              end
%              Sc = consensusmatrix(Idx);  % get consensus matrix
%              Keep = Sc == 1;  % keep only 100% clustered
% 
%              grpscon = [0 0]; grpctr = 1; blnTrans = 1;
%             for iN = 1:numel(X)
%                 if ~any(iN == grpscon(:,1))  % if not already sorted
%                     thisgrp = [iN; find(Keep(iN,:) == 1)']; % all nodes connected to this one are in the same group
%                     % if any of this group are already in the storage matrix,
%                     % then this is not transitive...
%                     for iG = 1:numel(thisgrp)
%                         if any(thisgrp(iG) == grpscon(:,1))
%                             blnTrans = 0; break
%                         else
%                             blnTrans = 1;
%                         end
%                     end
%                     if blnTrans == 0
%                         break
%                     else
%                         grpscon = [grpscon; thisgrp ones(numel(thisgrp),1)*grpctr]; % save this...
%                         grpctr = grpctr + 1;
%                     end
%                 end
%             end
%             grpscon(1,:) = []; % remove blank entry
%             [~,iG] = sort(grpscon(:,1));  % sort into index order
%             grpscon = grpscon(iG,:);  % sort into index order
%             if ~blnTrans
%                 error('Not transitive???')
%             end
%             grpIDs = unique(grpscon(:,2));
%             nGrps = numel(grpIDs);
%             Gsizes = arrayfun(@(x) sum(grpscon(:,2) == x),grpIDs);
%             
%             % now get members and stats of each group
%             ctr = 0; 
%             for iG = 1:nGrps
%                 % check if group later than close trajectory
%                 thisG = find(grpscon(:,2) == grpIDs(iG));
%                 ts = X(thisG);
%                 % keep if all after min time AND enough entries
%                 if numel(ts) >= pars.Nmin
%                     ctr = ctr+1;
%                     ixNotNan = find(~isnan(RcrPtData(iR).Tdelay(:,iT)));
%                     RcrPtData(iR).Orbit(iT).Period(ctr).IDs = ixNotNan(thisG); % get IDs of these cells
%                     RcrPtData(iR).Orbit(iT).Period(ctr).mTime = median(ts);
%                     RcrPtData(iR).Orbit(iT).Period(ctr).nPoints = numel(ts); 
%                 end
% 
%             end
%             if ctr > 1
%                 % sort period into time order (Period 1, Period 2 etc..)
%                [~,iP] = sort([RcrPtData(iR).Orbit(iT).Period(:).mTime]);
%                RcrPtData(iR).Orbit(iT).Period = RcrPtData(iR).Orbit(iT).Period(iP);  
%             end
%             RcrPtData(iR).Orbit(iT).nPeriods = ctr;
%             
%          else
%              RcrPtData(iR).Orbit(iT).nPeriods = 0;
%          end

   end 
   toc
    
    % keyboard
    
    figure
    imagesc(pars.bins,RcrPtData(iR).Theta,RcrPtData(iR).hstRecurrence); 
    colormap(hot); colormap(gray)
    xlabel('Time delay (s)')
    ylabel('Threshold (spikes/s)')
        
   
   %% estimate stability     
   % get all recurrence points of Period 1 etc 
   
   RcrPtData(iR).winDT = pars.winL / Data(iR).GaussQt;  % time-steps per window
   stepDT = pars.winstep / Data(iR).GaussQt; % time-steps per window shift
   RcrPtData(iR).winstrt = RcrPtData(iR).ix0:stepDT:RcrPtData(iR).ixEnd-RcrPtData(iR).winDT;  % indices of window start times
   RcrPtData(iR).winmids = round(RcrPtData(iR).winstrt + RcrPtData(iR).winDT/2); 
   
   for iT = 1:numel(RcrPtData(iR).Theta) 
       
       % sliding window of recurrence-times
       for iW = 1:numel(RcrPtData(iR).winstrt)
            % allIDs runs from 1 to number of elements in Tdelay; 
            % Windows run from time-step of checking after stimulation
            % (i.e. pars.T0)
            ixPts = find(allIDs > RcrPtData(iR).winstrt(iW)-RcrPtData(iR).ix0 & allIDs <= RcrPtData(iR).winstrt(iW)+RcrPtData(iR).winDT-RcrPtData(iR).ix0);
            RcrPtData(iR).Orbit(iT).AllRecurrentPoint.mDelay(iW) = nanmean(RcrPtData(iR).Tdelay(ixPts,iT));
            RcrPtData(iR).Orbit(iT).AllRecurrentPoint.sdDelay(iW) = nanstd(RcrPtData(iR).Tdelay(ixPts,iT));
            RcrPtData(iR).Orbit(iT).AllRecurrentPoint.pRecur(iW) = sum(~isnan(RcrPtData(iR).Tdelay(ixPts,iT))) ./ RcrPtData(iR).winDT;  % should agree with NoRecurrence....
       end
       
       % for each period
       if RcrPtData(iR).Orbit(iT).nPeriods > 0
           % cmap = flipud(cbrewer('seq','Reds',RcrPtData(iR).Orbit(iT).nPeriods)); % most dense first
%            hf = figure; hold on
           for iP = 1:RcrPtData(iR).Orbit(iT).nPeriods  % for each recurrence period
               tic
               RcrPtData(iR).Orbit(iT).Period(iP).frstRcr = Data(iR).stimbins(min(RcrPtData(iR).Orbit(iT).Period(iP).IDs)+RcrPtData(iR).ix0); % first point that recurs
               [RcrPtData(iR).Orbit(iT).Period(iP).Emax,RcrPtData(iR).Orbit(iT).Period(iP).allEgs,RcrPtData(iR).Orbit(iT).Period(iP).RotPeriod,RcrPtData(iR).Orbit(iT).Period(iP).Npts,...
                   RcrPtData(iR).Orbit(iT).Period(iP).Cond,RcrPtData(iR).Orbit(iT).Period(iP).R2] = LocalLinearDynamics(RcrPtData(iR).Orbit(iT).Period(iP).IDs,...
                            RcrPtData(iR).currPCs,RcrPtData(iR).ix0,RcrPtData(iR).Theta(iT),pars.NeighbourD,Data(iR).GaussQt);
               toc
               
               
%                figure(hf)
%                plot(Data(iR).stimbins(RcrPtData(iR).Orbit(iT).Period(iP).IDs),RcrPtData(iR).Orbit(iT).Period(iP).RotPeriod,'o','Color',cmap(iP,:))
%                axis([Data(iR).stimbins(ix0),Data(iR).stimbins(ixEnd),0,max(RcrPtData(iR).Orbit(iT).Period(iP).RotPeriod)]) 
               
               % now compute time-series of eigenvalues and periods
               % keyboard
               for iW = 1:numel(RcrPtData(iR).winstrt)
                   % find all time-point IDs that ar within the current
                   % window: shift window edges to match time-base for the
                   % Tdelay array
                   ixPts = find(RcrPtData(iR).Orbit(iT).Period(iP).IDs > RcrPtData(iR).winstrt(iW)-RcrPtData(iR).ix0 & RcrPtData(iR).Orbit(iT).Period(iP).IDs <= RcrPtData(iR).winstrt(iW)+RcrPtData(iR).winDT-RcrPtData(iR).ix0);
                   RcrPtData(iR).Orbit(iT).Period(iP).nWinPoints(iW) = numel(ixPts);  % how many points in this window for this period
                   if isempty(ixPts)
                       RcrPtData(iR).Orbit(iT).Period(iP).mDelay(iW) = nan;
                       RcrPtData(iR).Orbit(iT).Period(iP).nWinUsedPoints(iW) = 0;
                       RcrPtData(iR).Orbit(iT).Period(iP).mRotPeriod(iW) = nan;
                       RcrPtData(iR).Orbit(iT).Period(iP).Tstable(iW) = nan;
                       RcrPtData(iR).Orbit(iT).Period(iP).mEigR(iW) = nan;
                       RcrPtData(iR).Orbit(iT).Period(iP).mEigI(iW) = nan;
                       RcrPtData(iR).Orbit(iT).Period(iP).mEigAbs(iW) = nan;
                   else
                       RcrPtData(iR).Orbit(iT).Period(iP).mDelay(iW) = nanmean(RcrPtData(iR).Tdelay(ixPts,iT));
                       theseN = RcrPtData(iR).Orbit(iT).Period(iP).Npts(ixPts);
                        usedPts = ixPts(theseN >= pars.minNeighbours);
                        RcrPtData(iR).Orbit(iT).Period(iP).nWinUsedPoints(iW) = numel(usedPts);
                        theseRotations = RcrPtData(iR).Orbit(iT).Period(iP).RotPeriod(usedPts);
                        RcrPtData(iR).Orbit(iT).Period(iP).mRotPeriod(iW) = mean(theseRotations(theseRotations > 0));  % local estimate of current rotational period
                        RcrPtData(iR).Orbit(iT).Period(iP).Tstable(iW) = sum(theseRotations == 0) * Data(iR).GaussQt;
                        theseEmax = RcrPtData(iR).Orbit(iT).Period(iP).Emax(usedPts);
                        RcrPtData(iR).Orbit(iT).Period(iP).mEigRrot(iW) = mean(real(theseEmax(theseRotations > 0)));  % mean of real part, of Eigs of rotation
                        RcrPtData(iR).Orbit(iT).Period(iP).mEigIrot(iW) = mean(imag(theseEmax(theseRotations > 0)));  % mean of imaginary part, of Eigs of rotation                                              
                        RcrPtData(iR).Orbit(iT).Period(iP).mEigR(iW) = mean(real(theseEmax));  % mean of real part
                        RcrPtData(iR).Orbit(iT).Period(iP).mEigI(iW) = mean(imag(theseEmax));  % mean of imaginary part               
                        RcrPtData(iR).Orbit(iT).Period(iP).mEigAbs(iW) = mean(abs(theseEmax));  %abs of all Eigs

                   end
               end
%                figure
%                plot(Data(iR).stimbins(RcrPtData(iR).winmids),RcrPtData(iR).Orbit(iT).Period(iP).mRotPeriod,'o','Color',cmap(iP,:))
%                %plot(Data(iR).stimbins(RcrPtData(iR).winmids),Orbit(iT).Period(iP).Tstable,'o','Color',cmap(iP,:))
%                %plot(Data(iR).stimbins(RcrPtData(iR).winmids),Orbit(iT).Period(iP).mEigAbs,'o','Color',cmap(iP,:))
%                ylabel('Rotation period')
           end
           RcrPtData(iR).Orbit(iT).frstRcr = min([RcrPtData(iR).Orbit(iT).Period(:).frstRcr]); % first point that recurs in all Periods

       end
       %% non-periodic and non-recurrent 
       % for the non-period points - get local linear dynamics
       [RcrPtData(iR).Orbit(iT).NoPeriod.Emax,RcrPtData(iR).Orbit(iT).NoPeriod.allEgs,RcrPtData(iR).Orbit(iT).NoPeriod.RotPeriod,RcrPtData(iR).Orbit(iT).NoPeriod.Npts,...
           RcrPtData(iR).Orbit(iT).NoPeriod.Cond,RcrPtData(iR).Orbit(iT).NoPeriod.R2] = LocalLinearDynamics(RcrPtData(iR).Orbit(iT).NoPeriod.IDs,...
                            RcrPtData(iR).currPCs,RcrPtData(iR).ix0,RcrPtData(iR).Theta(iT),pars.NeighbourD,Data(iR).GaussQt);
       
      % for the non-recurrent points - get local linear dynamics                  
       [RcrPtData(iR).Orbit(iT).NoRecurrence.Emax,RcrPtData(iR).Orbit(iT).NoRecurrence.allEgs,RcrPtData(iR).Orbit(iT).NoRecurrence.RotPeriod,RcrPtData(iR).Orbit(iT).NoRecurrence.Npts,...
           RcrPtData(iR).Orbit(iT).NoRecurrence.Cond,RcrPtData(iR).Orbit(iT).NoRecurrence.R2] = LocalLinearDynamics(RcrPtData(iR).Orbit(iT).NoRecurrence.IDs,...
                            RcrPtData(iR).currPCs,RcrPtData(iR).ix0,RcrPtData(iR).Theta(iT),pars.NeighbourD,Data(iR).GaussQt);
                       
       % average in windows 
        for iW = 1:numel(RcrPtData(iR).winstrt)
           ixPts = find(RcrPtData(iR).Orbit(iT).NoPeriod.IDs > RcrPtData(iR).winstrt(iW)-RcrPtData(iR).ix0 & RcrPtData(iR).Orbit(iT).NoPeriod.IDs <= RcrPtData(iR).winstrt(iW)+RcrPtData(iR).winDT-RcrPtData(iR).ix0);
           RcrPtData(iR).Orbit(iT).NoPeriod.nWinPoints(iW) = numel(ixPts);  % how many points in this window for non-period
           if isempty(ixPts)
               RcrPtData(iR).Orbit(iT).NoPeriod.mDelay(iW) = nan;
               RcrPtData(iR).Orbit(iT).NoPeriod.nWinUsedPoints = 0;               
               RcrPtData(iR).Orbit(iT).NoPeriod.mRotPeriod(iW) = nan;
               RcrPtData(iR).Orbit(iT).NoPeriod.Tstable(iW) = nan;
               RcrPtData(iR).Orbit(iT).NoPeriod.mEigR(iW) = nan;
               RcrPtData(iR).Orbit(iT).NoPeriod.mEigI(iW) = nan;
               RcrPtData(iR).Orbit(iT).NoPeriod.mEigAbs(iW) = nan;
           else
               RcrPtData(iR).Orbit(iT).NoPeriod.mDelay(iW) = nanmean(RcrPtData(iR).Tdelay(ixPts,iT));
                theseN = RcrPtData(iR).Orbit(iT).NoPeriod.Npts(ixPts);
                usedPts = ixPts(theseN >= pars.minNeighbours);
                RcrPtData(iR).Orbit(iT).NoPeriod.nWinUsedPoints(iW) = numel(usedPts);
                theseRotations = RcrPtData(iR).Orbit(iT).NoPeriod.RotPeriod(usedPts);
                RcrPtData(iR).Orbit(iT).NoPeriod.mRotPeriod(iW) = mean(theseRotations(theseRotations > 0));  % local estimate of current rotational period
                RcrPtData(iR).Orbit(iT).NoPeriod.Tstable(iW) = sum(theseRotations == 0) * Data(iR).GaussQt;
                theseEmax = RcrPtData(iR).Orbit(iT).NoPeriod.Emax(usedPts);
%                 RcrPtData(iR).Orbit(iT).NoPeriod.mEigR(iW) = mean(real(theseEmax(theseRotations > 0)));  % mean of real part, of Eigs of rotation
%                 RcrPtData(iR).Orbit(iT).NoPeriod.mEigI(iW) = mean(imag(theseEmax(theseRotations > 0)));  % mean of imaginary part, of Eigs of rotation                                              
                RcrPtData(iR).Orbit(iT).NoPeriod.mEigR(iW) = mean(real(theseEmax));  % mean of real part
                RcrPtData(iR).Orbit(iT).NoPeriod.mEigI(iW) = mean(imag(theseEmax));  % mean of imaginary part                                             
                RcrPtData(iR).Orbit(iT).NoPeriod.mEigAbs(iW) = mean(abs(theseEmax));  %abs of all Eigs

           end
           
           ixPts = find(RcrPtData(iR).Orbit(iT).NoRecurrence.IDs > RcrPtData(iR).winstrt(iW)-RcrPtData(iR).ix0 & RcrPtData(iR).Orbit(iT).NoRecurrence.IDs <= RcrPtData(iR).winstrt(iW)+RcrPtData(iR).winDT-RcrPtData(iR).ix0);
           RcrPtData(iR).Orbit(iT).NoRecurrence.nWinPoints(iW) = numel(ixPts);  % how many points in this window for non-period
           if isempty(ixPts)
               RcrPtData(iR).Orbit(iT).NoRecurrence.nWinUsedPoints(iW) = 0;     
               RcrPtData(iR).Orbit(iT).NoRecurrence.mRotPeriod(iW) = nan;
               RcrPtData(iR).Orbit(iT).NoRecurrence.Tstable(iW) = nan;
               RcrPtData(iR).Orbit(iT).NoRecurrence.mEigR(iW) = nan;
               RcrPtData(iR).Orbit(iT).NoRecurrence.mEigI(iW) = nan;
               RcrPtData(iR).Orbit(iT).NoRecurrence.mEigAbs(iW) = nan;
           else
                theseN = RcrPtData(iR).Orbit(iT).NoRecurrence.Npts(ixPts);
                usedPts = ixPts(theseN >= pars.minNeighbours);
                RcrPtData(iR).Orbit(iT).NoRecurrence.nWinUsedPoints(iW) = numel(usedPts);
                theseRotations = RcrPtData(iR).Orbit(iT).NoRecurrence.RotPeriod(usedPts);
                RcrPtData(iR).Orbit(iT).NoRecurrence.mRotPeriod(iW) = mean(theseRotations(theseRotations > 0));  % local estimate of current rotational period
                RcrPtData(iR).Orbit(iT).NoRecurrence.Tstable(iW) = sum(theseRotations == 0) * Data(iR).GaussQt;
                theseEmax = RcrPtData(iR).Orbit(iT).NoRecurrence.Emax(usedPts);
%                 RcrPtData(iR).Orbit(iT).NoRecurrence.mEigR(iW) = mean(real(theseEmax(theseRotations > 0)));  % mean of real part, of Eigs of rotation
%                 RcrPtData(iR).Orbit(iT).NoRecurrence.mEigI(iW) = mean(imag(theseEmax(theseRotations > 0)));  % mean of imaginary part, of Eigs of rotation                                              
                RcrPtData(iR).Orbit(iT).NoRecurrence.mEigR(iW) = mean(real(theseEmax));  % mean of real part
                RcrPtData(iR).Orbit(iT).NoRecurrence.mEigI(iW) = mean(imag(theseEmax));  % mean of imaginary part                                          
                RcrPtData(iR).Orbit(iT).NoRecurrence.mEigAbs(iW) = mean(abs(theseEmax));  %abs of all Eigs

           end
           
           
       end      
   end
   % keyboard
%    % recurrence points of set (m,e) %    
%     figure 
%     plot(RcrPtData(iR).currPCs(RcrPtData(iR).Orbit(iT).Period(1).IDs+ix0,1),RcrPtData(iR).currPCs(RcrPtData(iR).Orbit(iT).Period(1).IDs+ix0,2),'k.')
  
    

end

save([prefix fname '_RecurrencePointsData_FilterPts'],'RcrPtData','pars'); % ,'-v7.3');
% save([fname '_RecurrencePointsData_Kmeans'],'RcrPtData','pars'); % ,'-v7.3');

clear('RcrPtData','pars');