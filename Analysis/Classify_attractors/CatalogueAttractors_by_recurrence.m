% script to catalogue attractors
% does periodic orbit recurrence ~ Lathrop & Kostelich (1989)
% checks recurrence plots - cf Marwan et al 2007
% ADD osc and burst to this: 

clear all; close all

addpath ../../Functions/

fname = 'da03';
prefix = [];


% get PCA projections
load(['../' prefix fname '_DataProperties_FunctionAndWindowSize'],'PCAdata','Data','SDF','FileTable')
nfiles = numel(PCAdata);

% % get cell types
% load(['../Ensembles/' prefix fname '_Analyses_Neurons_and_Groups'],'groupdata')
% load(['../Ensembles/' prefix fname '_Ensemble_Types'],'AcorrTypes')

% display
hstmap = cbrewer('seq','OrRd',10);

% parameters
stimstart = 30; % (s)
T0 = 35;  % start counting after stim (s)
maxWindow = 10;  % stop before end of time-series (s)
VarExplained = 0.80;    % how many PCs to keep

prctTheta = [1:1:10]; %  % of distance distribution to use as threshold; max distance to check (spikes/s)

% maxTheta = 10; % max distance to check (spikes/s)
% stepTheta = 0.5;
% Theta = stepTheta:stepTheta:maxTheta;

bstep = 1;  % s
bins = 0:bstep:50;  % recurrence time bins (s);

tmin = 5;  % s; minimum recurrence time to count as a peak
Nmin = 100;      % minimum number of points to use that peak
Kmax = 5;   % maximum number of k-means to try

minbin = bins(find(bins <= tmin,1,'last')); % lowest retained bin in histogram

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
%     % select subset of groups for PCA: specify here what to use, then save
%     % with appropriate names....
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
    ixSpon = find(Data(iR).bins < stimstart);
    ixStim = find(Data(iR).bins(1:Tend) >= stimstart);
    
%     % get PCA of stim period
%     % [PCAsubset(iR).coeffs, PCAsubset(iR).scores, PCAsubset(iR).eigvalues] = princomp(SDF(iR).spkfcn(ixStim,PCAsubset(iR).neuronIDs));
%     [PCAsubset(iR).coeffs, PCAsubset(iR).scores, PCAsubset(iR).eigvalues] = princomp(SDF(iR).spkfcn(:,PCAsubset(iR).neuronIDs));
%      PCAsubset(iR).prctVar = cumsum(PCAsubset(iR).eigvalues) / sum(PCAsubset(iR).eigvalues);
%      PCAsubset(iR).maxD = sum(PCAsubset(iR).prctVar <= VarExplained);
    
    % pick the subset of dimensions for whole thing 
    maxD = sum(PCAdata(iR).prctVar <= VarExplained);
    
    % get subset of PCs to act as state-space
    currPCs = PCAdata(iR).PCproj(:,1:maxD);  % from stimulation    
    % currPCs = PCAdata(iR).allPCproj(:,1:maxD);  % from start
    % currPCs = PCAsubset(iR).scores(:,1:PCAsubset(iR).maxD); % from start for oscillators 
        
    % get Theta
    dPs = pdist(currPCs,'euclidean');  % all distances
    Theta = prctile(dPs,prctTheta);
    
    % for every time-point, find nearest time-point within threshold
    ix0 = 1; % initial time point is start of recording
    % find(Data(iR).stimbins == T0); % initial time-point after stimulation
    ixEnd = find(Data(iR).stimbins == Data(iR).stimbins(end)-maxWindow);
    nTs = ixEnd - ix0;  % loop from start to before end
    tic
    for iN = 1:nTs
        pt0 = currPCs(ix0+iN,:);  % current x_i
        currixs = ix0+iN+1:numel(Data(iR).stimbins);
        D = sqrt(sum(bsxfun(@minus,currPCs(currixs,:),pt0).^2,2)); % Euclidean distance to all points ahead in time
        for iT = 1:numel(Theta)
            allpts = find(D <= Theta(iT));  % every point within this threshold
        
            % remove those immediately after T that are within theta -
            % short time-scale noise...
            
%             % Need the below for looking backward and forward in time
%             ix0pt = find(allpts == ix0);  % T0 in this subset
%             dpts = diff(allpts);
%             ulim = find(dpts(ix0pt+1:end) > 1,1,'first'); if isempty(ulim) ulim = numel(allpts); end
%             llim = find(dpts(1:ix0pt-1) > 1,1,'last'); if isempty(llim) llim = 0; end
%             allpts(llim+1:ulim) = [];
            
            dpts = diff(allpts);
            ulim = find(dpts > 1,1,'first'); if isempty(ulim) ulim = numel(allpts); end
            allpts(1:ulim) = [];
            
            dly = min(Data(iR).stimbins(currixs(allpts))) - Data(iR).stimbins(ix0+iN); % keep shortest delay
            if dly                 
                % store delay: any point with delay is a (m,e) point by
                % definition
                Tdelay(iN,iT) = dly;
            else
                % nothing
                Tdelay(iN,iT) = nan;
            end
            
        end
        
    end
    toc
    %% get (m,theta) points: find period 1, 2 etc of recurrence points
    tic
    hstRecurrence = zeros(numel(Theta),numel(bins));
    for iT = 1:numel(Theta)
         hst = histc(Tdelay(:,iT),bins);
         hstRecurrence(iT,:) = hst / sum(hst); % approx PDF
         % find histogram peaks to get points in each period
         % (1) histogram edge finding
         % keyboard
  
        filledbins = [find(hstRecurrence(iT,minbin:end) > 0)];
        if ~isempty(filledbins)
            length = diff(filledbins);
            hstedges = [[filledbins(1); filledbins(find(length > 1)+1)'] [filledbins(length > 1)'; filledbins(end)]] + minbin-1;
            npeaks = size(hstedges,1);
            ctr = 0;
            for i=1:npeaks
                % get all points within that histogram peak
                ixts = find(Tdelay(:,iT) >= bins(hstedges(i,1)) & Tdelay(:,iT) <= bins(hstedges(i,2))+1);

                if numel(ixts) >= Nmin
                    % if enough points, keep as candidate Period N set
                    ts = Tdelay(ixts,iT);
                    ctr = ctr+1;
                    Orbit(iT).Period(ctr).IDs = ixts; % get IDs of these cells
                    Orbit(iT).Period(ctr).mTime = median(ts);
                    Orbit(iT).Period(ctr).nPoints = numel(ts); 
                end

            end
            Orbit(iT).nPeriods = ctr;
            
         else
             Orbit(iT).nPeriods = 0;
        end
    end
    toc
    
%          %% (2) kmeans - wroks very well, but MUCH slower, and more
%          complex to explain
%          X = Tdelay(~isnan(Tdelay(:,iT)),iT);  % all recurrence times 
%          X(X < tmin) = []; % remove all that likely tangenital motion 
%          if ~isempty(X)
%              Idx = zeros(numel(X),Kmax);
%              for k= 1:Kmax
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
%             ctr = 0; 
%             for iG = 1:nGrps
%                 % check if group later than close trajectory
%                 thisG = find(grpscon(:,2) == grpIDs(iG));
%                 ts = X(thisG);
%                 % keep if all after min time AND enough entries
%                 if numel(ts) >= Nmin
%                     ctr = ctr+1;
%                     ixNotNan = find(~isnan(Tdelay(:,iT)));
%                     Orbit(iT).Period(ctr).IDs = ixNotNan(thisG); % get IDs of these cells
%                     Orbit(iT).Period(ctr).mTime = median(ts);
%                     Orbit(iT).Period(ctr).nPoints = numel(ts); 
%                 end
% 
%             end
%             if ctr > 1
%                 % sort period into time order (Period 1, Period 2 etc..)
%                [~,iP] = sort([Orbit(iT).Period(:).mTime]);
%                Orbit(iT).Period = Orbit(iT).Period(iP);  
%             end
%             Orbit(iT).nPeriods = ctr;
%             
%          else
%              Orbit(iT).nPeriods = 0;
%          end
%    end 
%    toc
    
    % keyboard
    
    % plot histogram, get number of periods etc

    figure
    imagesc(bins,Theta,hstRecurrence); 
    colormap(hstmap);
    title(['No. PCs = ' num2str(maxD)])
    % title(['No. PCs = ' num2str(PCAsubset(iR).maxD)])
  
    xlabel('Time delay (s)')
    ylabel('Threshold (spikes/s)')
    for iT = 1:numel(Theta)
        text(52, Theta(iT),num2str(Orbit(iT).nPeriods),'fontsize',6)
    end
    
    figname = ['Prep_' num2str(iR) '_' num2str(fname) '_' FileTable{iR}(1:strfind(FileTable{iR},'.mat')-1) '_RecurHst'];
    exportPPTfig(gcf,figname,[5 5 8 6])
   
   
    %% recurrence plot 
    N = size(currPCs,1);
    k = round(tmin / Data(iR).GaussQt); % number of bins with local dynamics only

    for iT=1:numel(Theta)
         % find recurrence points
        RecurPlot(iT).Rp = squareform(dPs <= Theta(iT));  % automatically puts 0s on diagonal, so no LOI
    end
   
   
    % pick one to plot and analyse
    ixPlot = numel(Theta);
    figure
    imagesc(Data(iR).stimbins,Data(iR).stimbins,RecurPlot(ixPlot).Rp)
    % imagesc(Data(iR).bins,Data(iR).bins,RecurPlot(ixPlot).Rp)

    xlabel('Time (s)')
    ylabel('Time (s)')
    title(['Theta = ' num2str(Theta(ixPlot)) ' spikes/s'])
    colormap(gray)
    figname = ['Prep_' num2str(iR) '_' num2str(fname) '_' FileTable{iR}(1:strfind(FileTable{iR},'.mat')-1) '_RecurPlot'];
    exportPPTfig(gcf,figname,[5 5 8 8])

    % keyboard
   
   
    

end