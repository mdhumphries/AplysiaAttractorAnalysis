%% are sequential attractors actually the same?
%% could there be more than one stable limit cycle/spiral in system? 
%
% NB this also checking across cases where there is clear transition in the
% system while it's running!

clear all; close all
addpath ../Functions/

% analysis parameters
spars.stimstart = 30;  % stimulation time
spars.VmaxD = 3;  % number of dimensions to visualise
spars.minDens = 0; % min proportion of recurrence points with period
spars.ixT = 3; % which threshold to use
spars.minRecur = 90; % minimum percentage of recurrence points per time-window
spars.manifoldT = 5; % distance to be on same manifold
spars.VarExplained = 0.8; % dimensions to retain
spars.divPeriod = 0;  % radius on checking min distance as division of oscillation period (0 = no radius; 1 = period; 2 = half period etc)
spars.nShuffles = 100;
% 
fname_Pre = '';
load ../da01_DataProperties_FunctionAndWindowSize FileTable
nfiles = numel(FileTable);
load RecurrenceStats RcrStats
fnames = {'da01','da02','da03'};

ixExample = 1;

% uses Parallel Distributed Toolbox if available
feature('numCores');
if isempty(gcp('nocreate'))
    parpool('local',10);
end

%% compare attractors
Tctr = 1;
for iR = 1:nfiles
    tic
    All = []; conSDF = []; pre = 0;
   for iF = 1:numel(fnames) 
       % load stuff and store across all 3 programs
       load(['../' fname_Pre fnames{iF} '_DataProperties_FunctionAndWindowSize'],'SDF','Data') % to get density functions and basic data
       load([fname_Pre fnames{iF} '_RecurrencePointsData_FilterPts']); % to get recurrence: density per window, and winmids etc
       All(iF).n = numel(Data(iR).rates); % number of neurons
       
       %  store for all comparisons
       All(iF).Data = Data(iR);
       All(iF).winDT = RcrPtData(iR).winDT; All(iF).winmids = RcrPtData(iR).winmids; All(iF).winstrt = RcrPtData(iR).winstrt;  
       All(iF).PrctWnone = 100*RcrPtData(iR).Orbit(spars.ixT).NoRecurrence.nWinPoints/RcrPtData(iR).winDT;  % density of any recurrence point per window
       All(iF).ixWpts = find(All(iF).PrctWnone <= 100 - spars.minRecur);  % find all windows with high recurrence density
       
       % find set of admissible time-points in each attractor
        ixbreaks = [0 find(diff(All(iF).ixWpts) > 1) numel(All(iF).ixWpts)]; % discontinuties
        ixAdPts = [];
        % get all unbroken sequences of recurrence points
        for iB = 1:numel(ixbreaks)-1
            strts = All(iF).winstrt(All(iF).ixWpts(ixbreaks(iB)+1)); % start of first window in this unbroken sequence
            ends = All(iF).winstrt(All(iF).ixWpts(ixbreaks(iB+1)))+All(iF).winDT; % end of last window in this unbroken sequence
            ixAdPts = [ixAdPts strts:ends];
        end
        All(iF).ixAdPts = ixAdPts;
       % keyboard
       
       % concatenate SDFs for global PCA
       All(iF).ixts = find(Data(iR).bins < spars.stimstart);
       All(iF).ixStim = find(Data(iR).bins >= spars.stimstart);
       conSDF = [conSDF; SDF(iR).spkfcn(All(iF).ixStim,:)];
       Prep(iR).PIdx(iF).ts = pre + (1:numel(All(iF).ixStim)); pre = Prep(iR).PIdx(iF).ts(end); % track index of each program in concatenated PCA
       
   end
   % compute PCA
   [Prep(iR).PCA.coeffs, ~, Prep(iR).PCA.eigvalues] = pca(conSDF);
   Prep(iR).PCA.prctVar = cumsum(Prep(iR).PCA.eigvalues) / sum(Prep(iR).PCA.eigvalues);  % compute variance explained
   Prep(iR).PCA.nPCs = sum(Prep(iR).PCA.prctVar <= pars.VarExplained);  % parameter from recurrence analysis saved data

   % project unnormalised axes
   for iP = 1:Prep(iR).PCA.nPCs 
        Prep(iR).PCprojU(:,iP) = sum(bsxfun(@times,conSDF,Prep(iR).PCA.coeffs(:,iP)'),2);  % raw
   end

  
        
   % (1) plot PCA

    projmap = cbrewer('qual','Set2',3);
    % projmap = cbrewer('seq','Reds',3);
    hS = figure; 
    for iF = 1:numel(fnames)
        plot3(Prep(iR).PCprojU(Prep(iR).PIdx(iF).ts,1),Prep(iR).PCprojU(Prep(iR).PIdx(iF).ts,2),Prep(iR).PCprojU(Prep(iR).PIdx(iF).ts,3),'Color',projmap(iF,:)); hold on
    end    
    figure(hS); grid on; title(['Recording' num2str(iR) ' on single axes set'])
    
    % fix storage of plotting - for paper figures
    for iF = 1:numel(fnames)
        % store for plotting in example figure
        % Viz(iR,iF).ProjSponOnda01 = All(iF).SponOnda01(:,1:3); % for plotting only, so just 3-D!
        Viz(iR,iF).CommonProj = Prep(iR).PCprojU(Prep(iR).PIdx(iF).ts,1:spars.VmaxD);
    end
    
    %  get min. distance between points on
    % attractor; compare to min distance to shuffle control
    % IF all 3 are actually moving on the attractor sufficiently to
    % estimate this!
    
    if all([RcrStats(iR,:).densAllPeriod] > spars.minDens)
        % initialise storage
        tmax = numel(All(1).Data.stimbins);
        Prep(iR).Dmin = zeros(numel(fnames),numel(fnames),tmax) + nan;  % set as Nans = non admissible points
        DminCtrl = zeros(numel(fnames),numel(fnames),spars.nShuffles,tmax) + nan; 

        ctr = 1;  % just for plotting
        figure
        subplot(2,2,4); title(['Recording #' num2str(iR)]) 
        for i = 1:numel(fnames)
            % set of admissible recurrent points on target
            ixP = Prep(iR).PIdx(i).ts; ixP = ixP(All(i).ixAdPts);  % set of admissible points within projection
            targetPCs = Prep(iR).PCprojU(ixP,1:Prep(iR).PCA.nPCs);  % get subset of concatenated PCs that correspond to this progam
           
            for j = i+1:numel(fnames)  
                tic
                
                mPeriod = (RcrStats(iR,i).mRot + RcrStats(iR,j).mRot) / 2; % mean of period
                radius = mPeriod/spars.divPeriod; % radius to check for distance

                % set of admissible time-points in source program
                ixP = Prep(iR).PIdx(j).ts; ixP = ixP(All(j).ixAdPts);  % set of admissible points within projection
                srcPCs = Prep(iR).PCprojU(ixP,1:Prep(iR).PCA.nPCs); % get subset of concatenated PCs that correspond to this progam
                srcts = All(j).Data.stimbins(All(j).ixAdPts); % corresponding times of admissible points within the progam
                targetts = All(i).Data.stimbins(All(i).ixAdPts); % reinstate the corresponding times of target admissible points
                
                % find common start point of target and source programs
                cmnStrt = max(srcts(1),targetts(1));
                ixSrc = find(srcts >= cmnStrt,1);
                ixTgt = find(targetts >= cmnStrt,1);
                srcts(srcts < cmnStrt) = []; targetts(targetts < cmnStrt) = []; % remove all time points before joint convergence...
                
                % for every point in source attractor, find min distance to
                % any point in target attractor, within period...
                ix = zeros(numel(srcts),1);
                for iN = 1:numel(srcts)
                    pt0 = srcPCs(iN+ixSrc-1,:);  % current x_i
                    ptts = srcts(iN);  % current time in src program
                    ixts = targetts < ptts + radius & targetts > ptts - radius;  % find all times in target within defined period
                    D = sqrt(sum(bsxfun(@minus,targetPCs(ixts,:),pt0).^2,2)); % Euclidean distance to all points within time radius
                    ix(iN) = find(All(j).Data.stimbins == ptts);  % store index for use in shuffled batch
                    Prep(iR).Dmin(i,j,ix(iN)) = min(D); % above diagonal: source -> target
                   
                end
                
                Dalltemp = zeros(spars.nShuffles,numel(srcts));
                parfor iB = 1:spars.nShuffles
                    % get shuffled target
                    ixR = randperm(All(i).n);
                    ControlStimOnda01 = zeros(numel(All(i).ixStim),Prep(iR).PCA.nPCs);
                    for iP = 1:Prep(iR).PCA.nPCs
                        ControlStimOnda01(:,iP) = sum(bsxfun(@times,conSDF(Prep(iR).PIdx(i).ts,ixR),Prep(iR).PCA.coeffs(:,iP)'),2); % shuffle neuron indices, and project onto da01 axes
                    end
                    ctrltargetPCs = ControlStimOnda01(All(i).ixAdPts,1:Prep(iR).PCA.nPCs); % control source: itself, shuffled
                    Dtemp = zeros(1,numel(srcts));
                    for iN = 1:numel(srcts)
                        pt0 = srcPCs(iN+ixSrc-1,:);  % current x_i
                        ptts = srcts(iN);  % current time in src program
                        ixts = targetts < ptts + radius & targetts > ptts - radius;  % find all times in target within defined period
                        Dctrl = sqrt(sum(bsxfun(@minus,ctrltargetPCs(ixts,:),pt0).^2,2)); % Euclidean distance to all points within time radius
                        % DminCtrl(i,j,iB,ix(iN)) = min(Dctrl); % below diagonal: source -> target*
                        Dtemp(iN) = min(Dctrl); % below diagonal: source -> target* 
                    end
                    Dalltemp(iB,:) = Dtemp;  % dealing with "slicing" issues
                end
                DminCtrl(i,j,:,ix) = Dalltemp;  % map Dtemp to correct indices...
                
                % and vice versa: from target to source
                ix = zeros(numel(targetts),1);
                for iN = 1:numel(targetts)
                    % current time....
                    pt0 = targetPCs(iN+ixTgt-1,:);  % current x_i
                    ptts = targetts(iN);
                    ixts = srcts < ptts + radius & srcts > ptts - radius;
                    D = sqrt(sum(bsxfun(@minus,srcPCs(ixts,:),pt0).^2,2)); % Euclidean distance to all points within time radius
                    ix(iN) = find(All(i).Data.stimbins == ptts);  % store index for use in shuffled batch
                    Prep(iR).Dmin(j,i,ix(iN)) = min(D); % below diagonal: target -> source
                end
                
                Dalltemp = zeros(spars.nShuffles,numel(targetts));
                parfor iB = 1:spars.nShuffles
                    % get shuffled source
                    ixR = randperm(All(j).n);
                    ControlStimOnda01 = zeros(numel(All(j).ixStim),Prep(iR).PCA.nPCs);                    
                    for iP = 1:Prep(iR).PCA.nPCs
                        ControlStimOnda01(:,iP) = sum(bsxfun(@times,conSDF(Prep(iR).PIdx(j).ts,ixR),Prep(iR).PCA.coeffs(:,iP)'),2); % shuffle neuron indices, and project onto da01 axes
                    end
                    ctrlsrcPCs = ControlStimOnda01(All(j).ixAdPts,1:Prep(iR).PCA.nPCs); % control source: itself, shuffled
                    Dtemp = zeros(1,numel(targetts));
                    for iN = 1:numel(targetts)
                        pt0 = targetPCs(iN+ixTgt-1,:);  % current x_i
                        ptts = targetts(iN);
                        ixts = srcts < ptts + radius & srcts > ptts - radius;
                        Dctrl = sqrt(sum(bsxfun(@minus,ctrlsrcPCs(ixts,:),pt0).^2,2)); % Euclidean distance to all points within time radius
                        %DminCtrl(j,i,iB,ix(iN)) = min(Dctrl); % below diagonal: target -> source*
                        Dtemp(iN) = min(Dctrl); % below diagonal: source -> target* 
                    end
                    Dalltemp(iB,:) = Dtemp;  % dealing with "slicing" issues
                end
                DminCtrl(j,i,:,ix) = Dalltemp;  % map Dtemp to correct indices...
                
                % keyboard
                dd1 = squeeze(Prep(iR).Dmin(i,j,:));
                Prep(iR).MedDmin(i,j) = nanmedian(dd1);
                Prep(iR).MaxDmin(i,j) = nanmax(dd1);
                % propotion of distances that are within tolerance of being
                % on same manifold
                Prep(iR).prctDmin(i,j) = sum(dd1 <= spars.manifoldT)./numel(dd1);
               
                dd2 = squeeze(Prep(iR).Dmin(j,i,:));
                Prep(iR).MedDmin(j,i) = nanmedian(dd2);
                Prep(iR).MaxDmin(j,i) = nanmax(dd2);
                % propotion of distances that are within tolerance of being
                % on same manifold
                Prep(iR).prctDmin(j,i) = sum(dd2 <= spars.manifoldT)./numel(dd2);
             
                
                Prep(iR).HaussD(i,j) = max(nanmax(dd1),nanmax(dd2)); % maximum of minium distances between sets
                
                
                dctrl1 = squeeze(DminCtrl(i,j,:,:));
                Prep(iR).MedDminCtrl(i,j,:) = nanmedian(dctrl1,2);
                Prep(iR).MaxDminCtrl(i,j,:) = nanmax(dctrl1,[],2);
                Prep(iR).prctDminCtrl(i,j,:) = sum(dctrl1 <= spars.manifoldT,2)./numel(dctrl1);
                
                dctrl2 = squeeze(DminCtrl(j,i,:,:));
                Prep(iR).MedDminCtrl(j,i,:) = nanmedian(dctrl2,2);
                Prep(iR).MaxDminCtrl(j,i,:) = nanmax(dctrl2,[],2);
                Prep(iR).prctDminCtrl(j,i,:) = sum(dctrl2 <= spars.manifoldT,2)./numel(dctrl2);
                
                Prep(iR).HaussDCtrl(i,j,:) = max([nanmax(dctrl1,[],2),nanmax(dctrl2,[],2)],[],2); % maximum of minium distances between sets

                % plot of whole distribution of min distance
                [f,x] = ecdf(dd1);
                [f2,x2] = ecdf(dd2);
                [fctrl,xctrl] = ecdf(nanmean(dctrl1));   % mean over shuffles
                [fctrl2,xctrl2] = ecdf(nanmean(dctrl2));
               
                subplot(2,2,ctr),stairs(x,f,'r'); hold on
                stairs(x2,f2,'k')
                stairs(xctrl,fctrl,'r--')
                stairs(xctrl2,fctrl2,'k--')
                ctr = ctr+1;
            end
        end
    end
   toc 
end


%% save stuff:
save([fname_Pre 'RecurrenceManifoldStats_CommonAxes_AllPreps'],'Prep','Viz','spars') 


%% SUMMARISE

% cmap = cbrewer('qual','Paired',nfiles); 
cmap = brewermap(10,'Paired'); 
% summary plot: 
% (1) color code by prep to show consistency within prep
% (2) link within prep to draw eye to groups
% (3) markersize according to stim sequence: no patterns
% (4) markersize according to "sensitised"

Msize = [5 7 9];  % stim order
% Msize = [5 15 10];  % does sensitised prep show clear effect?

hDmin = figure; hold on
hDHD = figure; hold on

sym = 'o';

for iR = 1:nfiles
    % ID P10 preps using squares
%     if any(iR == P10preps')
%         sym = 's'; 
%     else
%         sym = 'o';
%     end
    % (1) plot all r for PC correlations (for each PC, and together
         
%      % plot min distance
%      if ~isempty(Prep(iR).prctDmin)
%          prctDdata = 100*Prep(iR).prctDmin; prctDcontrol = 100*Prep(iR).prctDminCtrl;
%          figure(hDmin)
%          plot(prctDcontrol(prctDcontrol > 0),prctDdata(prctDdata > 0),'o-',...
%             'Color',cmap(iR,:),'MarkerSize',8,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
%          text(prctDcontrol(2,3),prctDdata(2,3),num2str(iR),'Fontsize',8,'Color',[1 1 1])
%      end
     
     if ~isempty(Prep(iR).HaussD)
         ctrl = mean(Prep(iR).HaussDCtrl,3); ctrl = ctrl(ctrl > 0); 
         sem = std(Prep(iR).HaussDCtrl,0,3) / sqrt(spars.nShuffles); sem = sem(sem > 0);
         data = Prep(iR).HaussD(Prep(iR).HaussD > 0);
         figure(hDHD)
         plot(ctrl,data,'o',...
            'Color',cmap(iR,:),'MarkerSize',8,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
         line([ctrl-2*sem ctrl+2*sem]',[data data]','Linewidth',0.5,'Color',cmap(iR,:))
         % text(prctDcontrol(2,3),prctDdata(2,3),num2str(iR),'Fontsize',8,'Color',[1 1 1])
     end
end

figure(hDmin)
title('Percentage points on same manifold')
ylabel('Data (%)')
xlabel('Control (%)')
line([0 100],[0 100],'Color',[0.6 0.6 0.6])
axis([0 100 0 100])

figure(hDHD)
title('Haussdorf distance')
ylabel('Data (spikes/s)')
xlabel('Control (spikes/s)')
axis square
axis([0 30 0 30])
line([0 30],[0 30],'Color',[0.6 0.6 0.6])

