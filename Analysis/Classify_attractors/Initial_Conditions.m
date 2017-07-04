%% check initial conditions at stimulation onset
% (1) spontaneous state: distance between
% (2) which neurons were recruited

clear all; close all
 
addpath ../../Functions/


spikepath = '../../Data/Spikes/';

pars.stimstart = 30; % 30 s into recording
pars.stimend = 32.5;
pars.prestim = 27.5;

fnames = {'da01','da02','da03'};  %
nStims = numel(fnames);

%% load and get rates


% for each program: find its spike file (loop over da01 etc) - load
% FileTable!
%
% get the spikes
% compute rates during stimulation

SpkData = struct([]);  

for iStim = 1:nStims
    load(['../' fnames{iStim} '_DataProperties_FunctionAndWindowSize.mat'],'FileTable')
    nPreps = numel(FileTable);

    for iPrep = 1:nPreps
        load([spikepath FileTable{iPrep}]);
        
        SpkData(iPrep,iStim).n_trains = numel(unique(spks(:,1)));  % count number in the retained period

        SpkData(iPrep,iStim).rates = zeros(1,SpkData(iPrep,iStim).n_trains);
        for iN = 1:SpkData(iPrep,iStim).n_trains
            ts = spks(spks(:,1)==iN,2);
            SpkData(iPrep,iStim).rates.postStim(iN) = sum(ts >= pars.stimstart & ts <= pars.stimend); 
            SpkData(iPrep,iStim).rates.preStim(iN) = sum(ts >= pars.prestim & ts <= pars.stimstart); 
        end
    end
end

%% visualise
rcmap = brewermap(10,'*OrRd');
dcmap = brewermap(10,'RdBu');
for iPrep = 1:nPreps
    for iStim = 1:nStims
        Viz(iPrep).matPre(:,iStim) = SpkData(iPrep,iStim).rates.preStim;
        Viz(iPrep).matPost(:,iStim) = SpkData(iPrep,iStim).rates.postStim;
        Viz(iPrep).matPrePost(:,iStim) = SpkData(iPrep,iStim).rates.postStim - SpkData(iPrep,iStim).rates.preStim;
    end
    
%     figure(iPrep)
%     subplot(131), imagesc(1:3,1:SpkData(iPrep,iStim).n_trains,Viz(iPrep).matPre); colormap(rcmap); freezeColors;
%     subplot(132), imagesc(1:3,1:SpkData(iPrep,iStim).n_trains,Viz(iPrep).matPost); colormap(rcmap); freezeColors;
%     subplot(133), imagesc(1:3,1:SpkData(iPrep,iStim).n_trains,Viz(iPrep).matPrePost); colormap(dcmap); freezeColors;
   
    figure(iPrep)
    % linked scatter plot of Post-Pre change for every neuron across 3
    % stims: change in sign clearly show variable initial state and
    % response!!
    line([0.5 3.5],[0 0],'Color',[0 0 0],'Linewidth',0.5); hold on
    LinkedUnivariateScatterPlots(gca,1:3,Viz(iPrep).matPrePost,[0.7 0.7 0.7]);
    % LinkedUnivariateScatterPlots(gca,1:3,Viz(iPrep).matPre,[0.7 0.7 0.7]);
    
    % count proportion of neurons with at least one opposite sign change in
    % pre:post.... 
   
    pos = sum(Viz(iPrep).matPrePost > 0,2);
    neg = sum(Viz(iPrep).matPrePost < 0,2);
    IC(iPrep).prepost = sum(pos > 0 & neg > 0) ./ SpkData(iPrep,iStim).n_trains;
end

%% load common-axes PCA, and get distances
load('RecurrenceManifoldStats_CommonAxes_AllPreps','Prep','spars') 
load('RecurrenceStats_FilterPts')
dt = 0.01;

for iPrep = 1:nPreps
    % get coalescence time
    TimeToCoalesce = [RcrStats(iPrep,:).firstW] - pars.stimstart;  % seconds to reach coalscence compared to min possible time we can measure

    for iStim = 1:nStims
        ixP = Prep(iPrep).PIdx(iStim).ts; ixP = ixP(1:round(TimeToCoalesce(iStim)/dt));  % set of admissible points within projection
        targetPCs = Prep(iPrep).PCprojU(ixP,1:Prep(iPrep).PCA.nPCs); 
        for j = iStim+1:nStims  
            % set of admissible time-points in source program
            ixP = Prep(iPrep).PIdx(j).ts;  ixP = ixP(1:round(TimeToCoalesce(j)/dt));
            srcPCs = Prep(iPrep).PCprojU(ixP,1:Prep(iPrep).PCA.nPCs); % get subset of concatenated PCs that correspond to this progam
            IC(iPrep).sponD(j,iStim) = sqrt(sum((srcPCs(1,:) - targetPCs(1,:)).^2));  % distance at onset of stim
             
            for iN = 1:numel(ixP)
                pt0 = srcPCs(iN,:);  % current x_i
                D = sqrt(sum(bsxfun(@minus,targetPCs,pt0).^2,2)); % Euclidean distance to all points within time radius
                Dmax(iStim,j,iN) = max(D); % max distance from this point to the other trajectory
            end
            IC(iPrep).stimDmax(j,iStim) = max(squeeze(Dmax(iStim,j,:))); % max of all max distnaces
        end
        % keyboard
    end
end

% correlate spontaneous distance and max evoked distance
cmap = brewermap(10,'Paired'); 

figure; hold on
allDmax = [];
for iPrep = 1:nPreps
    Dspon = IC(iPrep).sponD(IC(iPrep).sponD > 0);
    Dmax = IC(iPrep).stimDmax(IC(iPrep).stimDmax > 0);
    plot(Dspon,Dmax,'o','MarkerFaceColor',cmap(iPrep,:),'MarkerEdgeColor','none','MarkerSize',5);
    xlabel('Distance at stimulation')
    ylabel('Maximum divergence')
    allDmax = [allDmax; Dmax];
end

% correlate spontaneous and Haussdorf on attractor
figure; hold on
allSpon = [];
allHdist = [];
for iPrep = 1:nPreps
    Dspon = IC(iPrep).sponD(IC(iPrep).sponD > 0);
    data = Prep(iPrep).HaussD(Prep(iPrep).HaussD > 0);
    plot(Dspon,data,'o','MarkerFaceColor',cmap(iPrep,:),'MarkerEdgeColor','none','MarkerSize',5);
    xlabel('Distance at stimulation')
    ylabel('Distance on manifold')
    
    allSpon = [allSpon; Dspon];
    allHdist = [allHdist; data];
end

[r,p] = corr(allSpon,allHdist);

% parallel plot
allDistances = [allSpon allDmax allHdist];
figure
LinkedUnivariateScatterPlots(gca,1:3,[allSpon allDmax allHdist],[0.7 0.7 0.7]);

save InitialConditions IC pars allDistances



