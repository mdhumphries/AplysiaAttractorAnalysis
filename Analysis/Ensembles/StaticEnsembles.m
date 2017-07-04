%% script to do static clustering over whole recording
% ADD STRUCT TO COMMENTS
% StaticEnsemblesAll
%
%
% Mark Humphries

clear all; close all

addpath('../../Functions/')


spikepath = '../../Data/Spikes/';


fname_Pre = []; % 'Control'; % [];
stimset = 'da01';  % 'da01': first; 'da02': second; and 'da03' : third (control)

stimstart = 30;
stimoffset = 32.5;

% load dataset properties
fname = ['../' fname_Pre stimset '_DataProperties_FunctionAndWindowSize'];
load(fname)

%% parameters
nfiles = numel(FileTable);

clstrstart = stimstart; % + 6; % 7 seconds after stimulation: mean time to recover phase-locking... (cf ComparePLandR)

% % load windowed ensembles
% fname = [fname_Pre stimset '_EnsemblesSlidingWindow'];
% load(fname,'WinEnsembles')

tic
%% do each recording
for iD = 1:nfiles
    iD
    spkdatafile = FileTable{iD};

    % load spike data
    load([spikepath spkdatafile]);

    % remove spontaneous period from this analysis
    % spks(spks(:,2) < stimstart,:) = [];

    % recording times
    startts = 0; % floor(DataTable(iD,2));
    endts = floor(DataTable(iD,3));

    % all spike-train IDs
    StaticEnsemblesAll(iD).allIDs = unique(spks(:,1));

    % consensus parameters - derived from stimulation onset only
    Qt = Data(iD).GaussQt;  
    OPTS.BLpars = Data(iD).GaussSD;


    %% static clustering: over whole program
    MPspks = spks; 
    % MPspks(MPspks(:,2) < stimstart,:) = []; 
    MPspks(MPspks(:,2) < clstrstart,:) = []; 
    
    StaticEnsemblesAll(iD).IDs = unique(MPspks(:,1));
    
    % use full function here to regenerate SDFs - allows for flexibility in
    % determing start point etc
    [GMAX,GCON,Sxy,spkfcn] = consensus_cluster_spike_data_binless(MPspks,StaticEnsemblesAll(iD).IDs,[clstrstart DataTable(iD,3)],Qt,OPTS);
    
%     B = SD * sqrt(12); 
%     opts.overlap = 0;
%     [GMAXb,GCONb,Sxyb,spkcntb] = consensus_cluster_spike_data_bins(MPspks,StaticEnsemblesAll(iD).IDs,B,[clstrstart DataTable(iD,3)],opts);
    
%     plot_clusters(spks,GCON.grps,GCON.ngrps,[0 DataTable(iD,3)],'3B')
%     plot_clusters(MPspks,GCONb.grps,GCONb.ngrps,[clstrstart DataTable(iD,3)],'3B')
    
    % keyboard
     
    % StaticEnsemblesAll(iD).Cxy = Sxy{1}; % correlation of original spike-train functions
    StaticEnsemblesAll(iD).Cxy = corrcoef(spkfcn{1});
    StaticEnsemblesAll(iD).grps = GCON.grps(:,2);
    StaticEnsemblesAll(iD).grpsizes = GCON.grpsizes;
    StaticEnsemblesAll(iD).Q = GCON.Q;
    StaticEnsemblesAll(iD).Gmax = GMAX;
    StaticEnsemblesAll(iD).grpIDs = unique(StaticEnsemblesAll(iD).grps); 
    StaticEnsemblesAll(iD).ngrps = numel(StaticEnsemblesAll(iD).grpIDs);

    % properties of static ensembles
    for iS = 1:StaticEnsemblesAll(iD).ngrps

        theseIDs = StaticEnsemblesAll(iD).IDs(StaticEnsemblesAll(iD).grps == iS); 
        % compute phase-locking of whole MP
        EnsmblSpks = [];
        for iN = 1:numel(theseIDs)
            EnsmblSpks = [EnsmblSpks;  find(MPspks(:,1) == theseIDs(iN))];
        end

       %% PL to main oscillation? phase of spikes....
       ixs = arrayfun(@(x) find(Data(iD).stimbins <= x,1,'last'),MPspks(EnsmblSpks,2)); % find all phases at spike-times
       StaticEnsemblesAll(iD).PC1_PL_MP(iS) = mean(PCAdata(iD).phase(ixs)); % take mean of these
       
       %% compute windowed firing rates
       
     
       
    end
    
    
%     %% track network: relationship to static clustering
%     % save with this struct as match here is property of new static
%     % clusterin, not of windowed clustering....
%     for iW = 1:WinEnsembles(iD).nwins
%         [commonIDs,I1,I2] = intersect(StaticEnsemblesAll(iD).IDs,WinEnsembles(iD).cellIDs{iW}); % find common nodes between static and window
%         StaticEnsemblesAll(iD).MIstatic(iW) = MIpartitions(StaticEnsemblesAll(iD).grps(I1),WinEnsembles(iD).wgrpscon{iW}(I2));
%         StaticEnsemblesAll(iD).VIstatic(iW) = VIpartitions(StaticEnsemblesAll(iD).grps(I1),WinEnsembles(iD).wgrpscon{iW}(I2)) ./ log(numel(commonIDs));
%     %     MIstatic2(iW) = MIpartitions(Sgrps,WinGCON(iW).grps(:,2));
%     %     VIstatic2(iW) = VIpartitions(Sgrps,WinGCON(iW).grps(:,2)) ./ log(numel(WinGCON(iW).grps(:,2)));
%     end
% 
% % for iD = 1:nfiles
%     % omit windows pre-stim onset that finish after stim onset....
%     winends = WinEnsembles(iD).winstrt+WinEnsembles(iD).winsize;
%     cens = ~(WinEnsembles(iD).winstrt < stimstart & winends > stimstart);
%     
% figure
%     subplot(211),plot(WinEnsembles(iD).winstrt(cens),WinEnsembles(iD).MIstatic(cens)); hold on
%     line([stimstart stimstart],[0.5 0.8],'Color',[0.5 0.5 0.5])
%     %plot(MIstatic2,'r')
%     ylabel('NMI window vs static')
%     subplot(212),plot(WinEnsembles(iD).winstrt(cens),WinEnsembles(iD).VIstatic(cens)); hold on
%     line([stimstart stimstart],[0.2 0.5],'Color',[0.5 0.5 0.5])
%     %plot(VIstatic2,'r')
%     ylabel('VI window vs static')
end



%% SAVE all data
save([fname_Pre stimset '_StaticEnsembles'],'StaticEnsemblesAll')

clear('StaticEnsemblesAll')