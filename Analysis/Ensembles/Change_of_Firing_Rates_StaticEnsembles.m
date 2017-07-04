%% script to assess firing rate changes at level of ensemble...

clear all; close all

addpath('../../Functions/')

spikepath = '../../Data/Spikes/';

fname_Pre = []; % 'Control'; % [];
stimset = {'da01','da02','da03'};  % 'da01': first; 'da02': second; and 'da03' : third (control)

% load dataset properties
load(['../' fname_Pre stimset{1} '_DataProperties_FunctionAndWindowSize'],'FileTable');

nPreps = numel(FileTable);

% stimulation parameters
pars.stimstart = 30;

% window sizes
pars.win = 10:5:40;  % window sizes in seconds
pars.wstep = 5;
pars.corrType = 'Pearson';

%% get rates
for iPrep = 1:nPreps
   for iStim = 1:numel(stimset)
        load([stimset{iStim} '_StaticEnsembles.mat'],'StaticEnsemblesAll')

        % load spike data
        spkdatafile = FileTable{iPrep};
        load([spikepath spkdatafile]);

        MPspks = spks; 
        MPspks(MPspks(:,2) < pars.stimstart,:) = [];    
        EnsembleRates(iPrep,iStim).Ntotal = numel(unique(MPspks(:,1)));  % total number of neurons
        
        for iS = 1:StaticEnsemblesAll(iPrep).ngrps
            theseIDs = StaticEnsemblesAll(iPrep).IDs(StaticEnsemblesAll(iPrep).grps == iS);
            
            EnsembleRates(iPrep,iStim).Ensemble(iS).N = numel(theseIDs); % number of neurons in this ensemble
            
            % get all spikes
            EnsmblSpks = [];
            for iN = 1:numel(theseIDs)
                EnsmblSpks = [EnsmblSpks;  MPspks(MPspks(:,1) == theseIDs(iN),2)];
            end
            tEnd = max(EnsmblSpks);
            % keyboard
            for iW = 1:numel(pars.win)
                % for a range of window sizes, compute spike count per window
                EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).winstart = pars.stimstart:pars.wstep:tEnd - pars.win(iW);
                EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).winend = EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).winstart + pars.win(iW);           
                EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).winmids = EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).winstart + pars.win(iW)/2;
                
                for iCnt = 1:numel(EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).winstart)
                    EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).counts(iCnt) = sum(EnsmblSpks > EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).winstart(iCnt)...
                                                                                & EnsmblSpks <= EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).winend(iCnt));
                end
                  
                
                % and regress with time
                EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).R = corr(EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).winmids',EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).counts',...
                                                                'Type',pars.corrType);
            end
         end
   end
end

save Rates_StaticEnsembles EnsembleRates pars