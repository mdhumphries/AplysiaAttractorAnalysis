%% load data-set and get all stats on it
% Choose convolution and window size:
% (1) Convolution size: ISI stats should dictate ~ 500 ms
% (2) Window size linked to primary oscillation: access using PCA analysis
%
% OUTPUT:
% for all recordings:
%   FileTable: cell array of all filenames
%   DataTable: summary table of each recording's properties: [#trains T_start T_end median rate]
%   struct DataStats
%       DataStats.Hks : matrix; rejection of null hypothesis {0,1} for all pairs of firing rate distributions (KS test)
%       DataStats.Pks : corresponding P-value
%       DataStats.KSSTAT : corresponding KS-statistic
%       DataStats.BH_H : array of still-significant rejections at alpha=0.05 after Benjamini-Hochberg test
%
% per recording R,
%       struct Data:
% Data(iR).rates : firing rate of every train, from stimulus onset
% Data(iR).medrate : median rate of recording
% Data(iR).fpeak : peak frequency of mean spike-train spectrum (Chronux) 
% Data(iR).ISIprct : inter-spike intervals at the [2.5 25 50 75 97.5] % of the ISI distribution
% Data(iR).GaussQt: quantisation step for Gaussian convolution
% Data(iR).GaussSD: standard deviation of convolved Gaussian
% Data(iR).bins : time-stamps of convolved functions for all recording (also used for PCA projections)
% Data(iR).stimbins : time-stamps of convolved functions from stimulation onset only (also used for PCA projections)     
% Data(iR).ixts: time indexes of sub-set used for fitting de-trend regression, and for peak detection
%
%       struct SDF:
% SDF(iR).spkfcn: the complete set of Gaussian-convolved spike-density functions for each recording
%
%       struct PCAdata
% PCAdata(iR).allcoeffs : principal components (eigenvectors) - whole recording
% PCAdata(iR).allPCproj : projection of time-series onto each PC axis (weighted sum) - whole recording
% PCAdata(iR).alleigvalues : eigenvalues of each PC - whole recording
% PCAdata(iR).coeffs : principal components (eigenvectors) - from stim
% PCAdata(iR).PCproj : projection of time-series onto each PC axis (weighted sum) - from stim
% PCAdata(iR).eigvalues : eigenvalues of each PC - from stim
% PCAdata(iR).prctVar : cumulative distribution of variance accounted for (%)
% PCAdata(iR).nPCs : number of PCs accounting for 95% variance (closest number)
% PCAdata(iR).PCdetrend : detrended clock PC (usually PC2) - removed linear regression
% PCAdata(iR).filterPC: low-pass filter of detrended PC (1/4 of period)  
% PCAdata(iR).f : frequency axis for power spectrum in [0 10] Hz
% PCAdata(iR).S : power spectrum
% PCAdata(iR).fpeak : peak frequency of PC1 oscillation  (cycle period =1/fpeak)
% PCAdata(iR).phase : instantaneous phase of PC1 oscillation in [0, 2*pi]
% PCAdata(iR).Pksraw : detected peaks in detrended PC1 (index into time-stamps)
% PCAdata(iR).PksrawFlip : detected troughs in detrended PC1 (index into time-stamps) - flip around y-axis to find first cycle 
% PCAdata(iR).Pksfilter : detected peaks in filtered detrended PC1 (index into time-stamps)
% PCAdata(iR).PksrawFlip : detected troughs in filtered detrended PC1 (index into time-stamps) - flip around y-axis to find first cycle 
% PCAdata(iR).tpeaks : time of peak in each cycle (seconds; starting from
% earliest peak across both normal and inverted filtered PC1 signals)
% PCAdata(iR).tepochs : time-points of each epoch (seconds; [0 stim-start Peak#1 Peak#2... Peak#N]
%
% NEEDS: Chronux Toolbox to compute spike-train spectra: http://chronux.org/
%
% Mark Humphries 2/2/2015



clear all
close all

addpath ../Functions/

% path to Chronux toolbox: edit here %%%%%%%%%%%%%%%%%%%
if ismac
    gpath = genpath('/Users/mqbssmhg/Dropbox/My Toolboxes/Other tools/chronux/');
else
    gpath = genpath('C:\Users\mqbssmhg.DS\Dropbox\My Toolboxes\Other tools\chronux\');  % PC
end

addpath(gpath)

% expects that data-files always contain variable "spks"
% 3-stim data
spikepath = '../Data/Spikes/'; % path
fname_Pre = '';  % prefix
stimset = 'da03';  % 'da01': first; 'da02': second; and 'da03' : third


% data-set parameters
stimstart = 30; % 30 s into recording
stimoff = 32.5; % 2.5 s stimulation


%% Analysis parameters
% Gaussian convolution
Qt = 0.01;

% PCA
VarExplained = 0.95;  % 95%

% MT spectrum properties - spikes
params.Fs = 1000; % sampling time-step is 1 ms
params.trialave = 1; % average over all spike-trains
params.fpass = [0 1];
params.err = [2 0.05];

% MT spectrum properties - PC1
paramsPC.Fs = 1/Qt; % sampling time-step is 10 ms for convolved functions
paramsPC.trialave = 0; % average over all spike-trains
paramsPC.fpass = [0 10];



%% loop and analyse data
dirP = dir(spikepath); dirP = dirP(3:end);   % discard top two entries: always cd/ commands

nfiles = numel(dirP);
Dctr = 0;
DataTable = []; allrates = []; allisis = [];
for loop = 1:nfiles
    SpkData = struct([]);  % need spike-times in struct to pass to Chronux functions
    if any(strfind(dirP(loop).name,'.mat')) && any(strfind(dirP(loop).name,stimset)) 
        Dctr = Dctr + 1; % number of recordings in data-set
        load([spikepath '/' dirP(loop).name]);
        if ~exist('spks','var') error(['Data-file ' dirP(loop).name ' does not contain spks variable']); end

        % T_start_recording = min(spks(:,2));
        T_start_recording = stimstart;  % start from 30 s in
        T_end_recording = max(spks(:,2));
        
        FileTable{Dctr} = dirP(loop).name;  % keep list of filenames
        
        % remove spontaneous period from this analysis
        spks(spks(:,2) < stimstart,:) = [];
        n_trains = numel(unique(spks(:,1)));  % count number in the retained period

        Data(Dctr).rates = zeros(1,n_trains);
        isis = [];
        for iN = 1:n_trains
            SpkData(iN).ts = spks(spks(:,1)==iN,2);
            Data(Dctr).rates(iN) = sum(spks(:,1)==iN) ./ (T_end_recording - T_start_recording); 
            isis = [isis; diff(SpkData(iN).ts)];
        end
        allrates = [allrates; Data(Dctr).rates'];
        allisis = [allisis; isis];
        figure(1); hold on
        ecdf(Data(Dctr).rates)
        % title(['Dataset ' num2str(Dctr)])
        Data(Dctr).medRate = median(Data(Dctr).rates);
        
        %% get spectrum
        [Data(Dctr).Spectra,Data(Dctr).fspectra,rate,Data(Dctr).Serr] = mtspectrumpt(SpkData,params,1);
        
        figure
        shadedErrorBar(Data(Dctr).fspectra,Data(Dctr).Spectra,flipud(Data(Dctr).Serr),'r')

        
        % get peak frequency
        ixHigh = find(Data(Dctr).Spectra > mean(Data(Dctr).Spectra)+std(Data(Dctr).Spectra)*2);
        Shigh = Data(Dctr).Spectra(ixHigh) / sum(Data(Dctr).Spectra(ixHigh));
        Data(Dctr).fpeak = Data(Dctr).fspectra(ixHigh) * Shigh;  % (dot-product) weighted-mean value of f
        
        %% get all ISI stats
        Data(Dctr).ISIprct = prctile(allisis,[2.5 25 50 75 97.5]);
        
        % quick reference table
        DataTable = [DataTable; n_trains T_start_recording T_end_recording Data(Dctr).medRate];
    end
end
medianISI = arrayfun(@(x) x.ISIprct(3),Data);

medianAll = median(allisis);  % median ISI for whole data-set

%%  check for significant differences between firing rate distributions - quantify variation in equivalent programs
for r = 1:Dctr
    for c = r+1:Dctr
        [Hks(c,r),Pks(c,r),KSSTAT(c,r)] = kstest2(Data(r).rates,Data(c).rates);
    end
end

DataStats.Hks = [Hks zeros(Dctr,1)]; DataStats.Pks = [Pks zeros(Dctr,1)]; DataStats.KSSTAT = [KSSTAT zeros(Dctr,1)];

% check for multiple comparisons
DataStats.Z = squareform(DataStats.Pks);
[DataStats.BH_H,DataStats.BH_T] = benjaminihochberg(DataStats.Z,0.05);
DataStats.Z(DataStats.BH_H==1)

%% convolve and do PCA
% SD = 0.5;
for iR = 1:numel(FileTable)
    load([spikepath  FileTable{iR}]);

    IDs = unique(spks(:,1));
    
    Data(iR).GaussSD = Data(iR).ISIprct(3);  % median ISI for this recording
    Data(iR).GaussQt = Qt;
    
    % convolve with Gaussian
    sig = Data(iR).GaussSD / Data(iR).GaussQt; % SD in time-steps
    x = [-5*sig:1:5*sig]';  % x-axis values of the discretely sampled Gaussian, out to 5xSD
    h = (1/(sqrt(2*pi*sig^2)))*exp(-((x.^2*(1/(2*sig^2))))); % y-axis values of the Gaussian
    h = h ./sum(h); % make sure kernel has unit area, then can use for rate functions
    shiftbase = floor(numel(x)/2); 
    
    T = [0 DataTable(iR,3)];
    Data(iR).bins = T(1):Data(iR).GaussQt:T(2);  % all time-points to evaluate
    
    % compute spike-train convolution functions
    [SDF(iR).spkfcn,idxs] = convolve_spiketrains(spks,h,shiftbase,IDs,Data(iR).bins,Data(iR).GaussQt,T,'Gaussian');
    
    [PCAdata(iR).allcoeffs, PCAdata(iR).allPCproj, PCAdata(iR).alleigvalues] = princomp(SDF(iR).spkfcn);

    % remove spontaneous part of this analysis
    ixStim = find(Data(iR).bins >= stimstart);
    Data(iR).stimbins = Data(iR).bins(ixStim);  % from stimulation onset only
    % returns: eigenvectors (coefficients); projection onto each
    % eigenvector (weighted sum, minus mean); eigenvalues
    [PCAdata(iR).coeffs, PCAdata(iR).PCproj, PCAdata(iR).eigvalues] = princomp(SDF(iR).spkfcn(ixStim,:));
    
     % compute variance explained
     PCAdata(iR).prctVar = cumsum(PCAdata(iR).eigvalues) / sum(PCAdata(iR).eigvalues);
     PCAdata(iR).nPCs = sum(PCAdata(iR).prctVar <= VarExplained);

     % also get non-normlaused PCs as linear sums over spike functions
%      PCAdata(iR).PCprojU = zeros(numel(ixStim),PCAdata(iR).nPCs);  % unnormalised PCs
%      for iP = 1:PCAdata(iR).nPCs 
%         PCAdata(iR).PCprojU(:,iP) = sum(bsxfun(@times,SDF(iR).spkfcn(ixStim,:),PCAdata(iR).coeffs(:,iP)'),2);  % raw
%         % PCAdata(iR).PCproj(:,iP) = sum(bsxfun(@times,bsxfun(@minus,SDF(iR).spkfcn,mean(SDF(iR).spkfcn)),PCAdata(iR).coeffs(:,iP)'),2); % normalised to mean rate
%      end
    
     % keyboard
end

%% get PC1 and PC2 and get freq; 

for iR = 1:numel(FileTable)
    figure(101); clf; plot(PCAdata(iR).PCproj(:,1)); hold on
    plot(PCAdata(iR).PCproj(:,2),'r');
     title(['Recording ' num2str(iR)])

    % detrend with linear regression - using time after stimulation
    % offset...

    Data(iR).ixts = find(Data(iR).stimbins > stimoff+2.5);  % portion of time used to fit detrend and do peak detection
    [B,BINT,R,RINT,STATS] = regress(PCAdata(iR).PCproj(Data(iR).ixts,1),[ones(numel(Data(iR).ixts),1) Data(iR).bins(Data(iR).ixts)']);
    [B2,BINT2,R2,RINT2,STATS2] = regress(PCAdata(iR).PCproj(Data(iR).ixts,2),[ones(numel(Data(iR).ixts),1) Data(iR).bins(Data(iR).ixts)']);

    % detrend whole thing
    PC1detrend = PCAdata(iR).PCproj(:,1) - (B(1) + B(2) * Data(iR).stimbins)';
    PC2detrend = PCAdata(iR).PCproj(:,2) - (B2(1) + B2(2) * Data(iR).stimbins)';
           
    % power spectrum
    % mPC1 = mean(PCAdata(iR).PCproj(Data(iR).ixts,1));
    % [PCAdata(iR).S,PCAdata(iR).f] = mtspectrumc(PCAdata(iR).PCproj(:,1) - mPC1,paramsPC);
    [S,f] = mtspectrumc(PC1detrend,paramsPC);  % of detrended data
    [S2,f2] = mtspectrumc(PC2detrend,paramsPC);  % of detrended data

     figure(110), plot(f,S); hold on
     plot(f2,S2,'r');
     
    % peak frequency of PCs
    ixHigh1 = find(S > mean(S)+std(S)*2);  % region 2 SDs above mean
    Shigh1 = S(ixHigh1) / sum(S(ixHigh1));  % p(power) in that region
    fpeak1 = f(ixHigh1) * Shigh1;   % (dot-product) weighted-mean value of f
    
    ixHigh2 = find(S2 > mean(S2)+std(S2)*2); % region 2 SDs above mean
    Shigh2 = S2(ixHigh2) / sum(S2(ixHigh2));
    fpeak2 = f2(ixHigh2) * Shigh2;  % (dot-product) weighted-mean value of f
    
    % CHOOSE: if weighted frequency is higher, then likely to have greater
    % power at f>0...
    if fpeak1 > fpeak2
        PCAdata(iR).PCclock = 1;  
        PCAdata(iR).fpeak = fpeak1; 
        PCAdata(iR).PCdetrend = PC1detrend;
        PCAdata(iR).S = S;
        PCAdata(iR).f = f;
    else
        PCAdata(iR).PCclock = 2;  
        PCAdata(iR).fpeak = fpeak2;
        PCAdata(iR).PCdetrend = PC2detrend;
        PCAdata(iR).S = S2;
        PCAdata(iR).f = f2;
    end
    
    % instantaneous phase of PC
    H = hilbert(PCAdata(iR).PCdetrend - mean(PCAdata(iR).PCdetrend));
    PCAdata(iR).phase = angle(H) + pi; % phase on [0,2*pi]
    
    % find peaks in PC: use to define epochs
    winsize = round(1./PCAdata(iR).fpeak / Qt / 4);
    fcutoff = 1./(winsize * Qt);
    if ~rem(winsize,2)
        winsize = winsize+1;
    end
    PCAdata(iR).filterPC = movingAverage(PCAdata(iR).PCdetrend(Data(iR).ixts),winsize);

    pad = zeros(round(10/Qt),1);  % pad end by 10 seconds to find end peak
    PCAdata(iR).Pksraw = findpeaksAMPD2([PCAdata(iR).PCdetrend(Data(iR).ixts); pad]);
    PCAdata(iR).Pksfilter = findpeaksAMPD2([PCAdata(iR).filterPC; pad]);
    pksMat = findpeaks([PCAdata(iR).filterPC; pad]);
    pksMat = pksMat.loc;
    
    % negative peaks...
    PCAdata(iR).PksrawFlip = findpeaksAMPD2([PCAdata(iR).PCdetrend(Data(iR).ixts)*-1; pad]);
    PCAdata(iR).PksfilterFlip = findpeaksAMPD2([PCAdata(iR).filterPC*-1; pad]);
    
    pksMatUD = findpeaks([PCAdata(iR).filterPC*-1; pad]);
    pksMatUD = pksMatUD.loc;
    
    figure(102); clf; 
    plot(Data(iR).stimbins,PCAdata(iR).PCdetrend,'k'); hold on
    plot(Data(iR).stimbins(Data(iR).ixts),PCAdata(iR).filterPC,'r')
    plot(Data(iR).stimbins(Data(iR).ixts(PCAdata(iR).Pksraw)),PCAdata(iR).PCdetrend(Data(iR).ixts(PCAdata(iR).Pksraw)),'ks')
    plot(Data(iR).stimbins(Data(iR).ixts(PCAdata(iR).Pksfilter)),PCAdata(iR).filterPC(PCAdata(iR).Pksfilter),'rs')
    plot(Data(iR).stimbins(Data(iR).ixts(pksMat)),PCAdata(iR).filterPC(pksMat),'rd')
    title(['Recording ' num2str(iR)])
    
    figure(103); clf
    plot(Data(iR).stimbins,-1*PCAdata(iR).PCdetrend,'k'); hold on
    plot(Data(iR).stimbins(Data(iR).ixts),-1*PCAdata(iR).filterPC,'r')
    plot(Data(iR).stimbins(Data(iR).ixts(PCAdata(iR).PksrawFlip)),-1*PCAdata(iR).PCdetrend(Data(iR).ixts(PCAdata(iR).PksrawFlip)),'ks')
    plot(Data(iR).stimbins(Data(iR).ixts(PCAdata(iR).PksfilterFlip)),-1*PCAdata(iR).filterPC(PCAdata(iR).PksfilterFlip),'rs')
    plot(Data(iR).stimbins(Data(iR).ixts(pksMatUD)),-1*PCAdata(iR).filterPC(pksMatUD),'rd')
    title(['Recording ' num2str(iR)])
    
   % find which to keep
   if PCAdata(iR).Pksfilter(1) < PCAdata(iR).PksfilterFlip(1) 
       PCAdata(iR).tpeaks = Data(iR).stimbins(Data(iR).ixts(PCAdata(iR).Pksfilter));
   else
       PCAdata(iR).tpeaks = Data(iR).stimbins(Data(iR).ixts(PCAdata(iR).PksfilterFlip));
   end
    
   PCAdata(iR).tepochs = [0 stimstart PCAdata(iR).tpeaks];  % divisions into epochs
end
% 
% %% periods;
% 
% dataperiod = 1./[PCAdata(:).fpeak]
% 
% % change in individual recordings....
% figure
% for iR = 1:numel(FileTable)
%     dfs = diff(PCAdata(iR).tpeaks);
%     % plot(PCAdata(iR).tpeaks(1:end-1) - PCAdata(iR).tpeaks(1),dfs,'.-'); hold on
%     
%     plot(1:numel(dfs),dfs,'.-'); hold on
% end
% % xlabel('Time since first peak (s)')
% xlabel('Cycle number')
% ylabel('Elapsed time since last peak (s)')

%% is #PCs to 95% variance simply a function of number of neurons or overall rate? NO

Npcs = [PCAdata(:).nPCs];

figure
subplot(211),plot(DataTable(:,1),Npcs,'k.')
xlabel('# neurons in recording')
ylabel('#PCs to reach 95%')
subplot(212),plot(DataTable(:,4),Npcs,'k.')
xlabel('Median spike rate of recording')
ylabel('#PCs to reach 95%')



%%
fname = [fname_Pre stimset '_DataProperties_FunctionAndWindowSize'];
save(fname,'Data','DataTable','FileTable','PCAdata','DataStats','SDF')

clear('Data','DataTable','FileTable','PCAdata','DataStats','SDF')