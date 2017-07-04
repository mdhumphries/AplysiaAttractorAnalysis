% script to use P10 firing to ID putative motorneurons
% ID any with really strong match (look at match distributions etc) 
% ALSO: must be consistent between same prep programs!!

clear all; close all

datapath = '../../Data/P10/';
spikepath = '../../Data/Spikes/';

addpath ../../Functions/

MNpars.stimstart = 30;
MNpars.endwin = 123;  % allow for overrun of recording
MNpars.strtMatch = 33;  % to avoid stimulation artifacts, start after stim ends...

MNpars.GaussQt = 0.001;  % (s); bin size and grid size for kernel density estimate 
MNpars.GaussSD = 0.002;

MNpars.Tvar = 0.002;  % size of +/- window around target spike
MNpars.d = [0:0.001:0.025];  % delay range to test in ms, 
MNpars.dbins = round(MNpars.d ./ MNpars.GaussQt); %converted into bins

%% find P10s in filesets
files = dir(datapath); files(1:2) = []; 

ixda01 = find(arrayfun(@(x) any(findstr(x.name,'da01')),files)); % index of all da02 files
ixda02 = find(arrayfun(@(x) any(findstr(x.name,'da02')),files)); % index of all da02 files
ixda03 = find(arrayfun(@(x) any(findstr(x.name,'da03')),files)); % index of all da03 files

ixall = [ixda01; ixda02; ixda03];
seqall = [ones(numel(ixda01),1); ones(numel(ixda02),1)+1; ones(numel(ixda03),1)+2];

nfiles = numel(ixall); % numel(ixda01) + numel(ixda02) + numel(ixda03);

load ../da03_DataProperties_FunctionAndWindowSize FileTable
da03Files = FileTable;

load ../da02_DataProperties_FunctionAndWindowSize FileTable
da02Files = FileTable;

load ../da01_DataProperties_FunctionAndWindowSize FileTable
da01Files = FileTable;

%% get P10 and match spikes
for iF = 1:nfiles
    % filename of P10 file
    P10data(iF).name =  files(ixall(iF)).name;
    % get P10 spikes
    load([datapath P10data(iF).name]);
    P10data(iF).spks = spks(:,2);

    % keyboard
    
    % matching index of recorded data file
    if seqall(iF) == 1
        FTable = da01Files;
    elseif seqall(iF) == 2
        FTable = da02Files;
    else
        FTable = da03Files;
    end
    % get data-file spikes
    P10data(iF).iR = find(cellfun(@(x) strncmp(P10data(iF).name,x,7),FTable));  % first 7 characters; index of matching dataset
    load([spikepath FTable{P10data(iF).iR}]);  % "spks" is now dataset...
    P10data(iF).IDs = unique(spks(:,1));
    
    % start/end
    P10data(iF).ixstrt = find(P10data(iF).spks >= MNpars.stimstart,1);
    P10data(iF).ixend = find(P10data(iF).spks <= MNpars.endwin,1,'last');
    T = MNpars.endwin - MNpars.stimstart;
        
%     % convolve P10 with Gaussian 
%     sig = GaussSD / GaussQt; % SD in time-steps
%     x = [-5*sig:1:5*sig]';  % x-axis values of the discretely sampled Gaussian, out to 5xSD
%     h = (1/(sqrt(2*pi*sig^2)))*exp(-((x.^2*(1/(2*sig^2))))); % y-axis values of the Gaussian
%     h = h ./sum(h); % make sure kernel has unit area, then can use for rate functions
%     shiftbase = floor(length(x)/2); 
%     
%     T = [stimstart endwin];
%     P10data(iF).bins = T(1):GaussQt:T(2);  % all time-points to evaluate
%     
%     % compute spike-train convolution functions
%     [P10data(iF).spkfcn,idxs] = convolve_spiketrains([ones(numel(spksP10),1) spksP10],h,shiftbase,1,P10data(iF).bins,GaussQt,T,'Gaussian');
%     figure
%     plot(P10data(iF).bins,P10data(iF).spkfcn)
    
%     S = zeros(numel(P10data(iF).IDs),numel(d)); 
    P10data(iF).Pstim = zeros(numel(P10data(iF).IDs),numel(MNpars.d));
    P10data(iF).Pspon = zeros(numel(P10data(iF).IDs),numel(MNpars.d));
    for iN = 1:numel(P10data(iF).IDs)
        ts = spks(spks(:,1) == P10data(iF).IDs(iN),2);
        ts(ts <= MNpars.strtMatch) = []; % omit stimulation period
%         ixBins = arrayfun(@(x) find(P10data(iF).bins >= x,1,'first'),ts); % quantise to underlying bins
        %ixMatch = 1:ceil(numel(ts)/2);
        ixMatch = 1:numel(ts); 
        P10data(iF).nspikes(iN) = numel(ixMatch);
        if ~isempty(ts)  % i.e. some spikes available to match!
            for iD = 1:numel(MNpars.d)
                % sum of all P10 convolution...
                % S(iN,iD) = sum(P10data(iF).spkfcn(ixBins(ixMatch) + dbins(iD)));
                % NO: should be sum of spike-train convolution: weighted by each spike (and then normalised) 

                % OR: just p(spikeP10|spike_x) within +/- Tranmission variation ms
                Ps = arrayfun(@(x) any(P10data(iF).spks >= x+MNpars.d(iD)-MNpars.Tvar & P10data(iF).spks <= x+MNpars.d(iD)+MNpars.Tvar),ts(ixMatch));
                P10data(iF).Pstim(iN,iD) = sum(Ps) / numel(Ps);
            end
        end
        
        % and do spon period
        ts = spks(spks(:,1) == P10data(iF).IDs(iN),2);
        ts(ts >= MNpars.stimstart) = []; % keep only spontaneous period
        %ixMatch = 1:ceil(numel(ts)/2);
        ixMatch = 1:numel(ts); 
        P10data(iF).nspikesSpon(iN) = numel(ixMatch);
        if ~isempty(ts)  % i.e. some spikes available to match!
            for iD = 1:numel(MNpars.d)
                % OR: just p(spikeP10|spike_x) within +/- Tranmission variation ms
                Ps = arrayfun(@(x) any(P10data(iF).spks >= x+MNpars.d(iD)-MNpars.Tvar & P10data(iF).spks <= x+MNpars.d(iD)+MNpars.Tvar),ts(ixMatch));
                P10data(iF).Pspon(iN,iD) = sum(Ps) / numel(Ps);
            end
        end
       
    end
    
%     figure
%     colormap(cbrewer('seq','YlOrRd',10))
%     imagesc(d,IDs,S);
 
% uncomment to get example for one population
%     figure
%     colormap(cbrewer('seq','YlOrRd',10))
%     %subplot(121),
%     imagesc(MNpars.d,P10data(iF).IDs,P10data(iF).Pstim);
%     xlabel('Spike-to-P10 Delay (s)');
%     ylabel('Neuron');    
%     % title(['Recording ' num2str(P10data(iF).iR) ': stim'])
%     colorbar
%     % exportPPTfig(gcf,'Example_P10_P_spike')
%     
%     subplot(122),
%     imagesc(MNpars.d,P10data(iF).IDs,P10data(iF).Pspon);
%     xlabel('Spike-to-P10 Delay (s)');
%     ylabel('Neuron');    
%     title(['Recording ' num2str(P10data(iF).iR) ': spon'])
%     colorbar    
    %keyboard
   
    
    % look at distribution of max(P) for each neuron (i.e. over all delays)
    P10data(iF).maxPstim = max(P10data(iF).Pstim,[],2);
    maxP_delay = bsxfun(@eq,P10data(iF).Pstim,P10data(iF).maxPstim); % matrix of matches to max P value: often more than one per neuron
    P10data(iF).min_Delay_stim = arrayfun(@(x) MNpars.d(min(find(maxP_delay(x,:)==1))),P10data(iF).IDs);  % get just earliest one per neuron
    
    P10data(iF).maxPspon = max(P10data(iF).Pspon,[],2);
    maxP_delay = bsxfun(@eq,P10data(iF).Pspon,P10data(iF).maxPspon); % matrix of matches to max P value: often more than one per neuron
    P10data(iF).min_Delay_spon = arrayfun(@(x) MNpars.d(min(find(maxP_delay(x,:)==1))),P10data(iF).IDs);  % get just earliest one per neuron
 
    figure
    subplot(221),hist(P10data(iF).maxPstim,30); xlabel('P(spike_{P10} | spike-stim)'); ylabel('No. neurons')
    subplot(222), plot(P10data(iF).maxPstim,P10data(iF).min_Delay_stim,'k.','MarkerSize',10); 
    xlabel('P(spike_{P10} | spike-stim)'); ylabel('Min. delay (s)')
    subplot(223),hist(P10data(iF).maxPspon,30); xlabel('P(spike_{P10} | spike-spon)'); ylabel('No. neurons')
    subplot(224), plot(P10data(iF).maxPspon,P10data(iF).min_Delay_spon,'k.','MarkerSize',10); 
    xlabel('P(spike_{P10} | spike-spon)'); ylabel('Min. delay (s)')

    figure
    subplot(121),
    plot(P10data(iF).maxPspon,P10data(iF).maxPstim,'r+')
    xlabel('P(spike_{P10} | spike-spon)'); ylabel('P(spike_{P10} | spike-stim)');
    subplot(122),
    plot(P10data(iF).min_Delay_spon,P10data(iF).min_Delay_stim,'b+')
    xlabel('Min. delay(s) - spon'); ylabel('Min. delay(s) - stim')
   
    % keyboard
    
end

save P10_MN_Test P10data MNpars

%% compare sensitised and rest matches
da02files = [P10data(seqall==2).iR];
da03files = [P10data(seqall==3).iR];
i2 = find(seqall==2);
i3 = find(seqall==3);

cmap = cbrewer('qual','Paired',nfiles); 
figure
for iC = 1:numel(ixda02)
    match3 = find(da02files(iC) == da03files);
    
    % da02 result
    maxP_2 = P10data(i2(iC)).maxPstim;
    maxP_2spon = P10data(i2(iC)).maxPspon;
    
    % matching da03 result
    maxP_3 = P10data(i3(match3)).maxPstim;
    maxP_3spon = P10data(i3(match3)).maxPspon;
    
    Nspon = (P10data(i3(match3)).nspikesSpon + P10data(i2(iC)).nspikesSpon)/2;
    Nstim = (P10data(i3(match3)).nspikes + P10data(i2(iC)).nspikes)/2;
    
    
    % plot rest vs sensitised for each preparation
    
   
    scatter(maxP_2(maxP_2 > 0 & maxP_3 > 0),maxP_3(maxP_2 > 0 & maxP_3 > 0),...
        100*Nstim(maxP_2 > 0 & maxP_3 > 0)/max(Nstim),'MarkerEdgeColor',cmap(iC,:)); hold on
    % title('Bubble size = #spikes')
    
%     figure
%     subplot(121),scatter(maxP_2,maxP_3,100*Nstim/max(Nstim))
%     xlabel('da02: max P stim'); ylabel('da03: max P stim'); 
%     title('Bubble size = #spikes')
%     subplot(122),scatter(maxP_2spon,maxP_3spon,100*Nspon/max(Nspon))
%     xlabel('da02: max P spon'); ylabel('da03: max P spon');    
%     title('Bubble size = #spikes')

end

xlabel('P(spike_{P10}|spike) : Stim.2'); ylabel('P(spike_{P10}|spike) : Stim.3');
axis([0 1 0 1])
exportPPTfig(gcf,'Pop_Not_P10MNs')


