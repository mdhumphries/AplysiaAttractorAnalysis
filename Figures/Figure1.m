%% script to assemble sub-panels for Figure 1: raster and rates
clear all; close all

filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
spikepath = 'C:\Users\mqbssmhg.DS\Dropbox\paper2\datasets\deconcatenated_spks\experimental\';

run figure_properties.m

flag = 'B3'; % clusters only

offset = 2; % additional colormap entries to make sure the lightest ones are omitted

load([filepath 'da02_DataProperties_FunctionAndWindowSize'],'FileTable')
load([filepath '/Ensembles/da02_Ensemble_Types.mat'],'AcorrTypes')
load([filepath '/Ensembles/da02_Analyses_Neurons_and_Groups.mat'],'GroupList')
load([filepath '/Ensembles/da02_StaticEnsembles.mat'],'StaticEnsemblesAll')

%% check spread of types in each recording
% find out which programs had nice mix of ensemble types...
     
distTypes = zeros(numel(FileTable),4);
for iF = 1:numel(FileTable)
    % look up in GroupList
    ixRec = GroupList(:,1) == iF;
    recTypes = AcorrTypes(ixRec);
        
    % get spread of types...
    distTypes(iF,:) = hist(recTypes,1:4);
end

% examine distribution of types, and pick one with mixture
ixExample = 3;

%% raster plot of all ensembles in one example recording: ordered by TYPE first
% so:
% pick one, and 
% fine-tune by hand to get phase-sweep of oscillators?

% load spike-trains
load([spikepath FileTable{ixExample}]);

% use indices into types to sort group IDs
recTypes = AcorrTypes(GroupList(:,1) == ixExample);
spkIDs = StaticEnsemblesAll(ixExample).grps;  % group IDs of all neurons
grpctr = 1; Cctr=1;
cmapGrp = []; 

for iT = [1,4,3,2]  % order from bottom to top of raster 
    % find groups that are this type
    ixR = find(recTypes == iT);
    
    % keyboard
    % renumber 
    for iG = ixR'  % only works if is a row vector...
        spkIDs(StaticEnsemblesAll(ixExample).grps == iG) = grpctr;
        grpctr = grpctr+1;
    end
    
    % built type colormap
    Cmap2 = brewermap(2+offset,Typeclrs{Cctr});
    thisCmap = repmat(Cmap2(end-1:end,:),ceil(numel(ixR)/2),1);  % alternating shade colormap
    thisCmap = thisCmap(1:numel(ixR),:);
    % thisCmap = brewermap(numel(ixR)+offset,clrs{Cctr}); % add offset (+1) to stay away from lightest colours...
    % keyboard
    % cmapGrp = [cmapGrp; thisCmap(offset+1:end,:)];
    cmapGrp = [cmapGrp; thisCmap];
    Cctr = Cctr +1;
    
end

spkIDs = [unique(spks(:,1)) spkIDs];  % 2 column vector: neuron ID and grou

% plot clusters
[H,GS,I,O] = plot_clusters(spks,spkIDs,numel(recTypes),[0 max(spks(:,2))+1],'3',[],cmapGrp,M);

set(H,'Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 8 8]); 
title('');
axis([20 120 -5 size(spkIDs,1)])
line([110,120],[-5,-5],'Linewidth',plotlinewidth,'Color',[0 0 0])
text(115,-9,'10 s','FontSize',fontsize)

% axis([0.5 6.5 Tstart Tstart + Tlength])
axis off

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc Fig1c_raster
exportfig(gcf,'Fig1c_raster','Color',color,'Format',format,'Resolution',dpi)

%% panel d: network spectra

load([filepath 'da01_DataProperties_FunctionAndWindowSize'],'Data')
da01_Data = Data;
load([filepath 'da02_DataProperties_FunctionAndWindowSize'],'Data')
da02_Data = Data;
load([filepath 'da03_DataProperties_FunctionAndWindowSize'],'Data')
da03_Data = Data;


allspectra = zeros(numel(Data)*3,numel(Data(1).fspectra));
for iF = 1:numel(Data)
    allspectra(iF,:) = da01_Data(iF).Spectra;
    allspectra(iF+numel(Data),:) = da02_Data(iF).Spectra;
    allspectra(iF+numel(Data)*2,:) = da03_Data(iF).Spectra;
end

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); 
semilogx(Data(1).fspectra,allspectra,'Color',[0.7 0.7 0.7],'Linewidth',errorlinewidth); hold on
semilogx(Data(1).fspectra,median(allspectra,1),'Color',[0 0 0],'Linewidth',plotlinewidth)
%axis([0 125 0 6])
xlabel('Frequency (Hz)','FontSize',fontsize)
ylabel('Power','FontSize',fontsize)
axis([0.01 1 0 25])
set(gca,'XTick',[10^-2 10^-1 1],'XTickLabel',[0.01 0.1 1])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
print -depsc Fig1d_spectra


%% panel e: network rates

load([filepath 'Rates_and_Synchrony\da01_CorrelationEvolution.mat']);
da01_Netrate = reshape([CxyDataSet(:).networkrate],122,10);
load([filepath 'Rates_and_Synchrony\da02_CorrelationEvolution.mat']);
da02_Netrate = reshape([CxyDataSet(:).networkrate],122,10);

load([filepath 'Rates_and_Synchrony\da03_CorrelationEvolution.mat']);
da03_Netrate = reshape([CxyDataSet(:).networkrate],122,10);

allNetrate = [da01_Netrate da02_Netrate da03_Netrate];
t = 0:1:121;

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on
ptch = fill([30 32.5 32.5 30],[6 6 0 0],stimcolor);
set(ptch,'EdgeColor',stimcolor)
plot(t,allNetrate,'Color',[0.7 0.7 0.7],'Linewidth',errorlinewidth)
plot(t,median(allNetrate,2),'Color',[0 0 0],'Linewidth',plotlinewidth)
axis([0 125 0 6])
xlabel('Time (s)','FontSize',fontsize)
ylabel('Network rate (spikes/neuron/s)','FontSize',fontsize)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
print -depsc Fig1e_netrates



%% schematic loop

x = 0:0.01:50;
y1 = sin(x);
y2 = sin(x+2);

y1 = [linspace(1.5,y1(1),500).^2 y1];
y2 = [linspace(0,y2(1),500) y2];

y1 = y1+ randn(1,numel(y1))*0.05;
y2 = y2+ randn(1,numel(y2))*0.05;

figure
subplot(211),plot(y1); hold on
plot(y2);
subplot(212),
plot(y1,y2)








