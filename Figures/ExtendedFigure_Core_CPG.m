clear all; close all

if ispc
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
else
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';
end

run figure_properties
fnames = {'da01','da02','da03'};

% M = 4; % need bigger markers here 

load([filepath '/Classify attractors/ParticipationChanges_SamePCs']) % participation scores
load([filepath '/Ensembles/EnsembleTypeConsistency']) % neurons labelled by ensemble type

npreps = size(Coeffs,2);
nprogs = 3;

dist_panel = [4 4];

scatter_panel = [8 8];

offset = 2;
for i = 1:4
    Cmap2 = brewermap(2+offset,Typeclrs{i});
    typeCmap(i,:) = Cmap2(end,:);
end

%% distribution of persistence, across each program
YdataCell = cell(4,1);
Ydata = 100.*(reshape([TypeConsistency.AllPropPersistent],4,npreps))';
for i =1:4 YdataCell{i} = Ydata(:,i); end

strXlabel = {'Tonic','Osc.','Burst','Pause'};
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 dist_panel]); hold on
UnpairedUnivariateScatterPlots(gca,YdataCell,'MarkerEdgeColor',typeCmap,'strXlabel',strXlabel)
set(gca,'YLim',[0 100])
ylabel('Neurons consistently labelled (%)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc CoreCPG_consistency_spread

%% plot all participation stuff for persistent oscillators
% also plot all with per-prep colour map
allMaxP = []; allMaxDP = [];
hall = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 dist_panel]); hold on

for iPrep = 1:npreps
%     figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 scatter_panel]); hold on
    ixPlot = TypeConsistency(iPrep).PersistentTypeNeuronIDs{2}; % plot oscillators
%  
%     plot(Coeffs(iPrep).HighMaxP(ixPlot),Coeffs(iPrep).maxDMaxP(ixPlot),'ko')
%     axis([0 1 0 1])
%     xlabel('Highest Proportion of Maximum','FontSize',fontsize)
%     ylabel('Maximum change in proportion','FontSize',fontsize)
%     title(['Prep ' num2str(iPrep) ': consistent oscillators'])
%     set(gca,'FontName','Helvetica','FontSize',fontsize);
%     set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
    
    % all
    allMaxP = [allMaxP; Coeffs(iPrep).HighMaxP(ixPlot)];
    allMaxDP = [allMaxDP;  Coeffs(iPrep).maxDMaxP(ixPlot)];
    figure(hall)
    plot(100*Coeffs(iPrep).HighMaxP(ixPlot),100*Coeffs(iPrep).maxDMaxP(ixPlot),'o','MarkerEdgeColor',prep_Cmap(iPrep,:),'MarkerSize',M)
end

figure(hall)
axis([0 100 0 100])
xlabel('Maximum participation (%)','FontSize',fontsize)
ylabel('Maximum change in participation (%)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc CoreCPG_participation_Scatter_PrepColors

% plot all....
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 dist_panel]); hold on
plot(100*allMaxP,100*allMaxDP,'ko','MarkerSize',M)
axis([0 100 0 100])
xlabel('Maximum participation (%)','FontSize',fontsize)
ylabel('Maximum change in participation (%)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc CoreCPG_participation_Scatter


%% summarise per program
threshold = 90;  % find 75% percentile in data
for iPrep = 1:npreps
    ixPlot = TypeConsistency(iPrep).PersistentTypeNeuronIDs{2}; % get oscillators
    P = Coeffs(iPrep).HighMaxP(ixPlot);
    DP = Coeffs(iPrep).maxDMaxP(ixPlot);
    
    % get upper X-tile of Participation
    UpperQ = prctile(P,threshold);
    ixUpper = find(P > UpperQ);
    UpperP = P(ixUpper);
    UpperDP = DP(ixUpper);
    
    % summarise
    MedianUpperP(iPrep) = median(UpperP);
    MedianUpperDP(iPrep) = median(UpperDP);
end

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 dist_panel]); hold on
for i=1:npreps plot(100*MedianUpperP(i),100*MedianUpperDP(i),'o','MarkerEdgeColor',prep_Cmap(i,:)); end
axis([0 100 0 50])
xlabel('Median of maximum participation (%)','FontSize',fontsize)
ylabel('Median of maximum change (%)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc CoreCPG_participation_UpperQ



