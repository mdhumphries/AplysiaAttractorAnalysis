% Extended Figure summarising evidence that P10 does not contain
% motorneurons

clear all; close all

if ispc
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
    rawpath = 'C:\Users\mqbssmhg.DS\Dropbox\SpikeData\Angela Bruno Aplysia Dynamics Data\P10\';
else
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';
end

run figure_properties

ixExamples = 5; % prep #5, stims 2 and 3


%% examples of P10 firing recordings...

% load raw unit data
load([rawpath 'Jun0414B_raw']);
ts = (0:numel(da01)-1) .* (125 / numel(da01)); 
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 7 3]); hold on
plot(ts,da01,'Linewidth',errorlinewidth,'Color',[0 0 0])
axis tight
axis off
line([30 30],[-1300 1500],'Linewidth',axlinewidth,'Color',stimcolor)
line([110,120],[-1000,-1000],'Linewidth',plotlinewidth,'Color',[0 0 0],'Clipping','off')
text(115,-1100,'10 s','FontSize',fontsize)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_P10MN_da01_Jun0414B


load([rawpath 'Jun0914_Raw.mat']);
ts = (0:numel(da02)-1) .* (125 / numel(da02)); 
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 7 3]); hold on
plot(ts,da02,'Linewidth',errorlinewidth,'Color',[0 0 0])
axis tight
axis off
line([30 30],[-1300 1500],'Linewidth',axlinewidth,'Color',stimcolor)
line([110,120],[-1200,-1200],'Linewidth',plotlinewidth,'Color',[0 0 0],'Clipping','off')
text(115,-1300,'10 s','FontSize',fontsize)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_P10MN_da02_Jun0914


%% panel: summary of all programs
load([filepath 'P10 analysis/P10_MN_Test'])

seqall = [ones(3,1); ones(3,1)+1; ones(3,1)+2];

h12 = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 7 7]); hold on
h13 = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 7 7]); hold on

da01files = [P10data(seqall==1).iR];
da02files = [P10data(seqall==2).iR];
da03files = [P10data(seqall==3).iR];
i1 = find(seqall==1);
i2 = find(seqall==2);
i3 = find(seqall==3);

for iC = 1:3
    match2 = find(da01files(iC) == da02files);
    match3 = find(da01files(iC) == da03files);

    % da01 result
    maxP_1 = P10data(i1(iC)).maxPstim;
    maxP_1spon = P10data(i1(iC)).maxPspon;

    % da02 result
    maxP_2 = P10data(i2(match2)).maxPstim;
    maxP_2spon = P10data(i2(match2)).maxPspon;
    
    % matching da03 result
    maxP_3 = P10data(i3(match3)).maxPstim;
    maxP_3spon = P10data(i3(match3)).maxPspon;
   
    Nspon12 = (P10data(i2(match3)).nspikesSpon + P10data(i1(iC)).nspikesSpon)/2;
    Nstim12 = (P10data(i2(match3)).nspikes + P10data(i1(iC)).nspikes)/2;

    Nspon13 = (P10data(i3(match3)).nspikesSpon + P10data(i1(iC)).nspikesSpon)/2;
    Nstim13 = (P10data(i3(match3)).nspikes + P10data(i1(iC)).nspikes)/2;
    
    % 1 vs 2
    figure(h12)
    scatter(maxP_1(maxP_1 > 0 & maxP_2 > 0),maxP_2(maxP_1 > 0 & maxP_2 > 0),...
        100*Nstim12(maxP_1 > 0 & maxP_2 > 0)/max(Nstim12),'MarkerEdgeColor',prep_Cmap(da01files(iC),:)); hold on

    % 1 vs 3
    figure(h13)
    scatter(maxP_1(maxP_1 > 0 & maxP_3 > 0),maxP_3(maxP_1 > 0 & maxP_3 > 0),...
        100*Nstim13(maxP_1 > 0 & maxP_3 > 0)/max(Nstim13),'MarkerEdgeColor',prep_Cmap(da01files(iC),:)); hold on
end
figure(h12)
axis square
axis([0 1 0 1])
set(gca,'XTick',0:0.2:1,'YTick',0:0.2:1)
xlabel('P(spike_{P10} | spike) : stim. 1','FontSize',fontsize)
ylabel('P(spike_{P10} | spike) : stim. 2','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_P10MN_1vs2

figure(h13)
axis square
axis([0 1 0 1])
set(gca,'XTick',0:0.2:1,'YTick',0:0.2:1)
xlabel('P(spike_{P10} | spike) : stim. 1','FontSize',fontsize)
ylabel('P(spike_{P10} | spike) : stim. 3','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_P10MN_1vs3

