% Figure 3 extension: showing spontaneus jumps are perturbations of attractor,
% or not....
%
% (1) number of perturbations across all programs
% (2) number returned to attractor
% (3) if not, them the proportion of cycle until end of recording (i.e.
% could there be no recurrence because end was too close?)
% (4) If returned, then same or different attractor?

clear all; close all

if ismac
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';
else
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
end

run figure_properties

fnames = {'da01','da02','da03'};
ixT = 3;  % last threshold
M = 4; % need bigger markers here 
omitcolor = [0.7 0.7 0.7];

load([filepath 'Classify attractors/JumpAnalysis2'])  % get data on jumps
load([filepath 'Classify attractors/RecurrenceStats_FilterPts'])  % get data on final perturbations...
%% panel a: example density time-series with a jump
load([filepath 'da02_DataProperties_FunctionAndWindowSize'],'Data')
load([filepath '/Classify attractors/da02_RecurrencePointsData_FilterPts'])

ixExample = 9;

tStim = Data(ixExample).stimbins(1);
tEnd = Data(ixExample).stimbins(end) - tStim;

% tFinal = Data(ixExample).stimbins(end) - pars.maxWindow - tStim;  % stopping for checking recurrence before end
tFinal = Data(ixExample).stimbins(RcrPtData(ixExample).winmids(end)) - tStim;  % stopping for checking recurrence before end

% calculate percentage of recurrent points per window (any points)
PrctWper = 100*(1-RcrPtData(ixExample).Orbit(statpars.ixT).NoRecurrence.nWinPoints/RcrPtData(ixExample).winDT);

tsDensity = 100*RcrPtData(ixExample).Orbit(statpars.ixT).AllRecurrentPoint.pRecur;

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2]); hold on
% line([tFinal tFinal],[0 100],'Color',[0 0 0])

% % omit peri-stimulus
% ptch = fill([0 pars.T0-tStim pars.T0-tStim 0],[100 100 0 0],omitcolor);
% set(ptch,'EdgeColor',omitcolor)
% 
% % omit before end
% ptch = fill([tFinal tEnd tEnd tFinal],[100 100 0 0],omitcolor);
% set(ptch,'EdgeColor',omitcolor)
% text(80,105,'omitted','Fontsize',fontsize-1,'Color',omitcolor)

% stimulus
ptch = fill([0 2.5 2.5 0],[100 100 0 0],stimcolor);
set(ptch,'EdgeColor',stimcolor)
text(-10,115,'stimulus','Fontsize',fontsize-1,'Color',stimcolor)

% label jump
% tJumpSt = Jump(ixExample,2).Transition(1).tsJumpStart - tStim;
% tJumpEnd = Jump(ixExample,2).Transition(1).tsJumpEnd - tStim;
tJumpSt = Data(ixExample).stimbins(RcrPtData(ixExample).winmids(Jump(ixExample,2).Transition(1).ixWinJumpStart)) - tStim;
tJumpEnd = Data(ixExample).stimbins(RcrPtData(ixExample).winmids(Jump(ixExample,2).Transition(1).ixWinJumpEnd)) - tStim;

% bar over the top
% line([tJumpSt tJumpEnd],[110 110],'Linewidth',plotlinewidth,'Color',[0 0 0],'Clipping','off');
% line([0 90],[50 50],'Color',[0.7 0.7 0.7])  % threshold line

% shade jump
ptch = fill([tJumpSt tJumpEnd tJumpEnd tJumpSt],[100 100 0 0],[0.9 0.9 0.9]);
set(ptch,'EdgeColor',[0.9 0.9 0.9])

% label coalscence
%h = annotation('textarrow',[0 0 ],[0 1]);
%set(h,'parent',gca);
% set(h,'position',[RcrStats(ixExample,2).firstW-tStim 50 0 20]);
line([RcrStats(ixExample,2).firstW-tStim RcrStats(ixExample,2).firstW-tStim],[0 100],'Color',[0.6 0.6 0.6])

% data
hR = plot(Data(ixExample).stimbins(RcrPtData(ixExample).winmids) - tStim,tsDensity,'o');
set(hR,'Markersize',2,'Color',prep_Cmap(ixExample,:));
% set(gca,'YTick',[0 5 10 15])

% label checked window
tsChecked = Data(ixExample).stimbins(RcrPtData(ixExample).winmids(Jump(ixExample,2).Transition(1).ixWinBefore)) - tStim;
plot(tsChecked,tsDensity(Jump(ixExample,2).Transition(1).ixWinBefore),'o','Markersize',2,'MarkerFaceColor',[0 0 0])


axis([0 90 0 100])
set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap

% hl = line([0 2.5],[15.5 15.5],'Linewidth',1.5,'Color',stimcolor,'Clipping','off');
xlabel('Time (s)','FontSize',fontsize)
ylabel('Density (%)','FontSize',fontsize)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc Fig3_ExampleDensityTimeSeries

%% panel a: matching example recurrence delay time-series
tsDelay = RcrPtData(ixExample).Orbit(statpars.ixT).AllRecurrentPoint.mDelay;

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2]); hold on
% stimulus
ptch = fill([0 2.5 2.5 0],[50 50 0 0],stimcolor);
set(ptch,'EdgeColor',stimcolor)
% text(-10,55,'stimulus','Fontsize',fontsize-1,'Color',stimcolor)

% label coalscence
line([RcrStats(ixExample,2).firstW-tStim RcrStats(ixExample,2).firstW-tStim],[0 50],'Color',[0.6 0.6 0.6])

% shade jump
tJumpSt = Data(ixExample).stimbins(RcrPtData(ixExample).winmids(Jump(ixExample,2).Transition(1).ixWinJumpStart)) - tStim;
tJumpEnd = Data(ixExample).stimbins(RcrPtData(ixExample).winmids(Jump(ixExample,2).Transition(1).ixWinJumpEnd)) - tStim;
ptch = fill([tJumpSt tJumpEnd tJumpEnd tJumpSt],[50 50 0 0],[0.9 0.9 0.9]);
set(ptch,'EdgeColor',[0.9 0.9 0.9])

% data
hR = plot(Data(ixExample).stimbins(RcrPtData(ixExample).winmids) - tStim,tsDelay,'o');
set(hR,'Markersize',2,'Color',prep_Cmap(ixExample,:));
% set(gca,'YTick',[0 5 10 15])

% label checked window
tsChecked = Data(ixExample).stimbins(RcrPtData(ixExample).winmids(Jump(ixExample,2).Transition(1).ixWinBefore)) - tStim;
plot(tsChecked,tsDelay(Jump(ixExample,2).Transition(1).ixWinBefore),'o','Markersize',2,'MarkerFaceColor',[0 0 0])

axis([0 90 0 50])
set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap

% hl = line([0 2.5],[15.5 15.5],'Linewidth',1.5,'Color',stimcolor,'Clipping','off');
xlabel('Time (s)','FontSize',fontsize)
ylabel('Recurrence time (s)','FontSize',fontsize)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc Fig3_ExampleDelayTimeSeries



%% panel g: number of perturbations per program



Xvalues = [sum(~Summary.blnReturned); sum(Summary.blnReturned); Summary.allJumps];

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2]); hold on; 
barh(1:3,Xvalues(1:3),'Facecolor',[0 0 0])
axis([0 20 0.5 3.5])
xlabel('Occurrence','FontSize',fontsize);

text(Xvalues(1)-1,1,num2str(Xvalues(1)),'FontSize',fontsize-1,'Color',[1 1 1]);
text(Xvalues(2)-2,2,num2str(Xvalues(2)),'FontSize',fontsize-1,'Color',[1 1 1]);
text(Xvalues(3)-2,3,num2str(Xvalues(3)),'FontSize',fontsize-1,'Color',[1 1 1]);

% ylabel('Event','FontSize',fontsize)
% set(gca,'YTick',1:4,'YTickLabel',{'Not returned','Returned','Potential','All'});
set(gca,'YTick',1:3,'YTickLabel',{'Not returned','Returned','All'});
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc Fig3g_Types

% exportPPTfig(gcf,'TypeOfPerturbs',[10 15 5 7])

%% plot same/different
pAfter = Summary.pAfterEnd(~isnan(Summary.pAfterEnd));
[srtD,Idens] = sort(pAfter,'descend');  % sort into order of prop below threshold
%figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 2.5 4]); hold on; 
%line([50 50],[0.5 numel(Idens)+0.5],'Linewidth',axlinewidth,'Color',[0.7 0.7 0.7],'Linestyle','--');
%plot(100*Summary.pAfter(Idens),1:numel(Idens),'ko','Markersize',2); hold on
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2.5]); hold on; 
plot(1:numel(Idens),100*pAfter(Idens),'ko','Markersize',2); hold on
line([0.5 numel(Idens)+0.5],[50 50],'Linewidth',axlinewidth,'Color',[0.7 0.7 0.7],'Linestyle','--');

% for iN = 1:numel(Idens)
%     strTick{iN} = [num2str(Summary.JumpList(Idens(iN),1)) ',' num2str(Summary.JumpList(Idens(iN),2)) ',' num2str(Summary.JumpList(Idens(iN),3))];
% end
% % set(gca,'YTick',1:numel(Idens),'YTickLabel',strTick)        

set(gca,'Fontsize',fontsize)
axis([0.5 numel(Idens)+0.5 0 100])
xlabel('Divergence No.','FontSize',fontsize); 
ylabel('Same manifold (%)','FontSize',fontsize); 

% ylabel('Jump (Prep, Program, No.)','FontSize',fontsize)
% ylabel('Returned perturbation','FontSize',fontsize)
% set(gca,'XTick',[0 10 20 30],'YTick',[0 0.005 0.01]);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc Fig3h_Same

% exportPPTfig(gcf,'Jumps',[10 15 5 8])

%% plot difference from expected proportion
% expectedP = max(allSummary(:,10:11),[],2);  % proportion of rotation period
% expectedP(expectedP > 100) = 100;
% dExpect = allSummary(:,9)-expectedP; % difference in expected proportion
% [srtExpect,Iexp] = sort(dExpect);  % sort into order of difference in expected proportion
% 
% figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 8]); hold on; 
% plot(dExpect(Iexp),1:numel(Iexp),'ko'); hold on
% 
% for iN = 1:numel(Iexp)
%     strTick{iN} = [num2str(allSummary(Iexp(iN),4)) ',' num2str(allSummary(Iexp(iN),5)) ',' num2str(allSummary(Iexp(iN),6))];
% end
% set(gca,'YTick',1:numel(Iexp),'YTickLabel',strTick)        
% set(gca,'Fontsize',fontsize)
% % axis([0 30 0 0.01])
% xlabel('Deviation from expected points (%)','FontSize',fontsize); 
% ylabel('Jump (Prep, Program)','FontSize',fontsize)
% % set(gca,'XTick',[0 10 20 30],'YTick',[0 0.005 0.01]);
% set(gca,'FontName','Helvetica','FontSize',fontsize);
% set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
% 


% print -depsc ExtFig_RQA_Example_histogram

%% panel b: scatters of program quantities

