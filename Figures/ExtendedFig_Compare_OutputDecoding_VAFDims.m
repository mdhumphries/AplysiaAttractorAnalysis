%% script to assemble sub-panels for Figure 4: decoding
clear all; close all

if ispc
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
else
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';
end

run figure_properties

M = 4; % need bigger markers here

fname = {'da01','da02','da03'}; 

VAFsuffix = '_VAF_90'; % [];

% load overall P10 stats
load([filepath 'P10 analysis/P10_statespace_Stats'])
Stats80 = SummaryStats;
load([filepath 'P10 analysis/P10_statespace_Stats' VAFsuffix])
Stats90 = SummaryStats;

TrainColor = [0.4 0.4 0.8];
FitColor = [0.8 0.4 0.4];

% FitColor = [153,50,204] / 256;

MAEcolor = [0 0 0];
Rcolor = [0.4 0.8 0.4];

ErrorColor = [0.6 0.6 0.6];

stimcolor = [0.7 0.7 0.7];
tStim = 30;

%% panel: compare P10 dimension changes here...
D80 = [Stats80(:).D]; 
D90 = [Stats90(:).D];

strXlabel = {'80%','90%'};
lcolor = [0.7 0.7 0.7];

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 2 3]); hold on;
LinkedUnivariateScatterPlots(gca,1:2,[D80' D90'],lcolor,'strXlabel',strXlabel,'MarkerSize',M);
ylabel('Embedding dimensions','FontSize',fontsize); 
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc SIFig_GLM_VAF_compareD


%% panel a: example fit and forecast
iProg = 3; iPrep = 3;  % stim 3, prep 3

% load P10 data for that prep...
load([filepath 'P10 analysis/' fname{iProg} '_StateSpace_P10_GLMmodel' VAFsuffix]);

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 2.5]);
plot(P10data(iPrep).bins - tStim,P10data(iPrep).spkfcn,'k','Linewidth',plotlinewidth); hold on
plot(P10data(iPrep).bins(VizP10(iPrep,iProg).forecastNt(1:1000)) - tStim,...
    VizP10(iPrep,iProg).tsP10Model(numel(VizP10(iPrep,iProg).Nt)+1:numel(VizP10(iPrep,iProg).Nt)+1000),...
    'Color',FitColor,'Linewidth',errorlinewidth)
plot(P10data(iPrep).bins(VizP10(iPrep,iProg).Nt) - tStim,VizP10(iPrep,iProg).tsP10Model(1:numel(VizP10(iPrep,iProg).Nt)),...
    'Color',TrainColor,'Linewidth',errorlinewidth)
axis([0 95 0 350])

ptch = fill([0 2.5 2.5 0],[350 350 0 0],stimcolor);
set(ptch,'EdgeColor',stimcolor)
set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap

xlabel('Time (s)','Fontsize',fontsize)
ylabel('P10 firing (spikes/s)','Fontsize',fontsize)
line([P10data(iPrep).bins(VizP10(iPrep,iProg).Nt(1)),P10data(iPrep).bins(VizP10(iPrep,iProg).Nt(end))] - tStim,[350,350],...
    'Color',TrainColor,'Linewidth',plotlinewidth)
line([P10data(iPrep).bins(VizP10(iPrep,iProg).Nt(end)+1),P10data(iPrep).bins(VizP10(iPrep,iProg).Nt(end)+1000)] - tStim,[350,350],...
    'Color',FitColor,'Linewidth',plotlinewidth)
text(P10data(iPrep).bins(VizP10(iPrep,iProg).Nt(1)) - tStim,375,'Training fit','Color',TrainColor,'Fontsize',fontsize-1)
text(P10data(iPrep).bins(VizP10(iPrep,iProg).Nt(end)+1) - tStim,375,'Forecast fit','Color',FitColor,'Fontsize',fontsize-1)
text(75,175,'P10 data','Color',[0 0 0],'Fontsize',fontsize-1)

% exportPPTfig(gcf,['Prep' num2str(iP) '_Program' num2str(iF) ... 
%        '_example_forecast_R' num2str(Stats90(iP,iF).medianR) '_MAE' num2str(Stats90(iP,iF).medianMAE) '.png'],[10 15 6 4])

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc SIFig_GLM_VAF_exampleforecast

% exportPPTfig(gcf,'ExampleForecast',[10 15 6 4])

% error metrics
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 2.5]);
[ax,h1,h2] = plotyy(VizP10(iPrep,iProg).tsForecastStart - tStim,VizP10(iPrep,iProg).tsMAEforecast,...
    VizP10(iPrep,iProg).tsForecastStart - tStim,VizP10(iPrep,iProg).tsRforecast);

% cpos = get(ax(1),'Position');
set(ax,'Position',[0.17 0.28 0.65 0.61]);
set(ax,'XLim',[VizP10(iPrep,iProg).tsForecastStart(1) VizP10(iPrep,iProg).tsForecastStart(end)] - tStim)
set(ax(1),'YLim',[0 16],'YColor',[0 0 0],'YTick',[0 5 10 15])
set(ax(2),'YLim',[0 1],'YColor',[0 0 0],'YTick',[0 0.25 0.5 0.75 1])

set(h2,'Color',Rcolor,'Linewidth',plotlinewidth)
set(h1,'Color',MAEcolor,'Linewidth',plotlinewidth)

xlabel(ax(1),'Forecast start (s)','FontSize',fontsize)
ylabel(ax(1),'MAE (spikes/s)','FontSize',fontsize)
% ylabel(ax(2),'R','FontSize',fontsize,'Rotation',270)

text(65,17,'R','Color',Rcolor,'Fontsize',fontsize-1)
text(65,3,'MAE','Color',MAEcolor,'Fontsize',fontsize-1)

set(ax,'FontName','Helvetica','FontSize',fontsize);
set(ax,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc SIFig_GLM_VAF_example_metrics

% exportPPTfig(gcf,'ExampleMetrics',[10 15 6 4])


%% panel b: summary scatter

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on;
for iPrep = 1:3
    for iProg = 1:3
        % MAE error bar
        semMAE = Stats90(iPrep,iProg).stdMAE / sqrt(numel(VizP10(iPrep,iProg).tsForecastStart));
        line([Stats90(iPrep,iProg).meanMAE-2*semMAE,Stats90(iPrep,iProg).meanMAE+2*semMAE],[Stats90(iPrep,iProg).meanR2 Stats90(iPrep,iProg).meanR2],'Linewidth',errorlinewidth,'Color',ErrorColor)
        % R error bar
        semR2 = Stats90(iPrep,iProg).stdR2 / sqrt(numel(VizP10(iPrep,iProg).tsForecastStart));
        line([Stats90(iPrep,iProg).meanMAE,Stats90(iPrep,iProg).meanMAE],[Stats90(iPrep,iProg).meanR2-2*semR2 Stats90(iPrep,iProg).meanR2+2*semR2],'Linewidth',errorlinewidth,'Color',ErrorColor)
%         % MAE std bar
%         line([Stats90(iP,iF).meanMAE-Stats90(iP,iF).stdMAE,Stats90(iP,iF).meanMAE+Stats90(iP,iF).stdMAE],[Stats90(iP,iF).meanR2 Stats90(iP,iF).meanR2],'Linewidth',errorlinewidth,'Color',ErrorColor)
%         % R std bar
%         line([Stats90(iP,iF).meanMAE,Stats90(iP,iF).meanMAE],[Stats90(iP,iF).meanR2-Stats90(iP,iF).stdR2 Stats90(iP,iF).meanR2+Stats90(iP,iF).stdR2],'Linewidth',errorlinewidth,'Color',ErrorColor)
        
        % data-point
        plot([Stats90(iPrep,iProg).meanMAE],[Stats90(iPrep,iProg).meanR2],'o',...
            'MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k'); hold on

    end
end

% figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3 3]); hold on;
% % plot([Stats90(:).meanMAE],[Stats90(:).meanR].^2,'r.'); hold on
% plot([Stats90(:).medianMAE],[Stats90(:).medianR].^2,'k.','MarkerSize',M);
axis([0 40 0 1])
xlabel('MAE (spikes/s)','FontSize',fontsize); 
ylabel('R^2','FontSize',fontsize)

% set(gca,'XTick',[0 10 20 30])
set(gca,'XTick',[0 10 20 30 40])

%exportPPTfig(gcf,['Summary_R2_vs_MAE'],[10 15 5 5])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');


print -depsc SIFig_GLM_scatter

% exportPPTfig(gcf,'P10GLM_summary',[10 15 8 8])


%% compare 80 and 90% VAF
MAEmean80 = [Stats80(:).meanMAE]; 
MAEmean90 = [Stats90(:).meanMAE];
R2mean90 = [Stats90(:).meanR2];
R2mean80 = [Stats80(:).meanR2];

strXlabel = {'80%','90%'};
lcolor = [0.7 0.7 0.7];

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 2 3]); hold on;
LinkedUnivariateScatterPlots(gca,1:2,[MAEmean80' MAEmean90'],lcolor,'strXlabel',strXlabel,'MarkerSize',M);
ylabel('Mean MAE (spikes/s)','FontSize',fontsize); 
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
print -depsc SIFig_CompareMAE


figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 2 3]); hold on;
LinkedUnivariateScatterPlots(gca,1:2,[R2mean80' R2mean90'],lcolor,'strXlabel',strXlabel,'MarkerSize',M);
ylabel('Mean R^2','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
print -depsc SIFig_CompareR2

% %% panel b: example forecasts over rest of program
% 
% % best R^2
% iProg = 3; iPrep = 3;  % P10 dataset stim 3, prep 3 = Prep 5, stim 3 in full data-set
% load([filepath 'P10 analysis/' fname{iProg} '_StateSpace_P10_GLMmodel']);
% 
% figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3.5 2.1]);
% plot(P10data(iPrep).bins - tStim,P10data(iPrep).spkfcn,'k','Linewidth',plotlinewidth); hold on
% plot(P10data(iPrep).bins(VizP10(iPrep,iProg).forecastNt) - tStim,...
%     VizP10(iPrep,iProg).tsP10Model(numel(VizP10(iPrep,iProg).Nt)+1:end),...
%     'Color',FitColor,'Linewidth',errorlinewidth)
% plot(P10data(iPrep).bins(VizP10(iPrep,iProg).Nt) - tStim,VizP10(iPrep,iProg).tsP10Model(1:numel(VizP10(iPrep,iProg).Nt)),...
%     'Color',TrainColor,'Linewidth',errorlinewidth)
% axis([0 95 0 350])
% ptch = fill([0 2.5 2.5 0],[350 350 0 0],stimcolor);
% set(ptch,'EdgeColor',stimcolor)
% set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap
% 
% set(gca,'YTick',[0 100 200 300]);
% xlabel('Time (s)','Fontsize',fontsize-1)
% ylabel('P10 firing (spikes/s)','Fontsize',fontsize-1)
% line([P10data(iPrep).bins(VizP10(iPrep,iProg).Nt(1)),P10data(iPrep).bins(VizP10(iPrep,iProg).Nt(end))] - tStim,[350,350],...
%     'Color',TrainColor,'Linewidth',plotlinewidth)
% line([P10data(iPrep).bins(VizP10(iPrep,iProg).Nt(end)+1),P10data(iPrep).bins(VizP10(iPrep,iProg).forecastNt(end))] - tStim,[350,350],...
%     'Color',FitColor,'Linewidth',plotlinewidth)
% text(P10data(iPrep).bins(VizP10(iPrep,iProg).Nt(1)) - tStim,375,'Training fit','Color',TrainColor,'Fontsize',fontsize-1)
% text(P10data(iPrep).bins(VizP10(iPrep,iProg).Nt(end)+1) - tStim,375,'Forecast (full)','Color',FitColor,'Fontsize',fontsize-1)
% text(65,200,'P10 data','Color',[0 0 0],'Fontsize',fontsize-1)
% set(gca,'FontName','Helvetica','FontSize',fontsize-1);
% set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
% 
% print -depsc Fig4b_bestR2
% 
% exportPPTfig(gcf,'GLM_Best_R2',[10 15 6 3])
% 
% 
% % worst R2
% iProg = 1; iPrep = 2;  % P10 dataset stim 1, prep 2 = Prep 4, stim 1 in full data-set
% load([filepath 'P10 analysis/' fname{iProg} '_StateSpace_P10_GLMmodel']);
% 
% figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3.5 2.1]);
% plot(P10data(iPrep).bins - tStim,P10data(iPrep).spkfcn,'k','Linewidth',plotlinewidth); hold on
% plot(P10data(iPrep).bins(VizP10(iPrep,iProg).forecastNt) - tStim,...
%     VizP10(iPrep,iProg).tsP10Model(numel(VizP10(iPrep,iProg).Nt)+1:end),...
%     'Color',FitColor,'Linewidth',errorlinewidth)
% plot(P10data(iPrep).bins(VizP10(iPrep,iProg).Nt) - tStim,VizP10(iPrep,iProg).tsP10Model(1:numel(VizP10(iPrep,iProg).Nt)),...
%     'Color',TrainColor,'Linewidth',errorlinewidth)
% axis([0 95 0 800])
% ptch = fill([0 2.5 2.5 0],[800 800 0 0],stimcolor);
% set(ptch,'EdgeColor',stimcolor)
% set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap
% 
% set(gca,'YTick',[0 200 400 600 800]);
% xlabel('Time (s)','Fontsize',fontsize-1)
% ylabel('P10 firing (spikes/s)','Fontsize',fontsize-1)
% % line([P10data(iP).bins(VizP10(iP,iF).Nt(1)),P10data(iP).bins(VizP10(iP,iF).Nt(end))],[350,350],...
% %     'Color',TrainColor,'Linewidth',2)
% % line([P10data(iP).bins(VizP10(iP,iF).Nt(end)+1),P10data(iP).bins(VizP10(iP,iF).forecastNt(end))],[350,350],...
% %     'Color',FitColor,'Linewidth',2)
% % text(P10data(iP).bins(VizP10(iP,iF).Nt(1)),365,'Training fit','Color',TrainColor,'Fontsize',fontsize-1)
% % text(P10data(iP).bins(VizP10(iP,iF).Nt(end)+1),365,'Forecast (full)','Color',FitColor,'Fontsize',fontsize-1)
% % text(100,200,'P10 data','Color',[0 0 0],'Fontsize',fontsize-1)
% set(gca,'FontName','Helvetica','FontSize',fontsize-1);
% set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
% 
% print -depsc Fig4b_worstR2
% 
% exportPPTfig(gcf,'GLM_worst_R2',[10 15 6 3])
% 
% 
% % worst MAE
% iProg = 2; iPrep = 2;  % P10 dataset stim 2, prep 2 = Prep 4, stim 2 in full data-set
% load([filepath 'P10 analysis/' fname{iProg} '_StateSpace_P10_GLMmodel']);
% 
% figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3.5 2.1]);
% plot(P10data(iPrep).bins - tStim,P10data(iPrep).spkfcn,'k','Linewidth',plotlinewidth); hold on
% plot(P10data(iPrep).bins(VizP10(iPrep,iProg).forecastNt) - tStim,...
%     VizP10(iPrep,iProg).tsP10Model(numel(VizP10(iPrep,iProg).Nt)+1:end),...
%     'Color',FitColor,'Linewidth',errorlinewidth)
% plot(P10data(iPrep).bins(VizP10(iPrep,iProg).Nt) - tStim,VizP10(iPrep,iProg).tsP10Model(1:numel(VizP10(iPrep,iProg).Nt)),...
%     'Color',TrainColor,'Linewidth',errorlinewidth)
% axis([0 95 0 500])
% ptch = fill([0 2.5 2.5 0],[500 500 0 0],stimcolor);
% set(ptch,'EdgeColor',stimcolor)
% set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap
% 
% set(gca,'YTick',[0 100 200 300 400 500]);
% xlabel('Time (s)','Fontsize',fontsize-1)
% ylabel('P10 firing (spikes/s)','Fontsize',fontsize-1)
% % line([P10data(iP).bins(VizP10(iP,iF).Nt(1)),P10data(iP).bins(VizP10(iP,iF).Nt(end))],[350,350],...
% %     'Color',TrainColor,'Linewidth',2)
% % line([P10data(iP).bins(VizP10(iP,iF).Nt(end)+1),P10data(iP).bins(VizP10(iP,iF).forecastNt(end))],[350,350],...
% %     'Color',FitColor,'Linewidth',2)
% % text(P10data(iP).bins(VizP10(iP,iF).Nt(1)),365,'Training fit','Color',TrainColor,'Fontsize',fontsize-1)
% % text(P10data(iP).bins(VizP10(iP,iF).Nt(end)+1),365,'Forecast (full)','Color',FitColor,'Fontsize',fontsize-1)
% % text(100,200,'P10 data','Color',[0 0 0],'Fontsize',fontsize-1)
% set(gca,'FontName','Helvetica','FontSize',fontsize-1);
% set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
% 
% print -depsc Fig4b_worstMAE
% 
% exportPPTfig(gcf,'GLM_worst_MAE',[10 15 6 3])
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
