%% script to assemble sub-panels for Extended Figure showing quality of fit...
% NB, now using first panel in main figure of P10 decoding
clear all; close all

if ispc
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
else
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';
end

run figure_properties

M = 4; % need bigger markers here

fname = {'da01','da02','da03'}; 


% load overall P10 stats
load([filepath 'P10 analysis/P10_statespace_Stats'])

TrainColor = [0.4 0.4 0.8];
FitColor = [0.8 0.4 0.4];

% FitColor = [153,50,204] / 256;

MAEcolor = [0 0 0];
Rcolor = [0.4 0.8 0.4];

ErrorColor = [0.6 0.6 0.6];

stimcolor = [0.7 0.7 0.7];
tStim = 30;

fit_Scatter = [2 2]; % inset size; [4 4] for main SI Figure
fontsize = fontsize - 1;

%% panel: summary scatter of fit

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 fit_Scatter]); hold on; 

for iP = 1:3
    for iF = 1:3
        % MAE error bar
        semMAE = SummaryStats(iP,iF).Train.stdMAE / sqrt(numel(VizP10(iP,iF).tsForecastStart));
        line([SummaryStats(iP,iF).Train.meanMAE-2*semMAE,SummaryStats(iP,iF).Train.meanMAE+2*semMAE],[SummaryStats(iP,iF).Train.meanR2 SummaryStats(iP,iF).Train.meanR2],'Linewidth',errorlinewidth,'Color',ErrorColor)
        % R error bar
        semR2 = SummaryStats(iP,iF).Train.stdR2 / sqrt(numel(VizP10(iP,iF).tsForecastStart));
        line([SummaryStats(iP,iF).Train.meanMAE,SummaryStats(iP,iF).Train.meanMAE],[SummaryStats(iP,iF).Train.meanR2-2*semR2 SummaryStats(iP,iF).Train.meanR2+2*semR2],'Linewidth',errorlinewidth,'Color',ErrorColor)
%         % MAE std bar
%         line([SummaryStats(iP,iF).Train.meanMAE-SummaryStats(iP,iF).Train.stdMAE,SummaryStats(iP,iF).Train.meanMAE+SummaryStats(iP,iF).Train.stdMAE],[SummaryStats(iP,iF).Train.meanR2 SummaryStats(iP,iF).Train.meanR2],'Linewidth',errorlinewidth,'Color',ErrorColor)
%         % R std bar
%         line([SummaryStats(iP,iF).Train.meanMAE,SummaryStats(iP,iF).Train.meanMAE],[SummaryStats(iP,iF).Train.meanR2-SummaryStats(iP,iF).Train.stdR2 SummaryStats(iP,iF).Train.meanR2+SummaryStats(iP,iF).Train.stdR2],'Linewidth',errorlinewidth,'Color',ErrorColor)
        
        % data-point
        plot([SummaryStats(iP,iF).Train.meanMAE],[SummaryStats(iP,iF).Train.meanR2],'o',...
            'MarkerSize',2,'MarkerFaceColor',TrainColor,'MarkerEdgeColor',TrainColor); hold on

    end
end

% figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3 3]); hold on;
% % plot([SummaryStats(:).meanMAE],[SummaryStats(:).meanR].^2,'r.'); hold on
% plot([SummaryStats(:).medianMAE],[SummaryStats(:).medianR].^2,'k.','MarkerSize',M);
axis([0 40 0 1])
axis square
xlabel('MAE (spikes/s)','FontSize',fontsize); 
ylabel('R^2','FontSize',fontsize)

% set(gca,'XTick',[0 10 20 30])
set(gca,'XTick',[0 10 20 30 40])

%exportPPTfig(gcf,['Summary_R2_vs_MAE'],[10 15 5 5])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');


print -depsc ExtFigP10_fit_error_summary

% exportPPTfig(gcf,'P10GLM_summary',[10 15 8 8])


%% compare to forecasts: predicitve power?
allMAE = []; allR2 = [];
for iP = 1:3
    for iF = 1:3
            allMAE = [allMAE; SummaryStats(iP,iF).Train.meanMAE SummaryStats(iP,iF).meanMAE];
            allR2 = [allR2; SummaryStats(iP,iF).Train.meanR2 SummaryStats(iP,iF).meanR2];           
    end
end

% simple: compare means

% MAE
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on;
line([0 40],[0 40],'Linewidth',axlinewidth,'Color',[0 0 0]);
plot(allMAE(:,1),allMAE(:,2),'o',...
    'MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k'); hold on
axis square
axis([0 40 0 40])
xlabel('Fit: MAE (spikes/s)','FontSize',fontsize); 
ylabel('Forecast: MAE (spikes/s)','FontSize',fontsize); 

% set(gca,'XTick',[0 10 20 30 40])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFigP10_fitforecast_MAE

% R2
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on;
line([0 40],[0 40],'Linewidth',axlinewidth,'Color',[0 0 0]);
plot(allR2(:,1),allR2(:,2),'o',...
    'MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k'); hold on
axis square
axis([0 1 0 1])
xlabel('Fit: R^2','FontSize',fontsize); 
ylabel('Forecast: R^2','FontSize',fontsize); 

set(gca,'XTick',0:0.2:1,'YTick',0:0.2:1)

%exportPPTfig(gcf,['Summary_R2_vs_MAE'],[10 15 5 5])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
print -depsc ExtFigP10_fitforecast_R2

% complex: compare within-program over all windows, then compare between programs (but likely strongly correlates with
% means, right?)






