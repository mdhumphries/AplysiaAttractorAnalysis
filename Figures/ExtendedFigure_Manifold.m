% Extended Data Figure: more on same manifold...

clear all; close all
if ismac
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';
else
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
end
run figure_properties

%% panel a: project on to axes from concatenated time series
load([filepath '/Classify attractors/RecurrenceManifoldStats_CommonAxes_AllPreps.mat'])

npreps = size(Prep,2);

hD = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on

for iR = 1:npreps
     % plot Hausdorff distance
     if ~isempty(Prep(iR).HaussD)
         ctrl = mean(Prep(iR).HaussDCtrl,3); ctrl = ctrl(ctrl > 0); 
         sem = std(Prep(iR).HaussDCtrl,0,3) / sqrt(spars.nShuffles); sem = sem(sem > 0);
         data = Prep(iR).HaussD(Prep(iR).HaussD > 0);
         figure(hD)
         plot(data,ctrl,'o',...
            'Color',prep_Cmap(iR,:),'MarkerSize',M,'MarkerFaceColor',prep_Cmap(iR,:),'Linewidth',0.5);
         line([data data]',[ctrl-2*sem ctrl+2*sem]','Linewidth',errorlinewidth,'Color',prep_Cmap(iR,:))
         % text(prctDcontrol(2,3),prctDdata(2,3),num2str(iR),'Fontsize',8,'Color',[1 1 1])
     end
end

% title('Percentage points on same manifold')
xlabel('Data: distance between programs (a.u.)','FontSize',fontsize)
ylabel('Control: distance between programs (a.u.)','FontSize',fontsize)
axis square
axis([0 30 0 30])
line([0 30],[0 30],'Color',[0.6 0.6 0.6],'Linewidth',axlinewidth)


set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFigManifold_Concatenated_SameVolume

% exportPPTfig(gcf,'SameVolume_JointAxes',[10 15 8 8])

%% outliers

% Outlier from distance is Prep 4, program 1 (from pairs (1,2) and (1,3))
% Outlier from correlation is Prep 7, program 3  (from pairs (1,3) and (2,3))

load([filepath 'da03_DataProperties_FunctionAndWindowSize'],'PCAdata','Data','SDF')
point_bins = Data(1).stimbins-30;

exampleIDs = [4,1; 7,3];  % Recurrence plots for these outliers

ixT = 3; % max threshold
    
for iEx = 1:size(exampleIDs,1)

    load([filepath '/Classify attractors/da0' num2str(exampleIDs(iEx,2)) '_RecurrencePlotData.mat'])

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]);
    imagesc(point_bins,point_bins,RecurData(exampleIDs(iEx,1)).RecurPlot(ixT).Rp); hold on
    colormap(gray)

    set(gca,'XTick',[0 30 60 90],'YTick',[0 30 60 90])

    xlabel('Time post-stimulation (s)','Fontsize',fontsize)
    ylabel('Time post-stimulation (s)','Fontsize',fontsize)

    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

    %print -depsc Fig2d_recurplot
    exportfig(gcf,['ExampleRP_Prep_' num2str(exampleIDs(iEx,1)) '_Stim_' num2str(exampleIDs(iEx,2))],'Color',color,'Format',format,'Resolution',dpi)
end

