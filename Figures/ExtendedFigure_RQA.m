% Extended Figure? showing complexity of programs using recurrence
% quantification analysis
%
% (1) plot example histograms of line lengths
% (2) plot scatter of trapping and prediction time
% (3) plot scatter of prediction time vs complexity
% (4) plot at least 3 example recurrence plots across that range

clear all; close all

%filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';

run figure_properties

M = 4; % need bigger markers here 

ixExample = 5; % prep #5, stims 2 and 3

load([filepath '/Classify attractors/RecurPlotStats'])


%% panel a: example histogram of lengths
load([filepath '/Classify attractors/da02_RecurrencePlotData'])

L = RecurData(ixExample).RecurPlot(RPpars.ixT).dlengths .* 0.01; % .* Data(iR).GaussQt; % line lengths in seconds 

bins = min(L):0.01:max(L);
hL = hist(L,bins);

% raw histogram of line counts
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3 5]); hold on; 
bar(bins,hL./sum(hL),'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
line([RPpars.lminT RPpars.lminT],[0 0.01],'Linewidth',0.5,'Color',[0.8 0.2 0.4])
axis([0 30 0 0.01])
xlabel('Duration of line (s)','FontSize',fontsize); 
ylabel('P(duration)','FontSize',fontsize)
set(gca,'XTick',[0 10 20 30],'YTick',[0 0.005 0.01]);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

RPstats(ixExample,2).EntL

% % histogram of cumulative line counts - basis for all calculations
% figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on;
% bar(bins,RPstats(ixExample,2).Hl,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
% line([RPpars.lminT RPpars.lminT],[0.001 1e7],'Linewidth',0.5,'Color',[0.8 0.2 0.4])
% set(gca,'YScale','log')
% axis([0 30 min(RPstats(ixExample,2).Hd) 1e7])
% xlabel('Duration of line (s)','FontSize',fontsize);
% ylabel('Number of lines $\geq$ duration','FontSize',fontsize)
% set(gca,'XTick',[0 10 20 30],'YTick',[1e0 1e2 1e4 1e6]);
% set(gca,'FontName','Helvetica','FontSize',fontsize);
% set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_RQA_Example_histogram

%% panel b: scatters of program quantities

npreps = size(RPstats,1);
nstims = size(RPstats,2);

% Msize = [2 4 6];  % stim order
Msize = [5 10 15];  % stim order

hT = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 8 8]); hold on; 
hE = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 8 8]); hold on; 
hDV = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 8 8]); hold on; 


join ='o';  %'o-' to link by lines within prep
sym = 'o';
for iR = 1:npreps
    % compare diagonal and vertical trapping time
    figure(hDV)
%     subplot(211),plot([RPstats(iR,:).mPredTime],[RPstats(iR,:).mTrapTime],join,...
%         'Color',cmap(iR,:),'MarkerSize',2,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5); hold on
    plot([RPstats(iR,:).Lmax],[RPstats(iR,:).Vmax],join,...   
        'Color',prep_Cmap(iR,:),'MarkerSize',M,'MarkerFaceColor',prep_Cmap(iR,:),'Linewidth',0.5); hold on
    
    % compare max time on attractor with determinism of movement
    figure(hT)
    % plot([RPstats(iR,:).mPredTimeL],[RPstats(iR,:).DetL],join,...  
    plot([RPstats(iR,:).Lmax],[RPstats(iR,:).DetL],join,...   
        'Color',prep_Cmap(iR,:),'MarkerSize',M,'MarkerFaceColor',prep_Cmap(iR,:),'Linewidth',0.5);
    
    figure(hE) 
    % plot([RPstats(iR,:).mPredTime],[RPstats(iR,:).Ent],join,...
    plot([RPstats(iR,:).Lmax],[RPstats(iR,:).EntL],join,...
        'Color',prep_Cmap(iR,:),'MarkerSize',M,'MarkerFaceColor',prep_Cmap(iR,:),'Linewidth',0.5);
    
   
%      % Uncomment for increasing MarkerSize for sequence & labelled by
%      % number
%      for iF = 1:nstims
%         figure(hDV)
% %         subplot(211),plot(RPstats(iR,iF).mPredTime,RPstats(iR,iF).mTrapTime,sym,...
% %             'Color',cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',cmap(iR,:)); 
% %         if iF == 3
% %             text([RPstats(iR,iF).mPredTime],[RPstats(iR,iF).mTrapTime],num2str(iR),'Fontsize',8,'Color',[1 1 1])
% %         end
% 
%         % subplot(212)
%         plot([RPstats(iR,iF).Lmax],[RPstats(iR,iF).Vmax],join,...   
%             'Color',prep_Cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',prep_Cmap(iR,:)); 
%         if iF == 3
%             text([RPstats(iR,iF).Lmax],[RPstats(iR,iF).Vmax],num2str(iR),'Fontsize',8,'Color',[1 1 1])
%         end
%         
%         figure(hT)
%         % plot([RPstats(iR,:).Lmax],[RPstats(iR,:).DetL],join,...  
%         
%         plot([RPstats(iR,iF).Lmax],[RPstats(iR,iF).DetL],sym,...     
%              'Color',prep_Cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',prep_Cmap(iR,:));
%         if iF == 3
%             % text([RPstats(iR,iF).mPredTimeL],[RPstats(iR,iF).DetL],num2str(iR),'Fontsize',8,'Color',[1 1 1])
%             text([RPstats(iR,iF).Lmax],[RPstats(iR,iF).DetL],num2str(iR),'Fontsize',8,'Color',[1 1 1])
%         end       
%         
%     end
end
figure(hDV)
% subplot(211), xlabel('Mean prediction time (s)'); ylabel('Mean trapping time (s)')
%subplot(212), 
xlabel('Max prediction time (s)','FontSize',fontsize); 
ylabel('Max trapping time (s)','FontSize',fontsize)

% axis([0 50 -60 20]);  % recording is 90 seconds long, and longest coalescence is about 40 seconds so max period is 45 seconds
% line([0 50],[0 0],'Color',[0.5 0.5 0.5])


set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');


figure(hT)
xlabel('Max prediction time (s)','FontSize',fontsize); 
% ylabel(['Points in lines > ' num2str(RPpars.lminT) ' s (%)'],'FontSize',fontsize)
ylabel(['Proportion of recurrent points in lines of at least ' num2str(RPpars.lminT) ' s'],'FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
print -depsc ExtFig_RQA_MaxPred_Det


figure(hE)
xlabel('Max prediction time (s)','FontSize',fontsize)
ylabel('Entropy (bits)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

%% panel insets: example recurrence plots for illustrative cases

load([filepath 'da03_DataProperties_FunctionAndWindowSize'],'PCAdata','Data','SDF')
point_bins = Data(ixExample).stimbins - 30;

%exampleIDs = [7,3; 9,2; 3,2; 10,3; 5,2];  % (prep,stim)
% exampleIDs = [7,3; 3,2; 4,2; 5,2];  % (prep,stim) for max prediction time
exampleIDs = [7,3; 3,2; 10,3; 4,3; 5,2];  % (prep,stim) for max prediction time, Lmin = 5s

% exampleIDs = [7,3; 2,3; 6,1; 4,3];  % mean prediction time, Lmin = 1s

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


