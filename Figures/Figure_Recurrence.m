%% script to assemble sub-panels for Figure 2: raster and rates
clear all; close all

if ismac
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';
else
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
end

run figure_properties

M = 4; % need bigger markers here 

ixExample = 5; % prep #5, stim 3 used throughout

ixT = 3;  % which threshold? s

%% panel a: PCA projection of two recordings
Qt = 0.01; % time-step in ms
w = 201;
stimstart = 30;

stimbin = stimstart / Qt - (w-1)/2;

load([filepath 'da02_DataProperties_FunctionAndWindowSize'],'PCAdata','Data','SDF')

ixts = 1:stimbin;
ixStim = stimbin+1:numel(Data(ixExample).bins);

da02.sponPCprojU = zeros(numel(ixts),3);  % spon PCs on stim axes
da02.stimPCprojU = zeros(numel(ixStim),3);  % stim PCs on stim axes
for iP = 1:3 
    da02.sponPCprojU(:,iP) = sum(bsxfun(@times,SDF(ixExample).spkfcn(ixts,:),PCAdata(ixExample).coeffs(:,iP)'),2); % unnormalised 
    da02.stimPCprojU(:,iP) = sum(bsxfun(@times,SDF(ixExample).spkfcn(ixStim,:),PCAdata(ixExample).coeffs(:,iP)'),2); % unnormalised 
    da02.allPCprojU(:,iP) = sum(bsxfun(@times,SDF(ixExample).spkfcn,PCAdata(ixExample).coeffs(:,iP)'),2); % unnormalised 
end


load([filepath 'da03_DataProperties_FunctionAndWindowSize'],'PCAdata','Data','SDF')
da03_bins = Data(ixExample).bins;
ixts = 1:stimbin;
ixStim = stimbin+1:numel(da03_bins);

% do projection onto unnormalised PC axes
da03.sponPCprojU = zeros(numel(ixts),3);  % spon PCs on stim axes
da03.stimPCprojU = zeros(numel(ixStim),3);  % stim PCs on stim axes
for iP = 1:3 
    da03.sponPCprojU(:,iP) = sum(bsxfun(@times,SDF(ixExample).spkfcn(ixts,:),PCAdata(ixExample).coeffs(:,iP)'),2); % unnormalised 
    da03.stimPCprojU(:,iP) = sum(bsxfun(@times,SDF(ixExample).spkfcn(ixStim,:),PCAdata(ixExample).coeffs(:,iP)'),2); % unnormalised 
    da03.allPCprojU(:,iP) = sum(bsxfun(@times,SDF(ixExample).spkfcn,PCAdata(ixExample).coeffs(:,iP)'),2); % unnormalised 
end

%[r,~] = size(SDFSens(iR).spkfcn);

     
da02_1_smth = movingAverage(da02.allPCprojU(:,1),w); da02_2_smth = movingAverage(da02.allPCprojU(:,2),w);  da02_3_smth = movingAverage(da02.allPCprojU(:,3),w); 
da03_1_smth = movingAverage(da03.allPCprojU(:,1),w); da03_2_smth = movingAverage(da03.allPCprojU(:,2),w);  da03_3_smth = movingAverage(da03.allPCprojU(:,3),w); 

%% make figure panel
xmin = -7; xmax = 38; ymin = -15; ymax =15; zmin = -2; zmax = 30;

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]);
% draw axes lines at the back
axis([xmin xmax ymin ymax zmin zmax])
view([-74,14]); 
% view([-74,25])
xlabel('x')
ylabel('y')

line([xmax,xmax],[ymax,ymax],[zmin,zmax],'Linewidth',axlinewidth,'Color',[0 0 0]); hold on
line([xmin,xmax],[ymax,ymax],[zmin,zmin],'Linewidth',axlinewidth,'Color',[0 0 0])
line([xmax,xmax],[ymin,ymax],[zmin,zmin],'Linewidth',axlinewidth,'Color',[0 0 0]); 

%plot3(da02_1_smth(1:stimbin),da02_2_smth(1:stimbin),da02_3_smth(1:stimbin),'Color',sponcolor,'Linewidth',pclinewidth);  hold on
plot3(da03_1_smth(1:stimbin),da03_2_smth(1:stimbin),da03_3_smth(1:stimbin),'Color',sponcolor,'Linewidth',pclinewidth); hold on
%plot3(da02_1_smth(1),da02_2_smth(1),da02_3_smth(1),'o','MarkerFaceColor',[0 0 0],'MarkerSize',M,'MarkerEdgeColor',[0 0 0]);
plot3(da03_1_smth(1),da03_2_smth(1),da03_3_smth(1),'o','MarkerFaceColor',[0 0 0],'MarkerSize',M,'MarkerEdgeColor',[0 0 0]);

% plot3(da02_1_smth(stimbin+1:end),da02_2_smth(stimbin+1:end),da02_3_smth(stimbin+1:end),'Color',da02color,'Linewidth',pclinewidth);  
plot3(da03_1_smth(stimbin+1:end),da03_2_smth(stimbin+1:end),da03_3_smth(stimbin+1:end),'Color',prep_Cmap(ixExample,:),'Linewidth',pclinewidth);  
plot3(da03_1_smth(stimbin),da03_2_smth(stimbin),da03_3_smth(stimbin),'o','MarkerFaceColor',prep_Cmap(ixExample,:),'MarkerSize',M,'MarkerEdgeColor',prep_Cmap(ixExample,:));


text(xmin-3,ymax+5,zmin,'PC1','Fontsize',fontsize)
text(xmax-1,ymin-1,zmin,'PC2','Fontsize',fontsize)
text(xmax,ymax-1,zmax+1,'PC3','Fontsize',fontsize)

axis off

% print -depsc Fig2a_PCproj
exportfig(gcf,'Fig2a_PCproj','Color',color,'Format',format,'Resolution',dpi)

% exportPPTfig(gcf,'Fig2s_PCproj')


%% panel b cumulative distributions of %Variance...

load([filepath '/da01_DataProperties_FunctionAndWindowSize.mat'],'PCAdata')
nPreps = numel(PCAdata);
for i = 1:nPreps
    load([filepath '/da01_DataProperties_FunctionAndWindowSize.mat'],'PCAdata')
    PCA(i,1).prctVar = PCAdata(i).prctVar;
    load([filepath '/da02_DataProperties_FunctionAndWindowSize.mat'],'PCAdata')
    PCA(i,2).prctVar = PCAdata(i).prctVar;
    load([filepath '/da03_DataProperties_FunctionAndWindowSize.mat'],'PCAdata')
    PCA(i,3).prctVar = PCAdata(i).prctVar;
end

% 
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3.5 3.5]); hold on
line([0.75 200],[80 80],'Linewidth',axlinewidth,'Color',[0.8 0.3 0.3])
for iP = 1:nPreps
    for iS = 1:3
        nPC = numel(PCA(iP,iS).prctVar);
        semilogx(1:nPC,100*PCA(iP,iS).prctVar,'.-','Linewidth',0.25,'Color',stim_set(iS,:));
    end
end
set(gca,'YLim',[0 100],'XLim',[0.75 200],'XTick',[1 10 100],'XTickLabel',[1 10 100],'XScale','log')
xlabel('Embedding dimensions','Fontsize',fontsize)
ylabel('Variance explained (%)','Fontsize',fontsize)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc FigRecurrence_CumulVar


%% panel b inset: histogram of number PCs

load([filepath '/Classify attractors/da01_RecurrencePointsData.mat'])
da01_nPCs = [RcrPtData(:).maxD];

load([filepath '/Classify attractors/da02_RecurrencePointsData.mat'])
da02_nPCs = [RcrPtData(:).maxD];

load([filepath '/Classify attractors/da03_RecurrencePointsData.mat'])
da03_nPCs = [RcrPtData(:).maxD];

bins = 1:1:15;
hPCs = hist([da01_nPCs da02_nPCs da03_nPCs],bins);

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3 2]);
bar(bins,hPCs,'FaceColor',[0 0 0])
axis([1 10 0 10])
pos = get(gca,'Position');
set(gca,'Position',[pos(1)+0.1 pos(2)+0.1 0.5 0.4])
set(gca,'XTick',[2:2:10])
xlabel('Embedding dimensions','Fontsize',fontsize-1)
ylabel('Num. recordings','Fontsize',fontsize-1)

set(gca,'FontName','Helvetica','FontSize',fontsize-1);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc Fig2b_histogram

%exportPPTfig(gcf,'Fig2s_histogram')

%% panel d: recurrence plot example

% recurrence plot
load([filepath '/Classify attractors/da03_RecurrencePlotData.mat'])

ixT = 3; % max threshold

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]);
imagesc(Data(ixExample).stimbins-stimstart,Data(ixExample).stimbins-stimstart,RecurData(ixExample).RecurPlot(ixT).Rp); hold on
colormap(gray)
% line([0 da03_bins(end)],[30 30],'Linewidth',axlinewidth,'Color',stimcolor);
% line([30 30],[0 da03_bins(end)],'Linewidth',axlinewidth,'Color',stimcolor);

set(gca,'XTick',[0 40 80 120],'YTick',[0 40 80 120])

xlabel('Time post-stimulation (s)','Fontsize',fontsize)
ylabel('Time post-stimulation (s)','Fontsize',fontsize)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

%print -depsc Fig2d_recurplot
exportfig(gcf,'Fig2d_recurplot','Color',color,'Format',format,'Resolution',dpi)

%exportPPTfig(gcf,'Fig2_recurplot',[10 15 6 6])

%% panel X: theta histogram example: robustness to threshold

% all Theta histograms
cmap = brewermap(10,'OrRd');

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2.5]);
imagesc(pars.bins,pars.prctTheta,RcrPtData(ixExample).hstRecurrence);
colormap(cmap)

% axis([1 10 0 10])
set(gca,'YTick',[2,6,10],'YTickLabel',pars.prctTheta)

xlabel('Recurrence delay (s)','Fontsize',fontsize)
ylabel('Threshold (%)','Fontsize',fontsize)

pos=get(gca,'pos');
set(gca,'pos',[pos(1) pos(2) pos(3) pos(4)-0.075]);
pos=get(gca,'pos');
hc = colorbar('location','northoutside','position',[pos(1)+0.1 pos(2)+pos(4)+0.03 pos(3)*0.75 0.02])
set(hc,'xaxisloc','top','XTick',0:0.1:0.3,'LineWidth',axlinewidth)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(hc,'Fontsize',fontsize-1);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

%print -depsc Fig2d_thetahistogram
exportfig(gcf,'Fig2e_thetahistogram','Color',color,'Format',format,'Resolution',dpi,'FontMode','scaled','FontSize',1)

%exportPPTfig(gcf,'Fig2_thetahistogram',[10 15 6 4])

%% panel x: recurrence period estimates for all programs

for ixT = 1:3
    load([filepath '/Classify attractors/RecurrenceStats_FilterPts'])
    % (1) plot all histograms as one big panel
    % build combined matrix
    % sort by size of peak?
    allmaxP = zeros(30,1); allmaxPfltr = zeros(30,1); allMtime = zeros(30,1);
    allHsts = [];
    ctr = 1;
    for iF  = 1:size(RcrStats,1)
        load([filepath '/Classify attractors/da01_RecurrencePointsData_FilterPts.mat'])
        allHsts = [allHsts RcrPtData(iF).hstRecurrence(ixT,:)']; % max theta; column array
        % max peak in histogram: do we want this for only after 5s? 
        allmaxP(ctr) = max(RcrPtData(iF).hstRecurrence(ixT,:));     
        allmaxPfltr(ctr) = max(RcrPtData(iF).hstRecurrence(ixT,pars.bins>=pars.tmin));      
        nPoints = [RcrPtData(iF).Orbit(ixT).Period(:).nPoints];
        ixPrd = nPoints == max(nPoints);  % dominant period
        allMtime(ctr) =  RcrPtData(iF).Orbit(ixT).Period(ixPrd).mTime;
        allRecur(ctr) = 1 - RcrStats(iF,1).densNoRecur; % proportion of points that are recurrent
        allPropInPeriod(ctr) = RcrStats(iF,1).propInMaxPeriod;  % proportion of recurrent points in main period
        ctr=ctr+1;

        load([filepath '/Classify attractors/da02_RecurrencePointsData_FilterPts.mat'])
        allHsts = [allHsts RcrPtData(iF).hstRecurrence(end,:)']; % max theta; column array
        allmaxP(ctr) = max(RcrPtData(iF).hstRecurrence(end,:)); 
        allmaxPfltr(ctr) = max(RcrPtData(iF).hstRecurrence(ixT,pars.bins>=pars.tmin));      
        nPoints = [RcrPtData(iF).Orbit(ixT).Period(:).nPoints];
        ixPrd = nPoints == max(nPoints);  % dominant period
        allMtime(ctr) =  RcrPtData(iF).Orbit(ixT).Period(ixPrd).mTime;
        allRecur(ctr) = 1 - RcrStats(iF,2).densNoRecur; % proportion of points that are recurrent
        allPropInPeriod(ctr) = RcrStats(iF,2).propInMaxPeriod;  % proportion of recurrent points in main period
        ctr=ctr+1;

        load([filepath '/Classify attractors/da03_RecurrencePointsData_FilterPts.mat'])
        allHsts = [allHsts RcrPtData(iF).hstRecurrence(end,:)']; % max theta; column array
        allmaxP(ctr) = max(RcrPtData(iF).hstRecurrence(end,:)); 
        allmaxPfltr(ctr) = max(RcrPtData(iF).hstRecurrence(ixT,pars.bins>=pars.tmin));      
        nPoints = [RcrPtData(iF).Orbit(ixT).Period(:).nPoints];
        ixPrd = nPoints == max(nPoints);  % dominant period
        allMtime(ctr) =  RcrPtData(iF).Orbit(ixT).Period(ixPrd).mTime;
        allRecur(ctr) = 1 - RcrStats(iF,3).densNoRecur; % proportion of points that are recurrent
        allPropInPeriod(ctr) = RcrStats(iF,3).propInMaxPeriod;  % proportion of recurrent points in main period
        ctr=ctr+1;

    end

    [srt,IxP] = sort(allmaxPfltr,'descend');  % sort by size of peak: concentration of recurrent points
    % [srt,IxP] = sort(allRecur,'descend');  % sort by proportion of points that are recurrent
    % [srt,IxP] = sort(allPropInPeriod,'descend');  % sort by proportion of recurrent points in main Period

    cmap = brewermap(15,'OrRd');

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 4.5]);
    %imagesc(1:numel(allMtime),pars.bins,allHsts(:,IxP));
    %set(gca,'YDir','normal');
    imagesc(pars.bins,1:numel(allMtime),allHsts(:,IxP)');
    colormap(cmap)
    line([pars.tmin pars.tmin],[0.5 30.5],'Linewidth',axlinewidth,'Color',[0.5 0.5 0.5])

    xlabel('Recurrence delay (s)','Fontsize',fontsize)
    %ylabel('Rate of contraction (real eigenvalue)','FontSize',fontsize)
    ylabel('Program No.','FontSize',fontsize)

    % add colorbar over the top of the plot
    pos=get(gca,'pos');
    set(gca,'pos',[pos(1) pos(2) pos(3) pos(4)-0.075]);
    pos=get(gca,'pos');
    hc = colorbar('location','northoutside','position',[pos(1)+0.1 pos(2)+pos(4)+0.03 pos(3)*0.75 0.02]);
    set(hc,'xaxisloc','top','XTick',0:0.1:0.3)

    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

    % keyboard

    %print -depsc Fig2_rotation_period_histogram
    exportfig(gcf,['Fig2_rotation_period_histogram_Threshold' num2str(pars.prctTheta(ixT))],'Color',color,'Format',format,'Resolution',dpi)

    % (1b) histograms of density of recurrent points & density in peak - sorted as above
    cmap = brewermap(10,'Blues');
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 2 4.5]);
    imagesc(1:2,30,[allRecur(IxP)' allPropInPeriod(IxP)']);
    colormap(cmap)
    denspos = get(gca,'Position');
    set(gca,'Position',[0.25 0.2 0.2 pos(4)])
    axis off
    hc = colorbar('location','eastoutside','position',[0.5 0.2 0.025 pos(4)*0.75]);
    set(hc,'YTick',0:0.2:0.8)
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

    % print -depsc Fig2_recurrence_density
    exportfig(gcf,['Fig2_recurrence_density_Threshold' num2str(pars.prctTheta(ixT))],'Color',color,'Format',format,'Resolution',dpi)
    % exportfig(gcf,'Fig2_recurrence_density','Color',color,'Format','preview','Resolution',dpi)


    % (2) scatter of mean time of dominant period
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 2 4]); hold on;
    hs = scatter(ones(30,1),allMtime,'MarkerEdgeColor',[0.5 0.5 0.5],'SizeData',5); hold on

    % line([0.9 1.1],[mean(D_pre_postN),mean(D_pre_postN)],'Color','k','Linewidth',errorlinewidth)
    % line([1 1],[mean(D_pre_postN)-2*std(D_pre_postN)/sqrt(n),mean(D_pre_postN)+2*std(D_pre_postN)/sqrt(n)],'Color','k','Linewidth',errorlinewidth)

    line([0.75 1.25],[pars.tmin pars.tmin],'Color',[0.5 0.5 0.5],'Linewidth',axlinewidth)
    axis([0.75 1.25 0 30])

    ylabel('Mean period (s)','FontSize',fontsize); 
    set(gca,'XTick',[],'XColor','w')
    pos = get(gca,'Position');
    set(gca,'Position',[pos(1) pos(2) 0.3 pos(4)])

    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);

    print(['Fig2_rotation_period_scatter_Threshold' num2str(pars.prctTheta(ixT))],'-depsc')

    % exportPPTfig(gcf,'Period_Scatter',[10 15 3 6])
end

