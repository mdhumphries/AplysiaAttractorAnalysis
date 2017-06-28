%% script to assemble sub-panels for Figure 3: attractors!
clear all; close all

if ismac
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';
else
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
end

run figure_properties

% M = 4; % need bigger markers here 

%% panel a: coalescence and time on attractor
load([filepath '/Classify attractors/RecurrenceStats_FilterPts'])

npreps = size(RcrStats,1);
nstims = size(RcrStats,2);

% summary plot: 
% (1) color code by prep to show consistency within prep
% (2) link within prep to draw eye to groups
% (3) markersize according to stim sequence: no patterns
% (4) markersize according to "sensitised"

Msize = [2 4 6];  % stim order
% Msize = [5 7 9];  % for PPT

% Msize = [5 15 10];  % does sensitised prep show clear effect?

hF = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on

sym = 'o';
for iR = 1:npreps
    TimeToCoalesce = [RcrStats(iR,:).firstW] - 35;  % seconds to reach coalscence compared to min possible time we can measure
 
    figure(hF)
   plot(100*[RcrStats(iR,:).densAllPeriod],TimeToCoalesce,'o',...
        'Color',prep_Cmap(iR,:),'MarkerSize',min(Msize),'MarkerFaceColor',prep_Cmap(iR,:),'Linewidth',0.5);
    
%     % increasing MarkerSize for sequence
%     for iF = 1:nstims
%         figure(hF)
%         plot(100*RcrStats(iR,iF).densAllPeriod,TimeToCoalesce(iF),sym,...
%             'Color',prep_Cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',prep_Cmap(iR,:));
%         if iF == 3
%             text(100*RcrStats(iR,iF).densAllPeriod,TimeToCoalesce(iF),num2str(iR),'Fontsize',8,'Color',[1 1 1])
%         end
%     end
end

figure(hF)
% line([20 100],[35 35],'Color',[0.5 0.5 0.5],'Linewidth',axlinewidth);
xlabel('Density of recurrence points (%)','FontSize',fontsize)
ylabel('Time to coalscence (s)','FontSize',fontsize)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc Fig3a_coalescence

exportPPTfig(gcf,'Coalescence',[10 15 8 8])

%% panel b: example projections on to same axis?
%% not sure this is useful;

Qt = 0.01; % time-step in ms
w = 201;

% try: (i) shading each line
% (ii) finding good view, then replacing axes as in Fig 2a
% (iii) smoothing...

load([filepath '/Classify attractors/RecurrenceManifoldStats.mat'])

% ixExample = 10;
% cmap = brewermap(3,'Purples');

ixExample = 6;
cmap = brewermap(3,'Reds');

% http://stackoverflow.com/questions/596216/formula-to-determine-brightness-of-rgb-color
L =  0.299*prep_Cmap(ixExample,1) + 0.587*prep_Cmap(ixExample,2) + 0.114*prep_Cmap(ixExample,3);  

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); 
xmin = -20; xmax = 15; ymin = 0; ymax =35; zmin = 0; zmax = 45;

% draw axes lines at the back
axis([xmin xmax ymin ymax zmin zmax])
view([-135 26])

line([xmax,xmax],[ymin,ymin],[zmin,zmax],'Linewidth',axlinewidth,'Color',[0 0 0]); hold on % PC3
line([xmin,xmax],[ymin,ymin],[zmin,zmin],'Linewidth',axlinewidth,'Color',[0 0 0]);  % PC 1  
line([xmax,xmax],[ymin,ymax],[zmin,zmin],'Linewidth',axlinewidth,'Color',[0 0 0]);  % PC2


for iF = 1:3  
    pc1_smth = movingAverage(Viz(ixExample,iF).ProjStimOnda01(:,1),w); 
    pc2_smth = movingAverage(Viz(ixExample,iF).ProjStimOnda01(:,2),w);  
    pc3_smth = movingAverage(Viz(ixExample,iF).ProjStimOnda01(:,3),w); 
  
    pc1_Ctrlsmth = movingAverage(Viz(ixExample,iF).ProjControlStimOnda01(:,1),w); 
    pc2_Ctrlsmth = movingAverage(Viz(ixExample,iF).ProjControlStimOnda01(:,2),w);  
    pc3_Ctrlsmth = movingAverage(Viz(ixExample,iF).ProjControlStimOnda01(:,3),w); 

    prep_stim_color = cmap(iF,:);
    
    %plot3(Viz(ixExample,iF).ProjSponOnda01(:,1),Viz(ixExample,iF).ProjSponOnda01(:,2),Viz(ixExample,iF).ProjSponOnda01(:,3),'Color',[0 0 0]); hold on
    %plot3(Viz(ixExample,iF).ProjSponOnda01(1,1),Viz(ixExample,iF).ProjSponOnda01(1,2),Viz(ixExample,iF).ProjSponOnda01(1,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);        
    % plot3(Viz(ixExample,iF).ProjStimOnda01(:,1),Viz(ixExample,iF).ProjStimOnda01(:,2),Viz(ixExample,iF).ProjStimOnda01(:,3),'Color',prep_Cmap(ixExample,:)); hold on
    plot3(pc1_smth,pc2_smth,pc3_smth,'Color',prep_stim_color); hold on
    plot3(pc1_smth(1),pc2_smth(1),pc3_smth(1),'o','MarkerSize',4,'MarkerFaceColor',prep_stim_color,'MarkerEdgeColor','k'); % prep_stim_color);

    %plot3(Viz(ixExample,iF).ProjControlSponOnda01(:,1),Viz(ixExample,iF).ProjControlSponOnda01(:,2),Viz(ixExample,iF).ProjControlSponOnda01(:,3),'--','Color',[0.6 0.6 0.6]); hold on
    %plot3(Viz(ixExample,iF).ProjControlSponOnda01(1,1),Viz(ixExample,iF).ProjControlSponOnda01(1,2),Viz(ixExample,iF).ProjControlSponOnda01(1,3),'o','MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6]);        
    % plot3(Viz(ixExample,iF).ProjControlStimOnda01(:,1),Viz(ixExample,iF).ProjControlStimOnda01(:,2),Viz(ixExample,iF).ProjControlStimOnda01(:,3),'Color',[0.6 0.6 0.6]); hold on
    % plot3(pc1_Ctrlsmth,pc2_Ctrlsmth,pc3_Ctrlsmth,'Color',prep_Cmap(ixExample,:),'Color',[0.6 0.6 0.6]); hold on

end

text(xmin-3,ymin+5,zmin,'PC1','Fontsize',fontsize)
text(xmax+5,ymax+5,zmin,'PC2','Fontsize',fontsize)
text(xmax,ymin-1,zmax+1,'PC3','Fontsize',fontsize)

axis off

% print -depsc Fig3b_examplevolume
exportfig(gcf,'Fig3b_examplevolume','Color',color,'Format',format,'Resolution',dpi)

% exportPPTfig(gcf,'ExampleVolume',[10 15 8 8])

%% panel c: within same volume - distances between all points projected onto Program 1 axes

load([filepath '/Classify attractors/RecurrenceManifoldStats_AllPreps.mat'])

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
line([0 30],[0 30],'Color',[0 0 0],'Linewidth',axlinewidth)


set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc Fig3c_sameVolume

% exportPPTfig(gcf,'SameVolume_Prog1axes',[10 15 8 8])

%% panel D: initial conditions
load([filepath '/Classify attractors/InitialConditions.mat'])

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on
line([1.5 1.5],[0 80],'Color',[0 0 0],'Linewidth',axlinewidth,'Linestyle',':')
line([2.5 2.5],[0 80],'Color',[0 0 0],'Linewidth',axlinewidth,'Linestyle',':')

LinkedUnivariateScatterPlots(gca,1:3,allDistances,[0.7 0.7 0.7],...
    'Linewidth',plotlinewidth,'MarkerSize',M);

ylabel('Distance between programs (a.u.)','FontSize',fontsize)
axis square

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc FigAttractConditions_InitialConditions

%% panel d: example Sxy pair, and control?
ixExample = 5;

% sequential color map
Sxy_cmap = brewermap(10,'OrRd');

% order neurons by? Total Sxy in first map
SxyTot = sum(Viz(ixExample,2).Sxy);
[srt,Ix] = sort(SxyTot);

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3 3]); hold on
imagesc(Viz(ixExample,2).Sxy(Ix,Ix)); colormap(Sxy_cmap)
axis square;
axis off
% xlabel('Neuron','FontSize',fontsize-1)
% ylabel('Neuron','FontSize',fontsize-1)

print -depsc Fig3d_SxyStim2

exportPPTfig(gcf,'SxyStim2',[10 15 4 4])

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3 3]); hold on
imagesc(Viz(ixExample,3).Sxy(Ix,Ix)); colormap(Sxy_cmap)
axis square
axis off

print -depsc Fig3d_SxyStim3

% need this?
% figure
% imagesc(Viz(ixExample,2).Pxy(Ix,Ix)); colormap(Sxy_cmap)

% exportPPTfig(gcf,'SxyStim3',[10 15 4 4])


%% panel e: Sxy and control....
npreps = size(Prep,2);
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on

for iR = 1:npreps
     plot(Prep(iR).rhoSxy(Prep(iR).rhoSxy>0),Prep(iR).rhoPxy,'o',...
        'Color',prep_Cmap(iR,:),'MarkerSize',M,'MarkerFaceColor',prep_Cmap(iR,:),'Linewidth',0.5);
     % text(Prep(iR).rhoPxy(3),Prep(iR).rhoSxy(2,3),num2str(iR),'Fontsize',8,'Color',[1 1 1])
end

set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1],'YTick',[0 0.2 0.4 0.6 0.8 1])
xlabel('Data correlation','FontSize',fontsize)
ylabel('Control correlation','FontSize',fontsize)
line([0 1],[0 1],'Color',[0 0 0],'Linewidth',axlinewidth)
axis square
axis([0 1 0 1])

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc Fig3e_SxyStability

% exportPPTfig(gcf,'SxyStability',[10 15 8 8])

