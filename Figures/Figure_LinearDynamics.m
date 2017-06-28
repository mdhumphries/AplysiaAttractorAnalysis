% Figure of all eigenvalue results

clear all; close all

if ismac
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';
    spikepath = '/Users/mqbssmhg/Dropbox/paper2/datasets/deconcatenated_spks/experimental/';
else
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
    spikepath = 'C:\Users\mqbssmhg.DS\Dropbox\paper2\datasets\deconcatenated_spks\experimental\';
 % temp for Dell laptop
%     filepath = 'C:\Users\Mark\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
%     spikepath = 'C:\Users\Mark\Dropbox\SpikeData\Angela Bruno Aplysia Dynamics Data\datasets\deconcatenated_spks\experimental\';

end

run figure_properties

M = 4; % need bigger markers here 

ixExample = 5; % prep #5, stims 2 and 3

ixT = 3;  % which threshold? s

load([filepath '/Classify attractors/RecurrenceStats_FilterPts'])
npreps = size(RcrStats,1);
nstims = size(RcrStats,2);


%% panel 0: check that all programs had dominant periods: most points were in max period
% maxDom = zeros(npreps,nstims);
% for iR = 1:npreps
%     for iS = 1:nstims
%         maxDom(iR,iS) = max(RcrStats(iR,iS).propComparePeriods);
%     end
% end
% 
% bins = 0:5:100;
% figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2.5]); hold on
% h = histogram(100*[maxDom(:)],bins);
% set(h,'FaceColor',[0.6 0.6 0.6]);
% axis([50 100 0 20])
% 
% xlabel('Delays in largest peak (of all peaks) (%)','FontSize',fontsize)
% ylabel('No. programs','FontSize',fontsize)
% 
% set(gca,'FontName','Helvetica','FontSize',fontsize);
% set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');



%% panel a: example time-series of rotation and contraction from same example program
ixExample = 9;
load([filepath 'da03_DataProperties_FunctionAndWindowSize'],'Data')

load([filepath '/Classify attractors/da03_RecurrencePointsData_FilterPts'])

% dominant period
nPoints = [RcrPtData(ixExample).Orbit(statpars.ixT).Period(:).nPoints];
ixPrd = find(nPoints == max(nPoints));  
usedPts = RcrPtData(ixExample).Orbit(statpars.ixT).Period(ixPrd).Npts >= pars.minNeighbours;

% % time-series of raw points
% figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2]); hold on
% hR = plot(Data(ixExample).stimbins(RcrPtData(ixExample).ix0 + RcrPtData(ixExample).Orbit(statpars.ixT).Period(ixPrd).IDs),...
%     RcrPtData(ixExample).Orbit(statpars.ixT).Period(ixPrd).RotPeriod,'.');
% set(hR,'Markersize',2,'Color',prep_Cmap(ixExample,:));
% axis([30 120 0 20])
% xlabel('Time (s)','FontSize',fontsize)
% ylabel('Period (s)','FontSize',fontsize)
% set(gca,'FontName','Helvetica','FontSize',fontsize);
% set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
% 
% % print -depsc Fig2_rotation_timeseries
% 
% figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2]); hold on
% line([30 120],[0 0],'Color',[0 0 0],'Linewidth',axlinewidth)
% hE = plot(Data(ixExample).stimbins(RcrPtData(ixExample).ix0 + RcrPtData(ixExample).Orbit(statpars.ixT).Period(ixPrd).IDs),...
%     real(RcrPtData(ixExample).Orbit(statpars.ixT).Period(ixPrd).Emax),'.');
% set(hE,'Markersize',2,'Color',prep_Cmap(ixExample,:));
% axis([30 120 -0.02 0.02])
% xlabel('Time (s)','FontSize',fontsize)
% ylabel('Real Eig.','FontSize',fontsize)
% set(gca,'FontName','Helvetica','FontSize',fontsize);
% set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

% plot as sliding window time-series instead

% imaginary component or rotation?
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2.5]); hold on
% hR =  plot(Data(ixExample).stimbins(RcrPtData(ixExample).winmids),RcrPtData(ixExample).Orbit(statpars.ixT).Period(ixPrd).mRotPeriod,'o')
hR =  plot(Data(ixExample).stimbins(RcrPtData(ixExample).winmids) - Data(ixExample).stimbins(1),RcrPtData(ixExample).Orbit(statpars.ixT).Period(ixPrd).mEigI,'o')

set(hR,'Markersize',2,'Color',prep_Cmap(ixExample,:));
axis([0 90 0 0.02])
ptch = fill([0 2.5 2.5 0],[0.02 0.02 0 0],stimcolor);
set(ptch,'EdgeColor',stimcolor)
set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap

xlabel('Time (s)','FontSize',fontsize)
% ylabel('Period (s)','FontSize',fontsize)
ylabel('Imag. Eig.','FontSize',fontsize)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

% print -depsc Fig2_rotation_timeseries
% exportPPTfig(gcf,'Rotation_timeseries',[10 15 6 3])
print -depsc FigLinDyn_EigI_timeseries


% real eigenvalue
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2.5]); hold on
ptch = fill([0 2.5 2.5 0],[0.005 0.005 -0.02 -0.02],stimcolor);
set(ptch,'EdgeColor',stimcolor)
set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap

line([0 90],[0 0],'Color',[0 0 0],'Linewidth',axlinewidth)
hE =  plot(Data(ixExample).stimbins(RcrPtData(ixExample).winmids) - Data(ixExample).stimbins(1),RcrPtData(ixExample).Orbit(statpars.ixT).Period(ixPrd).mEigR,'o')
set(hE,'Markersize',2,'Color',prep_Cmap(ixExample,:));
axis([0 90 -0.02 0.005])
xlabel('Time (s)','FontSize',fontsize)
ylabel('Real Eig.','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc FigLinDyn_EigR_timeseries
% exportPPTfig(gcf,'EigR_timeseries',[10 15 6 3])

%% b: summary panel: plot just eigenvalues on plane
load([filepath '/Classify attractors/RecurrenceStats_FilterPts'])

npreps = size(RcrStats,1);
nstims = size(RcrStats,2);

Msize = [2 4 6];  % stim order
sym = 'o';


xmin = -0.0005; 
% complex component goes from 0 to infinity, for slow to fast rotation
% so xmax is fastest possible rotation = rotation in the minimum length of
% points submitted to dot(x) = Ax model
% 100 points = 1 second
xmax = 0.01*2*pi/1; % max imaginary eigenvalue: fastest possible rotation we can detect
xmax = 0.015;
ymin = -0.01; 
ymax = abs(ymin); %0.0025;

hm = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 5]); hold on
line([xmin xmax],[0 0],'Color',[0 0 0],'Linewidth',axlinewidth)

% do patches
patchinfo.Vertices = [xmin ymin; abs(xmin) ymin; abs(xmin) ymax; xmin ymax];
patchinfo.Faces = [1 2 3 4];
patchinfo.FaceColor = [0.7 0.7 0.7];
patch(patchinfo,'EdgeColor',[0.7 0.7 0.7]);

% draw horizontal width as edge of 2*SEMS that cross zero?
patchinfo.Vertices = [xmin 0.0001; xmax 0.0001; xmax -0.0001; xmin -0.0001];
patchinfo.Faces = [1 2 3 4];
patchinfo.FaceColor = [0.7 0.7 0.7];
patch(patchinfo,'EdgeColor',[0.7 0.7 0.7]);


for iR = 1:npreps
    % ID P10 preps using squares
%     if any(iR == P10preps')
%         sym = 's'; 
%     else
%         sym = 'o';
%     end
    
    semR = [RcrStats(iR,:).sdEigR] ./ sqrt([RcrStats(iR,:).nUsedPts]); 
    semI = [RcrStats(iR,:).sdEigI] ./ sqrt([RcrStats(iR,:).nUsedPts]); 
    
    figure(hm)
    plot([RcrStats(iR,:).mEigI],[RcrStats(iR,:).mEigR],'o',...
        'Color',prep_Cmap(iR,:),'MarkerSize',2,'MarkerFaceColor',prep_Cmap(iR,:),'Linewidth',0.5);

%     line([[RcrStats(iR,:).mEigI]; [RcrStats(iR,:).mEigI]],...
%         [[[RcrStats(iR,:).mEigR]-[RcrStats(iR,:).sdEigR]]; [[RcrStats(iR,:).mEigR]+[RcrStats(iR,:).sdEigR]]],'Color',prep_Cmap(iR,:))
    
    line([[RcrStats(iR,:).mEigI]; [RcrStats(iR,:).mEigI]],...
        [[RcrStats(iR,:).mEigR]-2*semR; [RcrStats(iR,:).mEigR]+2*semR],'Color',prep_Cmap(iR,:));

    line([[RcrStats(iR,:).mEigI]-2*semI; [RcrStats(iR,:).mEigI]+2*semI],...
        [[RcrStats(iR,:).mEigR]; [RcrStats(iR,:).mEigR]],'Color',prep_Cmap(iR,:));

    % increasing MarkerSize for sequence
    for iF = 1:nstims
        figure(hm)
        plot(RcrStats(iR,iF).mEigI,RcrStats(iR,iF).mEigR,sym,...
            'Color',prep_Cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',prep_Cmap(iR,:));
        % label one node of each
        if iF == 3
            text(RcrStats(iR,iF).mEigI,RcrStats(iR,iF).mEigR,num2str(iR),'Fontsize',8,'Color',[1 1 1])
        end
    end
end

set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap

% plot mean of all programs
allmeanR = [RcrStats.mEigR];
allmeanI = [RcrStats.mEigI];

plot(mean(allmeanI),mean(allmeanR),'o','Color',[0 0 0],'MarkerSize',3,'MarkerFaceColor',[0 0 0])

% scale axes
axis([xmin xmax  ymin ymax])

set(gca,'XTick',[0 0.005 0.01 0.015],'XTickLabel',[0 0.005 0.01 0.015])
set(gca,'YTick',[-0.01:0.0025:0.0025],'YTickLabel',[-0.01:0.0025:0.0025])
set(gca,'YTick',[-0.01:0.005:0.01],'YTickLabel',[-0.01:0.005:0.01])

xlabel('Eigenvalue: imaginary','FontSize',fontsize)
ylabel('Eigenvalue: real','FontSize',fontsize)


% keyboard

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc FigLinDyn_scatter_eigenvalues
% exportPPTfig(gcf,'Prep_eigenvalues',[10 15 8 8])


%% panel c: re-expressed as contraction and rotation
% summary plot: 
% (1) color code by prep to show consistency within prep
% (2) link within prep to draw eye to groups
% (3) markersize according to stim sequence: no patterns
% (4) markersize according to "sensitised"

Msize = [2 4 6];  % stim order
%Msize = [5 7 9];  % for PPT

% Msize = [5 15 10];  % does sensitised prep show clear effect?

hm = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 5]); hold on

sym = 'o';
join = 'o';
allContract = [];
allRotate = [];
for iR = 1:npreps
    Contract = exp(100*[RcrStats(iR,:).mEigR]) * 100 - 100;  % percentage per second
    % stdContract = exp(100*[RcrStats(iR,:).sdEigR]) * 100 - 100;  % percentage per second
    Rotate = 2*pi ./ [RcrStats(iR,:).mEigI] * 0.01;
    
    allContract = [allContract; Contract'];
    allRotate = [allRotate; Rotate'];
    
    figure(hm)
    plot(Rotate,Contract,join,...
        'Color',prep_Cmap(iR,:),'MarkerSize',2,'MarkerFaceColor',prep_Cmap(iR,:),'Linewidth',0.5);
%     line([[RcrStats(iR,:).mRot]; [RcrStats(iR,:).mRot]],...
%         [[Contract-stdContract]; [Contract+stdContract]],'Color',prep_Cmap(iR,:))

%     % increasing MarkerSize for sequence
%     for iF = 1:nstims
%         figure(hm)
%         plot(Rotate(iF),Contract(iF),sym,...
%             'Color',prep_Cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',prep_Cmap(iR,:));
%         % label one node of each
%         if iF == 3
%             text([RcrStats(iR,iF).mRot],Contract(iF),num2str(iR),'Fontsize',8,'Color',[1 1 1])
%         end
%     end
end


% add global mean
plot(mean(allRotate),mean(allContract),'o','Color',[0 0 0],'MarkerSize',3,'MarkerFaceColor',[0 0 0])

% % linear model estimates
axis([0 45 -60 60]);  % recording is 90 seconds long, and longest coalescence is about 40 seconds so max period is 45 seconds
line([0 45],[0 0],'Color',[0.5 0.5 0.5])
xlabel('Orbit rotation period (s)','FontSize',fontsize)
ylabel('Rate of change (%/s)','FontSize',fontsize)


set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc FigLinDyn_scatter_eigenvalues_converted
%exportPPTfig(gcf,'Prep_properties',[10 15 8 8])

%% panel d: change in rotation period during program: example
load([filepath 'da03_DataProperties_FunctionAndWindowSize'],'Data')

load([filepath '/Classify attractors/da03_RecurrencePointsData_FilterPts'])
load([filepath '/Classify attractors/PeriodChange'])

ixExample = 3;
% get recurrence density
tsDensity = RcrPtData(ixExample).Orbit(statpars.ixT).AllRecurrentPoint.pRecur;

% get time windows
ixW = RcrPtData(ixExample).winmids;
tsW = Data(ixExample).stimbins(ixW);

% find admissible windows
ixwins = 1:numel(ixW);
ixWinStrt = 1; ixWinEnd = numel(ixwins);
if oscpars.blnCoalesce  % then only check after coalescence
    ixC = find(Data(ixExample).stimbins == RcrStats(ixExample,3).firstW);
    ixWinStrt = find(ixW == ixC);
end
if oscpars.blnOff
    ixWinEnd = find(tsDensity > oscpars.pDensity,1,'last');  % final window with P% density of recurrence
end
ixwins = ixWinStrt:ixWinEnd; 
DelayWall = RcrPtData(ixExample).Orbit(statpars.ixT).AllRecurrentPoint.mDelay;

% dominant period
nPoints = [RcrPtData(ixExample).Orbit(statpars.ixT).Period(:).nPoints];
ixPrd = find(nPoints == max(nPoints));  


figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 3]); hold on
% hR = plot(Data(ixExample).stimbins(RcrPtData(ixExample).winmids) - Data(ixExample).stimbins(1),RcrPtData(ixExample).Orbit(statpars.ixT).Period(ixPrd).mRotPeriod,'o');
hR = plot(tsW(ixwins)-30,DelayWall(ixwins),'o');

set(hR,'Markersize',2,'Color',prep_Cmap(ixExample,:));
set(gca,'YTick',[0 5 10 15])
axis([0 90 0 17])
% HERE: FIX....
text(70,5,{'\rho = 0.73','P < 10^{-5}'},'FontSize',fontsize) 
line([RcrStats(ixExample,3).firstW RcrStats(ixExample,3).firstW]-30,[0 17],'Color',[0.7 0.7 0.7])

ptch = fill([0 2.5 2.5 0],[17 17 0 0],stimcolor);
set(ptch,'EdgeColor',stimcolor)
set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap
% hl = line([0 2.5],[15.5 15.5],'Linewidth',1.5,'Color',stimcolor,'Clipping','off');
text(3.5,17.5,'stimulus','Fontsize',fontsize-1,'Color',stimcolor)
xlabel('Time (s)','FontSize',fontsize)
ylabel('Orbital period (s)','FontSize',fontsize)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

% print -depsc Fig2_rotation_timeseries
% exportPPTfig(gcf,'Rotation_timeseries',[10 15 6 3])
print -depsc FigLinDyn_Rot_Timeseries_Example

%% panel e: the change in period over time

alpha = 0.01;
bins = -1:0.1:1;

%r = [PeriodChange.rho]; p =[PeriodChange.pval];
r = [PeriodChange.WrhoRecurAll]; p =[PeriodChange.WpvalRecurAll];  % weighted Spearman's rank

hNon = histc(r(p > alpha),bins);
hYes = histc(r(p <= alpha),bins);

hm = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 3]); hold on

hno = bar(bins,hNon,'histc'); hold on
hyes = bar(bins,hYes,'histc'); hold on
%shading flat
set(hno,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[1 1 1]);
set(hyes,'FaceColor',[0.8 0.2 0.2],'EdgeColor',[1 1 1]);
line([-1 1],[0 0],'Color',[0 0 0],'Linewidth',axlinewidth)
axis([-1 1 0 6])

xlabel('Correlation','FontSize',fontsize)
%ylabel('Rate of contraction (real eigenvalue)','FontSize',fontsize)
ylabel('No. programs','FontSize',fontsize)


set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc FigLinDyn_Rho_Period_vs_Time

%% panel f: correlation of change in period and mean contraction...

allmeanR = [RcrStats.mEigR];
r = [PeriodChange.WrhoRecurAll];

[RperiodEig,PperiodEig] = corr(r',allmeanR');

hm = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3 3]); hold on
plot(r,allmeanR,'ko','MarkerSize',2)
axis([-1 1 -1e-2 2e-3])

xlabel('Correlation (time vs period)','FontSize',fontsize)
ylabel('Eigenvalue: real','FontSize',fontsize)


set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

text(0.4,2e-3,['R = ' sprintf('%.2g',RperiodEig)],'FontSize',fontsize-1)

print -depsc FigLinDyn_Rho_vs_Eig


% %% examples that do not change....
% 
% rMat = reshape(r,npreps,nstims); pMat = reshape(p,npreps,nstims);
% 
% [iNonePrep,iNoneStim] = ind2sub([npreps,nstims], find(pMat > alpha));
% 
% for iN = 1:numel(iNonePrep)
%     strStim = ['da0' num2str(iNoneStim(iN))];
% 
%     % get basic data
%     load([filepath strStim '_DataProperties_FunctionAndWindowSize'],'Data','FileTable','DataTable')
%     % get spikes
%     load([spikepath FileTable{iNonePrep(iN)}]);
%     
%     % times
%     startts = 0; % floor(DataTable(iR,2));
%     % endts = floor(DataTable(ixP,3));
%     tStim = Data(iNonePrep(iN)).stimbins(1);
%     tEnd = floor(DataTable(iNonePrep(iN),3));
% 
%     % all spike-train IDs
%     allIDs = unique(spks(:,1));
% 
%     % load clusters
%     load([filepath '/Ensembles/' strStim '_StaticEnsembles.mat'],'StaticEnsemblesAll')
%     n = numel(StaticEnsemblesAll(iNonePrep(iN)).IDs);    
%     Sxy = StaticEnsemblesAll(iNonePrep(iN)).Cxy; Sxy(Sxy < 0) = 0;
%     Sxy(eye(n)==1) = 0;
%     
%     [NG,Sg,Sin,Sout] = sortbysimilarity([StaticEnsemblesAll(iNonePrep(iN)).IDs StaticEnsemblesAll(iNonePrep(iN)).grps],Sxy);
%     
%     % plot them
%     figure
%     [H,GS,I,O] = plot_clusters(spks,NG,StaticEnsemblesAll(iNonePrep(iN)).ngrps,[0 tEnd],'B3',[],M);
%     set(H,'Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 5.5 7]); 
%     title(['Prep ' num2str(iNonePrep(iN)) ', Stim: ' num2str(iNoneStim(iN)) ': \rho = ' num2str(rMat(iNonePrep(iN),iNoneStim(iN)))]);
%     line([110,120],[-1,-1],'Linewidth',plotlinewidth,'Color',[0 0 0],'Clipping','off')
%     text(115,-3,'10 s','FontSize',fontsize)
%     axis([tStim, max(get(gca,'XLim')), 0, n+1])
%     axis off
%     
%     pause
%  end



