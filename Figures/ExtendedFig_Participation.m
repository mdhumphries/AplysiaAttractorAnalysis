%% Extended Figure Participation

clear all; close all

if ispc
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
else
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';
end

run figure_properties
fnames = {'da01','da02','da03'};

% M = 4; % need bigger markers here 

load([filepath '/Classify attractors/ParticipationChanges_SamePCs'])   % participation scores
load([filepath '/Classify attractors/WhatIsParticipation_SamePCs']);  % rates and synchrony

load([filepath '/Classify attractors/RecurrenceManifoldStats_CommonAxes'],'spars','Prep'); % get analysis parameters

npreps = size(Coeffs,2);
nprogs = 3;

pnl_size = [2 2]; % for main figure; [4 4]; for SI Figure

pnl_size = [4 4];
%% does participation mean just rate?

% example rate versus Participation
iPrep = 1;
iProg = 1;

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 pnl_size]); hold on
plot(PartContribution(iPrep).rates(iProg,:)',100*Coeffs(iPrep).normWMax(iProg,:)','ko','Markersize',M)
axis square
axis([0 5 0 100])
set(gca,'XTick',[0:1:5])
text(0.5,95,sprintf('R^2=%.2g',PartContribution(iPrep).rRatePart(iProg)^2),'Fontsize',fontsize-1)
% ylabel('Participation (% of max)','FontSize',fontsize)
% xlabel('Firing rate (spike/s)','FontSize',fontsize)
ylabel('Participation (%)','FontSize',fontsize)
xlabel('Rate (spike/s)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFigPart_rateScatter_example

%exportPPTfig(gcf,'Example_Participation_Rate_Scatter',[10 15 6 6])

%% histogram of all R^2 for rate vs participation

for iPrep = 1:npreps
    for iProg = 1:nprogs
        allR2(iPrep,iProg) = PartContribution(iPrep).rRatePart(iProg)^2;
    end
end
binstep = 0.05;
binsR2 = 0:binstep:1;
hR2 = hist(allR2(:),binsR2);

mAllR2 = mean(allR2(:));
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 pnl_size]); hold on
bar(binsR2+binstep/2,hR2,1,'FaceColor',[1 1 1],'LineWidth',axlinewidth)
line([mAllR2 mAllR2],[0 7],'Color',[0.8 0.3 0.3],'Linewidth',axlinewidth)
axis([0 1 0 7])
axis square

% xlabel('R^2: participation','FontSize',fontsize)
% ylabel('Number of programs','FontSize',fontsize)
xlabel('R^2: rate','FontSize',fontsize)
ylabel('No. programs','FontSize',fontsize)
text(mAllR2+0.02,7.2,sprintf('%.2g',mAllR2),'Color',[0.8 0.3 0.3],'FontSize',fontsize-1);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

%exportPPTfig(gcf,'Participation_Rate_R_Hist',[10 15 6 4])

print -depsc ExtFigPart_rateR2Histogram

%% change in firing rate vs change in participation: example
iPrep = 2;
iChange = 1;

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 pnl_size]); hold on
plot(PartContribution(iPrep).dRates(:,iChange),100*Coeffs(iPrep).dMaxP(:,iChange),'ko','Markersize',M)
plot(PartContribution(iPrep).dRates(Coeffs(iPrep).VarPart,iChange),100*Coeffs(iPrep).dMaxP(Coeffs(iPrep).VarPart,iChange),'ro','Markersize',M)

axis square
axis([-1 2 -30 30])
set(gca,'XTick',[-1:1:2])
text(1,-25,sprintf('R^2=%.2g',PartContribution(iPrep).rDrateDpart(iChange)^2),'Fontsize',fontsize-1)
% ylabel('Participation (% of max)','FontSize',fontsize)
% xlabel('Firing rate (spike/s)','FontSize',fontsize)
ylabel('Change in participation (%)','FontSize',fontsize)
xlabel('Change in rate (spike/s)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFigPart_Change_rateScatter_example


%% histogram of all R^2 for change in rate vs change in participation

for iPrep = 1:npreps
    for iChange = 1:nprogs
        allR2change(iPrep,iChange) = PartContribution(iPrep).rDrateDpart(iChange)^2;
    end
end
binsR2 = 0:binstep:1;
hR2 = hist(allR2change(:),binsR2);

mAllR2 = mean(allR2change(:));
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 pnl_size]); hold on
bar(binsR2+binstep/2,hR2,1,'FaceColor',[1 1 1],'LineWidth',axlinewidth)
line([mAllR2 mAllR2],[0 7],'Color',[0.8 0.3 0.3],'Linewidth',axlinewidth)
axis([0 1 0 7])
axis square

% xlabel('R^2: participation','FontSize',fontsize)
% ylabel('Number of programs','FontSize',fontsize)
xlabel('R^2: change','FontSize',fontsize)
ylabel('No. programs','FontSize',fontsize)
text(mAllR2+0.02,7.2,sprintf('%.2g',mAllR2),'Color',[0.8 0.3 0.3],'FontSize',fontsize-1);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

% %% just maximum change in both
% for iPrep = 1:npreps
%     for iChange = 1:nprogs
%         allR2max(iPrep) = PartContribution(iPrep).rDrateDpart(iChange)^2;
%     end
% end
% binsR2 = 0:0.05:1;
% hR2 = hist(allR2change(:),binsR2);

print -depsc ExtFigPart_ChangeRateR2Histogram

%% AGAIN: for synchrony

% example rate versus Participation
iPrep = 4;
iProg = 1;

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 pnl_size]); hold on
plot(PartContribution(iPrep).TSxy(iProg,:)',100*Coeffs(iPrep).normWMax(iProg,:)','ko','Markersize',M)
axis square
axis([0 85 0 100])
set(gca,'XTick',[0:20:80])
text(5,95,sprintf('R^2=%.2g',PartContribution(iPrep).rSynchPart(iProg)^2),'Fontsize',fontsize-1)
% ylabel('Participation (% of max)','FontSize',fontsize)
% xlabel('Firing rate (spike/s)','FontSize',fontsize)
ylabel('Participation (%)','FontSize',fontsize)
% xlabel('Total synchrony ($\sum |Cxy|$)','FontSize',fontsize,'Interpreter','latex')
xlabel('Total synchrony (sum |Cxy|)','FontSize',fontsize);

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFigPart_synchScatter_example

%% histogram of R^2 for synch
for iPrep = 1:npreps
    for iProg = 1:nprogs
        allR2(iPrep,iProg) = PartContribution(iPrep).rSynchPart(iProg)^2;
    end
end
binstep = 0.05;
binsR2 = 0:binstep:1;
hR2 = hist(allR2(:),binsR2);

mAllR2 = mean(allR2(:));
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 pnl_size]); hold on
bar(binsR2+binstep/2,hR2,1,'FaceColor',[1 1 1],'LineWidth',axlinewidth)
line([mAllR2 mAllR2],[0 15],'Color',[0.8 0.3 0.3],'Linewidth',axlinewidth)
axis([0 1 0 15])
axis square

% xlabel('R^2: participation','FontSize',fontsize)
% ylabel('Number of programs','FontSize',fontsize)
xlabel('R^2: total synchrony','FontSize',fontsize)
ylabel('No. programs','FontSize',fontsize)
text(mAllR2+0.02,7.2,sprintf('%.2g',mAllR2),'Color',[0.8 0.3 0.3],'FontSize',fontsize-1);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

%exportPPTfig(gcf,'Participation_Rate_R_Hist',[10 15 6 4])

print -depsc ExtFigPart_synchR2Histogram

%% change in synch vs change in participation: example
iPrep = 6;
iChange = 1;

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 pnl_size]); hold on
plot(PartContribution(iPrep).dSynch(:,iChange),100*Coeffs(iPrep).dMaxP(:,iChange),'ko','Markersize',M)
plot(PartContribution(iPrep).dSynch(Coeffs(iPrep).VarPart,iChange),100*Coeffs(iPrep).dMaxP(Coeffs(iPrep).VarPart,iChange),'ro','Markersize',M)

axis square
axis([-20 10 -25 25])
set(gca,'XTick',[-20:10:10])
text(-18,20,sprintf('R^2=%.2g',PartContribution(iPrep).rDsynchDpart(iChange)^2),'Fontsize',fontsize-1)
% ylabel('Participation (% of max)','FontSize',fontsize)
% xlabel('Firing rate (spike/s)','FontSize',fontsize)
ylabel('Change in participation (%)','FontSize',fontsize)
xlabel('Change in total synchrony','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFigPart_Change_synchScatter_example

%% %% histogram of all R^2 for change in synch vs change in participation

for iPrep = 1:npreps
    for iChange = 1:nprogs
        allR2change(iPrep,iChange) = PartContribution(iPrep).rDsynchDpart(iChange)^2;
    end
end
binsR2 = 0:binstep:1;
hR2 = hist(allR2change(:),binsR2);

mAllR2 = mean(allR2change(:));
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 pnl_size]); hold on
bar(binsR2+binstep/2,hR2,1,'FaceColor',[1 1 1],'LineWidth',axlinewidth)
line([mAllR2 mAllR2],[0 9],'Color',[0.8 0.3 0.3],'Linewidth',axlinewidth)
axis([0 1 0 9])
axis square

% xlabel('R^2: participation','FontSize',fontsize)
% ylabel('Number of programs','FontSize',fontsize)
xlabel('R^2: change total synchrony','FontSize',fontsize)
ylabel('No. programs','FontSize',fontsize)
text(mAllR2+0.02,9,sprintf('%.2g',mAllR2),'Color',[0.8 0.3 0.3],'FontSize',fontsize-1);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFigPart_ChangeSynchR2Histogram
