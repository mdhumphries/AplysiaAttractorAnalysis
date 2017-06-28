%% script to assemble sub-panels for Figure X: single neuron contribution...
clear all; close all

%filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';

run figure_properties

load([filepath '/Classify attractors/RobustnessNeuron_vs_Program_Separate'])

%% panel a: distance between participation distributions - SI?

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on
line([0.1 0.45],[1 1],'Color','k','Linewidth',axlinewidth); hold on
plot(Robust.allH,Robust.allHD_MagCtrl,'o','Color',[0.8 0.2 0.2],'Markersize',M); 
plot(Robust.allH(Robust.allHD_MagCtrl > 1),Robust.allHD_MagCtrl(Robust.allHD_MagCtrl > 1),'ko','Markersize',M); 
axis square
axis([0.1 0.45 0.2 1.2])  
ylabel('Distance between programs (norm.)','FontSize',fontsize)
xlabel('P(participation) distance','FontSize',fontsize);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

% ADD r^2 here...
Robust.lm_NormDist_M
text(0.375,0.9,'R^2 = 0.06','FontSize',fontsize,'Color',[0.8 0.2 0.2])  % separate axes

% exportPPTfig(gcf,'ProgDist_vs_ParticipationDist',[10 15 6 6])

print -depsc FigVar_Distance

%% panel b: total participation change - main figure
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on
line([4 24],[1 1],'Color','k'); hold on
plot(Robust.allDpart,Robust.allHD_MagCtrl,'o','Color',[0.8 0.2 0.2],'Markersize',M); 
plot(Robust.allDpart(Robust.allHD_MagCtrl > 1),Robust.allHD_MagCtrl(Robust.allHD_MagCtrl > 1),'ko','Markersize',M); 
axis square
axis([4 24 0.2 1.2])  
ylabel('Distance between programs (norm.)','FontSize',fontsize)
xlabel('Total participation change (a.u.)','FontSize',fontsize)

% update R^2 from model:
Robust.lm_ChangePart_M
text(18,0.9,'R^2 = 0.03','FontSize',fontsize,'Color',[0.8 0.2 0.2])  % separate axes

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc FigVar_TotalChange

% exportPPTfig(gcf,'ProgDist_vs_TotalPartChange',[10 15 6 6])

%% panel c: proportion of variable neurons within that pair of programs
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on
line([0 0.3],[1 1],'Color','k'); hold on
plot(Robust.allPchange,Robust.allHD_MagCtrl,'o','Color',[0.8 0.2 0.2],'Markersize',M); 
plot(Robust.allPchange(Robust.allHD_MagCtrl > 1),Robust.allHD_MagCtrl(Robust.allHD_MagCtrl > 1),'ko','Markersize',M); 
axis square
axis([0 0.3 0.2 1.2])  
ylabel('Distance between programs (norm.)','FontSize',fontsize)
xlabel('Proportion of variable neurons','FontSize',fontsize)

Robust.lm_Pchange_M
text(0.2,0.9,'R^2 = 0.035','FontSize',fontsize,'Color',[0.8 0.2 0.2])  % separate axes

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

% print -depsc FigVar_PropNeurons


% exportPPTfig(gcf,'ProgDist_vs_NeuronProp',[10 15 6 6])
