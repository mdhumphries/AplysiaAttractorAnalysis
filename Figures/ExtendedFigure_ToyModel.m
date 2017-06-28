% Extended Figure? Toy model of emergent oscillation, and attractors

clear all; close all

if ispc
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
else
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';
end

run figure_properties

M = 4; % need bigger markers here 

load([filepath '/Classify attractors/Toy3Dmodel'])
time = sim.dt:sim.dt:sim.T;

clr1 = [0.8 0.6 0.2];
clr2 = [0.6 0.3 0.8];
clr3 = [0.3 0.8 0.6];

ratesize = [10 15 5 2];

% indices of each model
Basic = 1:3;
IC2 = 4:6;
IC3 = 7:9;
Noise1 = 10:12;
Noise2 = 13:15;
Noise3 = 16:18;
TransP = 19:21;
SusP = 22:24;
Decay = 25:27;
Unconnected = 28:30;  

%% panel ?: response of model without connections - neurons are not oscillators

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2]); hold on; 
plot(time,r(Unconnected(1),:),'Color',clr1,'Linewidth',plotlinewidth); hold on
plot(time,r(Unconnected(2),:),'Color',clr2,'Linewidth',plotlinewidth);
plot(time,r(Unconnected(3),:),'Color',clr3,'Linewidth',plotlinewidth);
axis([0 3 0 3])
xlabel('Time (s)','FontSize',fontsize); 
ylabel('Rate (Hz)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_ToyOsc_Rates_Unconnected


%% panel b: basic attractor: 3 neurons and PCA

coeffs1 = pca(r(Basic,:)');
pc1 = r(Basic,:)' * coeffs1(:,1);  % reconstruct PC1
pc2 = r(Basic,:)' * coeffs1(:,2);  % reconstruct PC2

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',ratesize); hold on; 
plot(time,r(Basic(1),:),'Color',clr1,'Linewidth',plotlinewidth); hold on
plot(time,r(Basic(2),:),'Color',clr2,'Linewidth',plotlinewidth);
plot(time,r(Basic(3),:),'Color',clr3,'Linewidth',plotlinewidth);
axis([0 3 0 3])
xlabel('Time (s)','FontSize',fontsize); 
ylabel('Rate (Hz)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_ToyOsc_basic_rates

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on; 
plot(pc1(2:end),pc2(2:end),'k','Linewidth',errorlinewidth); hold on
plot(pc1(2),pc2(2),'k.','Markersize',10)
axis square
xlabel('PC 1','FontSize',fontsize)
ylabel('PC 2','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_ToyOsc_basic_PCA

%% panel c: convergence

coeffs2 = pca(r(IC2,:)');
pc1IC2 = r(IC2,:)' * coeffs2(:,1);  % reconstruct PC1
pc2IC2 = r(4:6,:)' * coeffs2(:,2);  % reconstruct PC2

coeffs3 = pca(r(IC3,:)');
pc1IC3 = r(IC3,:)' * coeffs3(:,1);  % reconstruct PC1
pc2IC3 = r(IC3,:)' * coeffs3(:,2);  % reconstruct PC2

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',ratesize); hold on; 
plot(time,r(IC2(1),:),'Color',clr1,'Linewidth',plotlinewidth); hold on
plot(time,r(IC2(2),:),'Color',clr2,'Linewidth',plotlinewidth);
plot(time,r(IC2(3),:),'Color',clr3,'Linewidth',plotlinewidth);
axis([0 3 0 3])
xlabel('Time (s)','FontSize',fontsize); 
ylabel('Rate (Hz)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_ToyOsc_IC2_rates


figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',ratesize); hold on; 
plot(time,r(IC3(1),:),'Color',clr1,'Linewidth',plotlinewidth); hold on
plot(time,r(IC3(2),:),'Color',clr2,'Linewidth',plotlinewidth);
plot(time,r(IC3(3),:),'Color',clr3,'Linewidth',plotlinewidth);
axis([0 3 0 3])
xlabel('Time (s)','FontSize',fontsize); 
ylabel('Rate (Hz)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_ToyOsc_IC3_rates


figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 5]); hold on; 
plot(pc1(2:end),pc2(2:end),'k','Linewidth',errorlinewidth); hold on
plot(pc1(2),pc2(2),'k.','Markersize',10)
plot(pc1IC2(2:end),pc2IC2(2:end),'r','Linewidth',errorlinewidth); hold on
plot(pc1IC2(2),pc2IC2(2),'r.','Markersize',10)
plot(pc1IC3(2:end),pc2IC3(2:end),'b','Linewidth',errorlinewidth); hold on
plot(pc1IC3(2),pc2IC3(2),'b.','Markersize',10)
axis square
xlabel('PC 1','FontSize',fontsize)
ylabel('PC 2','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_ToyOsc_3ICs_PCA

%% panel d:  transient perturbation
coeffsP = pca(r(TransP,:)');
pc1P = r(TransP,:)' * coeffsP(:,1);  % reconstruct PC1
pc2P = r(TransP,:)' * coeffsP(:,2);  % reconstruct PC2

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',ratesize); hold on; 
plot(time,r(TransP(1),:),'Color',clr1,'Linewidth',plotlinewidth); hold on
plot(time,r(TransP(2),:),'Color',clr2,'Linewidth',plotlinewidth);
plot(time,r(TransP(3),:),'Color',clr3,'Linewidth',plotlinewidth);
line([1.5 1.6],[3.1 3.1],'Color',[0 0 0],'Linewidth',2)
axis([0 3.1 0 3.1])
xlabel('Time (s)','FontSize',fontsize); 
ylabel('Rate (Hz)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_ToyOsc_TransPerturb_rates

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on; 
plot(pc1P(2:sim.tsperturb(1)),pc2P(2:sim.tsperturb(1)),'k','Linewidth',errorlinewidth); hold on
plot(pc1P(2),pc2P(2),'k.','Markersize',10)
plot(pc1P(sim.tsperturb),pc2P(sim.tsperturb),'r','Linewidth',errorlinewidth)
plot(pc1P(sim.tsperturb(end):end),pc2P(sim.tsperturb(end):end),'Color',[0 0 1],'Linewidth',errorlinewidth); hold on
axis square
xlabel('PC 1','FontSize',fontsize)
ylabel('PC 2','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_ToyOsc_TransPerturb_PC


%% panel e: sustained perturbation
coeffsPP = pca(r(SusP,:)');
pc1PP = r(SusP,:)' * coeffsPP(:,1);  % reconstruct PC1
pc2PP = r(SusP,:)' * coeffsPP(:,2);  % reconstruct PC2


figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',ratesize); hold on; 
plot(time,r(SusP(1),:),'Color',clr1,'Linewidth',plotlinewidth); hold on
plot(time,r(SusP(2),:),'Color',clr2,'Linewidth',plotlinewidth);
plot(time,r(SusP(3),:),'Color',clr3,'Linewidth',plotlinewidth);
line([1.5 3],[3.1 3.1],'Color',[0 0 0],'Linewidth',2)
axis([0 3.1 0 3.1])
xlabel('Time (s)','FontSize',fontsize); 
ylabel('Rate (Hz)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_ToyOsc_SustainedPerturb_rates

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on; 
plot(pc1PP(2:1500),pc2PP(2:1500),'k','Linewidth',errorlinewidth); hold on
plot(pc1PP(2),pc2PP(2),'k.','Markersize',10)
plot(pc1PP(1500:end),pc2PP(1500:end),'Color',[0 0 1],'Linewidth',errorlinewidth); hold on
axis square
xlabel('PC 1','FontSize',fontsize)
ylabel('PC 2','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_ToyOsc_SustainedPerturb_PC


%% panel f: spiral
coeffsDec = pca(r(Decay,:)');
pc1Dec = r(Decay,:)' * coeffsDec(:,1);  % reconstruct PC1
pc2Dec = r(Decay,:)' * coeffsDec(:,2);  % reconstruct PC2

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',ratesize); hold on; 
plot(time,r(Decay(1),:),'Color',clr1,'Linewidth',plotlinewidth); hold on
plot(time,r(Decay(2),:),'Color',clr2,'Linewidth',plotlinewidth);
plot(time,r(Decay(3),:),'Color',clr3,'Linewidth',plotlinewidth);
axis([0 3 0 3])
xlabel('Time (s)','FontSize',fontsize); 
ylabel('Rate (Hz)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_ToyOsc_Spiral_rates

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on; 
plot(pc1Dec(2:end),pc2Dec(2:end),'k','Linewidth',errorlinewidth); hold on
plot(pc1Dec(2),pc2Dec(2),'k.','Markersize',10)
axis square
axis([-2 3 -2 3])
xlabel('PC 1','FontSize',fontsize)
ylabel('PC 2','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_ToyOsc_Spiral_PC
