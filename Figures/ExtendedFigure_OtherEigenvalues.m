%% figure panels with exciting spiral(s)

clear all; close all;
if ispc
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
    spikepath = 'C:\Users\mqbssmhg.DS\Dropbox\SpikeData\Angela Bruno Aplysia Dynamics Data\datasets\deconcatenated_spks\experimental\';
else
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';

end

load([filepath '/Ensembles/Rates_StaticEnsembles'])
load([filepath '/Classify attractors/CheckAllEigs.mat'])

iW = 3; % 20 s windows
Rlimit = 0.2; % to be considered changing...
[nPreps,nStims] = size(EnsembleRates);

run figure_properties
% M = 4; % need bigger markers here 

%% get all proportions up/no change/down

for iPrep = 1:nPreps
    for iStim = 1:nStims
        % sort into change of rate types
        ngrps = numel(EnsembleRates(iPrep,iStim).Ensemble);
        rateType = zeros(ngrps,1);
        for iS = 1:ngrps
            rateType(iS) = sign(EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).R) .* ...
                                (abs(EnsembleRates(iPrep,iStim).Ensemble(iS).Window(iW).R) >= Rlimit);
            Nneurons(iS) = EnsembleRates(iPrep,iStim).Ensemble(iS).N;
        end

        % proportions of each
        Prop(iPrep,iStim).Inc = sum(Nneurons(rateType == 1)) ./ EnsembleRates(iPrep,iStim).Ntotal;
        Prop(iPrep,iStim).NoChange = sum(Nneurons(rateType == 0)) ./ EnsembleRates(iPrep,iStim).Ntotal;
        Prop(iPrep,iStim).Dec = sum(Nneurons(rateType == -1)) ./ EnsembleRates(iPrep,iStim).Ntotal;

    end
end

% plot them
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); 
set(gca, 'ColorOrder', prep_Cmap,'NextPlot', 'replacechildren');
plot([1:3],reshape([Prop.Inc],nPreps,nStims),'.-','Linewidth',plotlinewidth); hold on

xlabel('Stimulation','FontSize',fontsize)
ylabel('Increasing neurons (%)','FontSize',fontsize)
axis([0.5 3.5 0 0.3])
% set(gca,'XTick',[10^-2 10^-1 1],'XTickLabel',[0.01 0.1 1])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

% print -depsc Fig1d_spectra
% exportfig(gcf,'Panel_Spiral_Colorbar','Color',color,'Format',format,'Resolution',dpi)
% 

% plot them
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); 
set(gca, 'ColorOrder', prep_Cmap,'NextPlot', 'replacechildren');
plot([1:3],reshape([Prop.Dec],nPreps,nStims),'.-','Linewidth',plotlinewidth,'MarkerSize',M); hold on

xlabel('Stimulation','FontSize',fontsize)
ylabel('Decreasing neurons (%)','FontSize',fontsize)
axis([0.5 3.5 0 1])
% set(gca,'XTick',[10^-2 10^-1 1],'XTickLabel',[0.01 0.1 1])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

%% plot set of eigenvalues

[nPreps,nStims] = size(AllEigs);

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); 
hold on;
for iPrep = 1:nPreps
   for iStim = 1:nStims
        plot(1:numel(AllEigs(iPrep,iStim).mUniqueEigMag),AllEigs(iPrep,iStim).mUniqueEigMag,'o-',...
            'Color',prep_Cmap(iPrep,:),'MarkerSize',M,'MarkerFaceColor',prep_Cmap(iPrep,:),'Linewidth',plotlinewidth);
        line([0.5 4.5],[0 0],'Color',[0 0 0],'Linewidth',axlinewidth)
   end
end
set(gca,'XLim',[0.2 4.5],'XTick',1:4)
xlabel('Dimensions with magnitude changes','FontSize',fontsize)
ylabel('Mean real eigenvalue','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_OtherEigenvalues
% exportfig(gcf,'Panel_Spiral_Colorbar','Color',color,'Format',format,'Resolution',dpi)


%% plot proportion vs size of (max) positive real eigenvalue
ctr = 1;
for iPrep = 1:nPreps
   for iStim = 1:nStims
       EgMax(ctr) = max(AllEigs(iPrep,iStim).mUniqueEigMag);
       EgMin(ctr) = min(AllEigs(iPrep,iStim).mUniqueEigMag);
       PInc(ctr) = Prop(iPrep,iStim).Inc;
       PDec(ctr) = Prop(iPrep,iStim).Dec;
       ctr = ctr + 1;
   end
end

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); 
hold on;
plot(PInc,EgMax,'.')
xlabel('Increasing rate (%)','FontSize',fontsize)
ylabel('Max. magnitude of eigenvalue','FontSize',fontsize)

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); 
hold on;
plot(PDec,EgMin,'.')
xlabel('Decreasing rate (%)','FontSize',fontsize)
ylabel('Min. magnitude of eigenvalue','FontSize',fontsize)




%% plot proportion of complex eigenvalues in dominant position: quantify consistency of rotations

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); 
set(gca, 'ColorOrder', prep_Cmap,'NextPlot', 'replacechildren');
plot([1:3],100*reshape([AllEigs.Pcomplex12],nPreps,nStims),'.-','Linewidth',plotlinewidth,'MarkerSize',M); hold on

% for iPrep = 1:nPreps
% plot([1:3],100*[AllEigs(iPrep,:).Pcomplex1],'.-','Color',prep_Cmap(iPrep,:),'Linewidth',plotlinewidth,'MarkerSize',M); hold on
% 
% end

xlabel('Stimulation','FontSize',fontsize)
ylabel('Time rotating (%)','FontSize',fontsize)
axis([0.5 3.5 0 100])
% set(gca,'XTick',[10^-2 10^-1 1],'XTickLabel',[0.01 0.1 1])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc ExtFig_TimeRotating
