% Extended Data Figure: spontaneous peturbations:
% Example of
% (1) return to attractor:  (10,3)?
% (2) both different attractor examples (3,2); (7,2) 
% 
% Recurrence, and raster?
% And try 2D PCA, to see it???


clear all; close all

if ismac
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';
    % spikepath = '/Users/mqbssmhg/Dropbox/paper2/datasets/deconcatenated_spks/experimental/';
    % FIX
else
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
    spikepath = 'C:\Users\mqbssmhg.DS\Dropbox\SpikeData\Angela Bruno Aplysia Dynamics Data\datasets\deconcatenated_spks\experimental\';
end


run figure_properties
ixT = 3; % max threshold

prtblinewidth = 1;

% need to get times of:
% (1) perturbation crossing
% (2) start and end of period (black bar on "jump" example)
load([filepath 'Classify attractors/JumpAnalysis2'])  % get data on jumps
load([filepath 'Classify attractors/RecurrenceStats_FilterPts'])  % use recurrence point statistics based on using more robust model fits - min 1s and contiguous points

% keyboard

examples = {[4,3],[7,1],[7,2],[2,1],[1,3]};  % [animal, program];  % fall off attractor; different attractor & fall off; different attractor; same attractor?

% USED in Supplemental Figure: 
% a: [4,3]: return to spontaneous
% b: [1,3]: return to different manifold
% c: [2,1]: return to same manifold

%%  plot and export all examples
for iE = 1:numel(examples)
    ixPrep = examples{iE}(1);
    ixProg = examples{iE}(2);
    
    strStim = ['da0' num2str(ixProg)];
    
    % recurrence plot
    load([filepath '/Classify attractors/' strStim '_RecurrencePlotData.mat'])
    load([filepath '/Classify attractors/' strStim '_RecurrencePointsData_FilterPts.mat'])  

    load([filepath strStim '_DataProperties_FunctionAndWindowSize'],'Data','FileTable','DataTable')
    tStim = Data(ixPrep).stimbins(1);
    bins = Data(ixPrep).stimbins - tStim;
    maxBins = bins(end);
    tEnd = floor(DataTable(ixPrep,3));
%     tEnd = maxBins-pars.maxWindow; % end of checked recurrence time
%     bins(bins >= tEnd) = [];
        

    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5.5 5.5]);
    imagesc(bins,bins,RecurData(ixPrep).RecurPlot(ixT).Rp); hold on
    colormap(gray)
    
    if ~isempty(Jump(ixPrep,ixProg).Transition)
        % add divergent+return lines
        tJumpSt = Data(ixPrep).stimbins(RcrPtData(ixPrep).winmids(Jump(ixPrep,ixProg).Transition(1).ixWinJumpStart)) - tStim;
        tJumpEnd = Data(ixPrep).stimbins(RcrPtData(ixPrep).winmids(Jump(ixPrep,ixProg).Transition(1).ixWinJumpEnd)) - tStim;

        tPerturb = RcrStats(ixPrep,ixProg).tsPerturb(1) - tStim;

%         line([tJumpSt bins(end)],[tJumpSt tJumpSt],'Linewidth',axlinewidth,'Color',jumpcolor);
%         line([tJumpSt tJumpSt],[tJumpSt bins(end)],'Linewidth',axlinewidth,'Color',jumpcolor);
% 
%         line([tJumpEnd bins(end)],[tJumpEnd tJumpEnd],'Linewidth',axlinewidth,'Color',jumpcolor);
%         line([tJumpEnd tJumpEnd],[tJumpEnd bins(end)],'Linewidth',axlinewidth,'Color',jumpcolor);

        line([tPerturb bins(end)],[tPerturb tPerturb],'LineStyle','-','Linewidth',prtblinewidth,'Color',jumpcolor);
        line([tPerturb tPerturb],[tPerturb bins(end)],'LineStyle','-','Linewidth',prtblinewidth,'Color',jumpcolor);
    end
    
    % if has "end" divergence, add that too
    blnEnd = any(Summary.JumpList(:,1) == ixPrep & Summary.JumpList(:,2) == ixProg & Summary.JumpList(:,3)==0);
    if blnEnd
        tOff = RcrStats(ixPrep,ixProg).tsOff(1) - tStim;
        line([tOff bins(end)],[tOff tOff],'Linewidth',prtblinewidth,'Color',endcolor);
        line([tOff tOff],[tOff bins(end)],'Linewidth',prtblinewidth,'Color',endcolor);

    end
    % keyboard

    
    set(gca,'XTick',[0 30 60 90],'YTick',[0 30 60 90])

    xlabel('Time post-stimulation (s)','Fontsize',fontsize)
    ylabel('Time post-stimulation (s)','Fontsize',fontsize)

    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
    
    % print(['ExtFig_Prep' num2str(ixPrep) '_' strStim '_RecurPlot'],'-depsc') 
    exportfig(gcf,['ExtFig_Prep' num2str(ixPrep) '_' strStim '_RecurPlot'],'Color',color,'Format',format,'Resolution',dpi)


    %% clustered raster, with perturbation time
    % load spike-trains
    load([spikepath FileTable{ixPrep}]);
    startts = 0; % floor(DataTable(iR,2));
    % endts = floor(DataTable(ixP,3));

    % all spike-train IDs
    allIDs = unique(spks(:,1));

    % load clusters
    load([filepath '/Ensembles/' strStim '_StaticEnsembles.mat'],'StaticEnsemblesAll')
    n = numel(StaticEnsemblesAll(ixPrep).IDs);    
    Sxy = StaticEnsemblesAll(ixPrep).Cxy; Sxy(Sxy < 0) = 0;
    Sxy(eye(n)==1) = 0;
    
    [NG,Sg,Sin,Sout] = sortbysimilarity([StaticEnsemblesAll(ixPrep).IDs StaticEnsemblesAll(ixPrep).grps],Sxy);
    
    % plot them
    figure
    [H,GS,I,O] = plot_clusters(spks,NG,StaticEnsemblesAll(ixPrep).ngrps,[0 tEnd],'B3',[],M);
    set(H,'Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 5.5 7]); 
    title('');
    line([110,120],[-1,-1],'Linewidth',plotlinewidth,'Color',[0 0 0],'Clipping','off')
    text(115,-3,'10 s','FontSize',fontsize)
    axis([tStim, max(get(gca,'XLim')), 0, n+1])
    axis off

    % event lines
    % line([30 30],[0 numel(StaticEnsemblesAll(ixPrep).IDs)+1],'Color',stimcolor)
    if ~isempty(Jump(ixPrep,ixProg).Transition)
        tJumpSt = Data(ixPrep).stimbins(RcrPtData(ixPrep).winmids(Jump(ixPrep,ixProg).Transition(1).ixWinJumpStart));
        tJumpEnd = Data(ixPrep).stimbins(RcrPtData(ixPrep).winmids(Jump(ixPrep,ixProg).Transition(1).ixWinJumpEnd));
        tPerturb = RcrStats(ixPrep,ixProg).tsPerturb(1);

%         line([tJumpSt tJumpSt],[0 n+1],'Color',jumpcolor)
%         line([tJumpEnd tJumpEnd],[0 n+1],'Color',jumpcolor)
        line([tPerturb tPerturb],[0 n+1],'LineStyle','-','Linewidth',prtblinewidth,'Color',jumpcolor);
    end
  
    if blnEnd
        tOff = RcrStats(ixPrep,ixProg).tsOff(1);
        line([tOff tOff],[0 n+1],'Linewidth',prtblinewidth,'Color',endcolor);
    end

    %ptch = fill([30 32.5 32.5 30],[6 6 0 0],stimcolor);
    %set(ptch,'EdgeColor',stimcolor)

    % print(['ExtFig_Prep' num2str(ixPrep) '_' strStim '_RasterPlot'],'-depsc') 
    exportfig(gcf,['ExtFig_Prep' num2str(ixPrep) '_' strStim '_RasterPlot'],'Color',color,'Format',format,'Resolution',dpi)
    
    %% PCA projection, with start/end times
end
