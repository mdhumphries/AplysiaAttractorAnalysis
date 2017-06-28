%% figure panels with exciting spiral(s)

clear all; close all;
if ispc
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
    spikepath = 'C:\Users\mqbssmhg.DS\Dropbox\SpikeData\Angela Bruno Aplysia Dynamics Data\datasets\deconcatenated_spks\experimental\';
end

w = 401;  % smoothing window

run figure_properties
M = 4; % need bigger markers here 
Mclr = [0.7 0.7 0.7];  % stimulation onset markers

MinClr = 200;

% which prep?
iPrep = 10; % prep #6 and prep #10

% raster
exampleStim = 1;            
offset = 2; % additional colormap entries to make sure the lightest ones are omitted
Rateclrs = {'Reds','Greys','Blues'};

iW = 3; % 20 s windows
Rlimit = 0.2; % to be considered changing...

load([filepath 'Classify attractors/RecurrenceManifoldStats_CommonAxes_AllPreps.mat'])

%% set up axes
% prep 6
if iPrep == 6
    xmin = -15; xmax = 19; ymin = -23; ymax =17; zmin = 0; zmax = 36;
    best_view = [-128 22];

elseif iPrep == 10
    % prep 10
    xmin = 0; xmax = 32; ymin = -15; ymax =6; zmin = -9; zmax = 7;
    best_view = [-164 58];
end

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 6 6]);
% draw axes lines at the back
axis([xmin xmax ymin ymax zmin zmax])
view(best_view)  
xlabel('x')
ylabel('y')

line([xmax,xmax],[ymin,ymin],[zmin,zmax],'Linewidth',axlinewidth,'Color',[0 0 0]); hold on % PC3
line([xmin,xmax],[ymin,ymin],[zmin,zmin],'Linewidth',axlinewidth,'Color',[0 0 0]);  % PC 1  
line([xmax,xmax],[ymin,ymax],[zmin,zmin],'Linewidth',axlinewidth,'Color',[0 0 0]);  % PC2

text(xmin-3,ymin+5,zmin,'PC1','Fontsize',fontsize)
text(xmax+5,ymax+5,zmin,'PC2','Fontsize',fontsize)
text(xmax,ymin-1,zmax+1,'PC3','Fontsize',fontsize)

axis off


%% now plot data
for iStim = 1:size(Viz,2)
    n = numel(Viz(iPrep,iStim).CommonProj(:,1));
    nDims = size(Viz(iPrep,iStim).CommonProj,2);
    smth = zeros(n,nDims);
    for iD = 1:nDims
        smth(:,iD) = movingAverage(Viz(iPrep,iStim).CommonProj(:,iD),w);
    end
    % scaling colours
    % blue to red
    cd = uint8([linspace(0,255,n); ...      % red in [0,255]
            zeros(1,n)+50; ...          % green
            linspace(255,0,n);...          % blue
            zeros(1,n)+255]);            % alpha
% light to dark        
%     cd = uint8([linspace(MinClr,0,n); ...      % red in [0,255]
%             linspace(MinClr,0,n); ...          % green
%             zeros(1,n)+255;...          % blue
%             zeros(1,n)+255]);            % alpha
        
    h1 = plot3(smth(:,1),smth(:,2),smth(:,3),'Linewidth',errorlinewidth); hold on
    plot3(smth(1,1),smth(1,2),smth(1,3),'o','MarkerSize',M,'MarkerFaceColor',Mclr,'MarkerEdgeColor',Mclr)
    
%     h1 = plot3(smth(:,2),smth(:,3),smth(:,1),'Linewidth',errorlinewidth); hold on
%     plot3(smth(1,2),smth(1,3),smth(1,1),'o','MarkerSize',M,'MarkerFaceColor',Mclr,'MarkerEdgeColor',Mclr)

    grid on 
    drawnow
    set(h1.Edge, 'ColorBinding','interpolated', 'ColorData',cd)

end
% print -depsc Fig3b_examplevolume
exportfig(gcf,['Panel_Spiral_Prep' num2str(iPrep)],'Color',color,'Format',format,'Resolution',dpi)

%% 2D plots
for i = 1:3
    for j = i+1:3
        figure
        h1 = plot(smth(:,i),smth(:,j),'Linewidth',errorlinewidth); hold on
        plot(smth(1,i),smth(1,j),'o','MarkerSize',M,'MarkerFaceColor',Mclr,'MarkerEdgeColor',Mclr)
        xlabel(['PC ' num2str(i)]);
        ylabel(['PC ' num2str(j)])
    end
end

%% and create a colorbar
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 5]);
h1 = plot(ones(1,n),1:n,'Linewidth',5); % thick line
drawnow
set(h1.Edge, 'ColorBinding','interpolated', 'ColorData',cd);  % same colors
axis([0 2 0 n+100])
axis off
dTick = 0.05;
ticks = 0:1000:9000;
line([ones(1,numel(ticks)); ones(1,numel(ticks))+dTick],[ticks; ticks],'Color',[0 0 0],'Linewidth',0.5)
for iT = 1:numel(ticks)
    text(1+dTick+0.025,ticks(iT),num2str(ticks(iT)/100),'Fontsize',fontsize);
end


% print -depsc Panel_Spiral_Colorbar
exportfig(gcf,'Panel_Spiral_Colorbar','Color',color,'Format',format,'Resolution',dpi)

%% plot clustered ensemble for one stim
load([filepath '/Ensembles/Rates_StaticEnsembles'])

% load spike-trains
strStim = ['da0' num2str(exampleStim)];
load([filepath strStim '_DataProperties_FunctionAndWindowSize'],'Data','FileTable','DataTable')
load([spikepath FileTable{iPrep}]);
startts = 0; % floor(DataTable(iR,2));
% endts = floor(DataTable(ixP,3));

tStim = Data(iPrep).stimbins(1);
tEnd = floor(DataTable(iPrep,3));

% all spike-train IDs
allIDs = unique(spks(:,1));

% load clusters
load([filepath '/Ensembles/' strStim '_StaticEnsembles.mat'],'StaticEnsemblesAll')
n = numel(StaticEnsemblesAll(iPrep).IDs);    
spkIDs = StaticEnsemblesAll(iPrep).grps;  % group IDs of all neurons

% sort into change of rate types
ngrps = StaticEnsemblesAll(iPrep).ngrps;
rateType = zeros(ngrps,1);
for iS = 1:ngrps
    rateType(iS) = sign(EnsembleRates(iPrep,exampleStim).Ensemble(iS).Window(iW).R) .* ...
                        (abs(EnsembleRates(iPrep,exampleStim).Ensemble(iS).Window(iW).R) >= Rlimit);
    
end

% re-order into groups
grpctr = 1; Cctr=1;
cmapGrp = []; 

for iT = [1,0,-1]  % order from bottom to top of raster 
    % find groups that are this type
    ixR = find(rateType == iT);
    
    % renumber 
    for iG = ixR'  % only works if is a row vector...
        spkIDs(StaticEnsemblesAll(iPrep).grps == iG) = grpctr;
        grpctr = grpctr+1;
    end
    
    % built type colormap
    Cmap2 = brewermap(2+offset,Rateclrs{Cctr});
    thisCmap = repmat(Cmap2(end-1:end,:),ceil(numel(ixR)/2),1);
    thisCmap = thisCmap(1:numel(ixR),:);
    cmapGrp = [cmapGrp; thisCmap];
    Cctr = Cctr +1;  
end

spkIDs = [allIDs spkIDs];  % 2 column vector: neuron ID and grou

% [~,ixG] = sort(rateType,'descend');
% NG = zeros(size(StaticEnsemblesAll(iPrep).grps));
% NG(:,1) = 1:numel(allIDs);
% tempG = StaticEnsemblesAll(iPrep).grps; 
% for j = 1:ngrps NG(tempG == ixG(j),2) = j; end  % remap group membership indices


% plot them
figure
[H,GS,I,O] = plot_clusters(spks,spkIDs,ngrps,[0 max(spks(:,2))+1],'3',[],cmapGrp,M);
set(H,'Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 5.5 7]); 
title('');
line([110,120],[-1,-1],'Linewidth',plotlinewidth,'Color',[0 0 0],'Clipping','off')
text(115,-3,'10 s','FontSize',fontsize)
axis([tStim, max(get(gca,'XLim')), 0, n+1])
axis off

exportfig(gcf,['Spiral_Prep' num2str(iPrep) '_' strStim '_RasterPlot'],'Color',color,'Format',format,'Resolution',dpi)





