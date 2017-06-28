%% script to assemble sub-panels for Figure X: single neuron contribution...
clear all; close all

if ispc
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\';
else
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/';
end

run figure_properties
fnames = {'da01','da02','da03'};

% M = 4; % need bigger markers here 

load([filepath '/Classify attractors/ParticipationChanges_SamePCs'])
load([filepath '/Classify attractors/RecurrenceManifoldStats_CommonAxes_AllPreps.mat'],'spars','Prep'); % get analysis parameters
% got to load the Common_Axes set here, as they have the "end program" ID
% points for the concatenated SDFs

npreps = size(Coeffs,2);
nprogs = 3;

dist_panel_sz = [4 4];

%% panel a: contribution to concatenated time-series
hbinsNorm = 0:0.05:1;

hHst = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 dist_panel_sz]); hold on
hECDF = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 dist_panel_sz]); hold on

for iPrep = 1:npreps
    PartNorm = Coeffs(iPrep).sumW ./ max(Coeffs(iPrep).sumW);
    hPartNorm = hist(PartNorm,hbinsNorm);
    pPartNorm = hPartNorm ./ sum(hPartNorm);
    mPartNorm(iPrep) = median(PartNorm);
    
    for j = iPrep+1:npreps
        PartNorm2 = Coeffs(j).sumW ./ max(Coeffs(j).sumW);
        [H,P(iPrep,j),ksstat(iPrep,j)] = kstest2(PartNorm,PartNorm2,'alpha',0.05);
    end
    
    figure(hECDF)
    [f,x] = ecdf(Coeffs(iPrep).sumW);
    stairs(x,f,'Color',prep_Cmap(iPrep,:))
    
    figure(hHst)
    % stairs(hbinsNorm,pPartNorm,'Color',prep_Cmap(iPrep,:))
    plot(hbinsNorm,pPartNorm,'Color',prep_Cmap(iPrep,:))
end

% line([20 100],[35 35],'Color',[0.5 0.5 0.5],'Linewidth',axlinewidth);
figure(hHst)
xlabel('Participation (% of max)','FontSize',fontsize)
ylabel('P(participation)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
% exportPPTfig(gcf,'Participation_Concat_Hist',[10 15 8 8])


figure(hECDF)
xlabel('Participation','FontSize',fontsize)
ylabel('P(participation) > X','FontSize',fontsize)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
%exportPPTfig(gcf,'Participation_Concat_ECDF',[10 15 8 8])

% print -depsc Fig3a_coalescence
%% Distributions of each program

hHst = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 dist_panel_sz]); hold on
hECDF = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 dist_panel_sz]); hold on

allHsts = zeros(npreps*nprogs,numel(hbinsNorm)); ctr = 1;
for iPrep = 1:npreps
    for iProg = 1:nprogs
        hPartNorm = hist(Coeffs(iPrep).normWMax(iProg,:),hbinsNorm);
        pPartNorm = hPartNorm ./ sum(hPartNorm);
        mPartNormAll(iPrep,iProg) = median(Coeffs(iPrep).normWMax(iProg,:));

        figure(hECDF)
        [f,x] = ecdf(Coeffs(iPrep).sumWProg(iProg,:));
        stairs(x,f,'Color',prep_Cmap(iPrep,:))

        figure(hHst)
        % stairs(hbinsNorm,pPartNorm,'Color',prep_Cmap(iPrep,:))
        plot(hbinsNorm,pPartNorm,'Color',prep_Cmap(iPrep,:))
        
        allHsts(ctr,:) = pPartNorm; ctr = ctr + 1;
    end
end

figure(hHst)
xlabel('Participation (% of max)','FontSize',fontsize)
ylabel('P(participation)','FontSize',fontsize)
set(gca,'Layer','top') % esnure that the axes are on top, so that plot objects don't overlap
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
%exportPPTfig(gcf,'Participation_Separate_Hist',[10 15 8 8])

figure(hECDF)
xlabel('Participation','FontSize',fontsize)
ylabel('P(participation) > X','FontSize',fontsize)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
%exportPPTfig(gcf,'Participation_Separate_ECDF',[10 15 8 8])

% plot all as greyscale, plus mean
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 dist_panel_sz]); hold on
plot(100*hbinsNorm,allHsts,'Color',[0.7 0.7 0.7],'Linewidth',errorlinewidth); hold on
plot(100*hbinsNorm,mean(allHsts),'Color',[0 0 0],'Linewidth',plotlinewidth);

xlabel('Participation (% of max)','FontSize',fontsize)
ylabel('P(participation)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
set(gca,'Layer','top') % esnure that the axes are on top, so that plot objects don't overlap

print -depsc FigPart_distributions

%% change between programs

% example program
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 6 4]); hold on
plot(Coeffs(3).HighMaxP,Coeffs(3).maxDMaxP,'ko')
axis([0 1 0 0.5])
xlabel('Highest Proportion of Maximum','FontSize',fontsize)
ylabel('Maximum change in proportion','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

%exportPPTfig(gcf,'Prep3_ChangeParticipation',[10 15 8 6])

% collapsed over all programs
allMaxP = []; allMaxDP = [];
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 4]); hold on
for iPrep = 1:npreps
    allMaxP = [allMaxP; Coeffs(iPrep).HighMaxP];
    allMaxDP = [allMaxDP;  Coeffs(iPrep).maxDMaxP];
    plot(100*Coeffs(iPrep).HighMaxP,100*Coeffs(iPrep).maxDMaxP,'ko','MarkerSize',M)
    % title(['Variable participation plot for all programs'])
    
end
axis([0 100 0 100])
xlabel('Maximum participation (%)','FontSize',fontsize)
ylabel('Maximum change in participation (%)','FontSize',fontsize)

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

print -depsc FigPart_PartChange

%exportPPTfig(gcf,'AllPrep_ChangeParticipation',[10 15 8 6])

nChangers = sum(allMaxP > 0.5 & allMaxDP > 0.25);

%% plot examples of strong and weak; and max and min change

ixExamplePart = [2,6];
ixExampleVar = [5,10]; 

nExamples = 2;  % how many of the top neurons from each prep

% just use grayscale to keep life simple....
progColor1 = [0.7 0.7 0.7;
            0.4 0.4 0.4;
            0 0 0];

% progColor1 = [0.7 0.4 0.4;
%             0.4 0.4 0.7;
%             0.4 0.7 0.4];

        
% progColor2 = [0.4 0.4 0.8;
%             0.2 0.2 0.8;
%             0 0 0.8];
% hStrong = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 6 15]); hold on
% hWeak = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 6 15]); hold on
% hMaxVar = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 6 15]); hold on
% hMinVar = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 6 15]); hold on

ctrPart = 1; ctrVar = 1;
for iPrep = 1:npreps
    conSDF = []; endProg = [];
    if any(iPrep == ixExamplePart) || any(iPrep == ixExampleVar)
        for iProg = 1:nprogs
           % load stuff and store across all 3 programs
           load([filepath fnames{iProg} '_DataProperties_FunctionAndWindowSize']) % to get density functions and basic data
           All(iProg).ixts = find(Data(iPrep).bins < spars.stimstart);
           All(iProg).ixStim = find(Data(iPrep).bins >= spars.stimstart);
           All(iProg).n = numel(Data(iPrep).rates);

           conSDF = [conSDF; SDF(iPrep).spkfcn(All(iProg).ixStim,:)];
           endProg(iProg) = Prep(iPrep).PIdx(iProg).ts(end);
        end

        % plot by global ranking
        mS = round(max(max(conSDF)))+1;
        endProg = [0 endProg];
        
        if any(iPrep == ixExamplePart)
            
            for iE = 1:nExamples
                maxPart = conSDF(:,Coeffs(iPrep).Rank(iE));
                % hOver = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2]);
                
                Y = round(max(maxPart) + 1);
                
                for iP = 2:nprogs+1
                    % figure(hOver); plot(maxPart(endProg(iP-1)+1:endProg(iP)),'Color',progColor1(iP,:)); hold on
                    
                    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3 1.5]);
                    plot(maxPart(endProg(iP-1)+1:endProg(iP)),'Color',progColor1(iP-1,:));
                    set(gca,'XTick',1:2000:10000,'XTickLabel',(0:2000:10000) .* 0.01);
                    set(gca,'YLim',[0 Y]);
                    xlabel('Time (s)','Fontsize',fontsize); ylabel('Rate (Hz)','Fontsize',fontsize)
                    set(gca,'FontName','Helvetica','FontSize',fontsize);
                    set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
                    set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap                 
                    print(['Strong_Participators_Prep' num2str(iPrep) '_Example' num2str(iE) '_Prog' num2str(iP-1)],'-depsc') 

                end
                
%                 % figure processing for overlap figure
%                 figure(hOver)
%                 set(gca,'XTick',1:2000:10000,'XTickLabel',(0:2000:10000) .* 0.01);
%                 axis tight
%                 xlabel('Time (s)','Fontsize',fontsize); ylabel('Rate (Hz)','Fontsize',fontsize)
%                 set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap
%                 set(gca,'FontName','Helvetica','FontSize',fontsize);
%                 set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
%                 print(['Strong_Participators_Prep' num2str(iPrep) '_Example' num2str(iE)],'-depsc') 
            end
                %exportPPTfig(gcf,['Strong_Participators_Prep_' num2str(iPrep)],[10 15 6 4])

%             figure(hWeak)
%             subplot(npreps,1,iPrep),plot(conSDF(:,Coeffs(iPrep).Rank(end-1:end))); hold on
%             line([endProg; endProg],[0 0 0; mS mS mS],'Color',[0 0 0],'Linewidth',axlinewidth)
%             set(gca,'FontName','Helvetica','FontSize',fontsize);
%             set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
        end
        
        if any(iPrep == ixExampleVar)
    
            [Varsrt,Varrnk] = sort(Coeffs(iPrep).maxDMaxP,'descend'); % from most to least variable
%             figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 6 4]); hold on
%             plot(conSDF(:,Varrnk(1:2))); hold on
%             axis tight
%             line([endProg; endProg],[0 0 0; mS mS mS],'Color',[0 0 0],'Linewidth',axlinewidth)
%             set(gca,'FontName','Helvetica','FontSize',fontsize);
%             set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
            for iE = 1:nExamples
                maxVar = conSDF(:,Varrnk(iE));
                hOver = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 2]);
                Y = round(max(maxVar) + 1);

                for iP = 2:nprogs+1
                    % figure(hOver); plot(maxVar(endProg(iP-1)+1:endProg(iP)),'Color',progColor1(iP,:)); hold on
                    
                    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3 1.5]);
                    plot(maxVar(endProg(iP-1)+1:endProg(iP)),'Color',progColor1(iP-1,:));
                    set(gca,'XTick',1:2000:10000,'XTickLabel',(0:2000:10000) .* 0.01);
                    set(gca,'YLim',[0 Y]);
                    xlabel('Time (s)','Fontsize',fontsize); ylabel('Rate (Hz)','Fontsize',fontsize)
                    set(gca,'FontName','Helvetica','FontSize',fontsize);
                    set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
                    set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap
                    print(['Variable_Participators_Prep' num2str(iPrep) '_Example' num2str(iE) '_Prog' num2str(iP-1)],'-depsc') 

                end
                
%                 % figure processing for overlap figure
%                 figure(hOver)
%                 set(gca,'XTick',1:2000:10000,'XTickLabel',(0:2000:10000) .* 0.01);
%                 axis tight
%                 xlabel('Time (s)','Fontsize',fontsize); ylabel('Rate (Hz)','Fontsize',fontsize)
%                 set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap
%                 set(gca,'FontName','Helvetica','FontSize',fontsize);
%                 set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
%                 print(['Variable_Participators_Prep' num2str(iPrep) '_Example' num2str(iE)],'-depsc') 
            end

%             figure(hMinVar)
%             subplot(npreps,1,iPrep),plot(conSDF(:,Varrnk(end-1:end))); hold on
%             line([endProg; endProg],[0 0 0; mS mS mS],'Color',[0 0 0],'Linewidth',axlinewidth)
%             set(gca,'FontName','Helvetica','FontSize',fontsize);
%             set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

            %if iR == 1 title('Minimum participation'); end
        end
    end
end

% figure(hWeak)
% exportPPTfig(gcf,'Weak_Participators',[10 15 8 20])
% 
% figure(hStrong)
% exportPPTfig(gcf,'Strong_Participators',[10 15 8 20])
% 
% figure(hMinVar)
% exportPPTfig(gcf,'Min_Changers',[10 15 8 20])
% 
% figure(hMaxVar)
% exportPPTfig(gcf,'Max_Changers',[10 15 8 20])


%% how many neurons genuinely change?

% example distribution etc
data = Coeffs(3).dMaxP(:,1:2);
bins = -0.4:0.025:0.5;
h = hist(data,bins);

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3 2]); hold on
bar(100*bins,h/sum(h),'FaceColor',[1 1 1],'LineWidth',axlinewidth);
axis([-30 50 0 0.2])
xlabel('Change in participation (%)','FontSize',fontsize)
ylabel('P(change)','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap

print -depsc FigPart_ExampleChangeDist

%exportPPTfig(gcf,'Prog3_example_change_histogram',[10 15 6 4])

% plot fit of model
x = bins(1):0.01:bins(end);
[fcdf,xcdf] = ecdf(data(:));
yc = cdf(Coeffs(3).Pars,x);

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3 2]); hold on
stairs(100*xcdf,fcdf,'Color',[0 0 0])
stairs(100*x,yc,'Color',[0.8 0.3 0.3]);
% label means and 3SD threshold
upper = 100*(Coeffs(3).Pars.mu + 3 * Coeffs(3).Pars.sigma);
lower = 100*(Coeffs(3).Pars.mu - 3 * Coeffs(3).Pars.sigma);
line([upper upper],[0 1],'Color',[0.7 0.7 0.7],'Linestyle','--','Linewidth',errorlinewidth)
line([lower lower],[0 1],'Color',[0.7 0.7 0.7],'Linestyle','--','Linewidth',errorlinewidth)

xlabel('Change in participation (%)','FontSize',fontsize)
ylabel('P(change)','FontSize',fontsize)
text(30,0.87,'Data','Color',[0 0 0],'FontSize',fontsize-1)
text(0,1.1,'Best-fit','Color',[0.8 0.3 0.3],'FontSize',fontsize-1)
axis([-30,45,0,1])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');
set(gca,'Layer','top') % esnure that the axes are on top, so the patch objects don't overlap

print -depsc FigPart_ExampleChangeModel

%exportPPTfig(gcf,'Prog3_example_fit_model',[10 15 6 4])

%% summarise number of changed neurons

varHist = zeros(npreps,1);
for iPrep = 1:npreps
    % proportion of neurons exceeding the change threshold at least once
    varHist(iPrep) = numel(Coeffs(iPrep).VarPart) ./ numel(Coeffs(iPrep).Rank); 
end

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 3 4]); hold on
barh(1:npreps,100*varHist,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);
axis([0 30 0.5 10.5])
ylabel('Animal','FontSize',fontsize)
xlabel('Variable neurons (%)','FontSize',fontsize)
% axis([0.5 10.5 0 0.3])
% xlabel('Animal','FontSize',fontsize)
% ylabel('Proportion variable neurons','FontSize',fontsize)
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'TickDir','out','LineWidth',axlinewidth,'Box','off');

%exportPPTfig(gcf,'Distribution_of_VariableNeurons',[10 15 6 4])
print -depsc FigPart_NeuronCount

