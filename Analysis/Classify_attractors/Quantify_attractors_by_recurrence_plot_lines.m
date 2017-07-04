%% quantifying attractor from recurrence plot properties
% Main issue: which is the correct histogram to use?
% The definition in Marwan et al 2007 (Eq. 45) is the number of unique
% lines: i.e. only counting the length of each clearly separate line
% The version in Webber & Zbilut 1994 suggests instead a cumulative
% histogram


clear all; close all
addpath ../../Functions/

% pick which threshold to look at
RPpars.ixT = 3; % threshold of 10%
RPpars.lminT = 5; % minimum length of line in seconds to be longer than accidental recapitulation of prior dynamics
RPpars.vminT = 1; % minimum length of line in seconds to be longer than accidental recapitulation of prior dynamics

load da01_RecurrencePlotData
nfiles = numel(RecurData);

fnames = {'da01','da02','da03'};

NoPerColor = [0.5 0.5 0.5];
PerColor = [1 0 0];

%% visualise...
for iR = 1:nfiles
    iR
   for iF = 1:numel(fnames) 
       iF
       % load 
       load([fnames{iF} '_RecurrencePlotData'],'RecurData');
       load(['../' fnames{iF} '_DataProperties_FunctionAndWindowSize'],'Data','FileTable')
       
       % dimensions
       RPstats(iR,iF).maxD = RecurData(iR).maxD;
       
       % set minimum in time-steps
       RPpars.lmin = round(RPpars.lminT ./ Data(iR).GaussQt); % 
       RPpars.vmin = round(RPpars.vminT ./ Data(iR).GaussQt); % 

       % keyboard
       
%        %% dimensionality of system - not large enough set of thresholds
%        % to do this?
%        % get all P(l)
%        
%        figure; hold on
%        cmapD2 = cbrewer('qual','Set2',3); ctr = 1;
%        
%        setL = 0.1:0.1:40; % effectively the "embedding dimension"
%        ixL = round(setL ./ Data(iR).GaussQt);
% 
%        for iT = 1:numel(RecurData(iR).Theta)
%             L = RecurData(iR).RecurPlot(iT).dlengths .* Data(iR).GaussQt; % line lengths in seconds 
%             Lbins = Data(iR).GaussQt:Data(iR).GaussQt:max(setL);
%             Lhist = hist(L,Lbins);
%             Dim(iT).Phist = Lhist ./ sum(Lhist);
%             for iL = 1:numel(setL)
%                 Dim(iT).C2(iL) = sum(Dim(iT).Phist(1:ixL(iL)));
%             end
%             Dim(iT).v = Dim(iT).C2
%             
%             subplot(211),plot(log2(setL),log2(Dim(iT).C2),'.','Color',cmapD2(iT,:)); hold on
%             subplot(212),
%             
%             plot(setL,log2(Dim(iT).C2)/log2(RecurData(iR).Theta(iT)),'.','Color',cmapD2(iT,:)); hold on
%             
%             ylabel('D2 estimate'); xlabel('l')
% 
%        end 
%        
%        % compute D2 estimator
%        figure; hold on
%        for iT1 = 1:numel(RecurData(iR).Theta)
%            for iT2 = iT1+1:numel(RecurData(iR).Theta)
%                num = log(RecurData(iR).Theta(iT1) / RecurData(iR).Theta(iT2)); 
%                 for iL = 1:numel(setL)
%                     p1 = sum(Dim(iT1).Phist(1:ixL(iL)));
%                     p2 = sum(Dim(iT2).Phist(1:ixL(iL)));
%                     D2(iT1,iT2,iL) =  log(p1/p2) / num;  % Marwan et al
%                     % 2007 Eq 79 - not clear what it means, as it is an
%                     % estimator at some length L
%                 end  
%                 plot(setL,squeeze(D2(iT1,iT2,:)),'.','Color',cmapD2(ctr,:))
%                 ctr = ctr+1;
%            end
%        end
%        legend('1st vs 5th','1st vs 10th','5th vs 10th')
        
       %% estimate correlation dimension using recurrence network... (Ch4 in
       % book & Donner et al 2011)
%        for iT = 1:numel(RecurData(iR).Theta)
%            A = RecurData(iR).RecurPlot(iT).Rp;  
%            NetTrans = clusttriang(real(A));
%            RPstats(iR,iF).D2(iT) = log(NetTrans) / log(3/4);
%        end
%        
       %% diagonal lines - for largest threshold
       % do histogram - in steps
       L = RecurData(iR).RecurPlot(RPpars.ixT).dlengths; % .* Data(iR).GaussQt; % line lengths in seconds 
 
       % max L; DIVERGENCE = 1/Lmax 
       RPstats(iR,iF).Lmax = max(L) .* Data(iR).GaussQt;

       RPstats(iR,iF).Lbins = 1:1:max(L);
       RPstats(iR,iF).Lhist = hist(L,RPstats(iR,iF).Lbins);  % histogram of unique occurrences
       
       % now make histogram of all occurrences
       RPstats(iR,iF).Hl = zeros(1,max(L));
       for iV = 2:max(L)
            if RPstats(iR,iF).Lhist(iV)
                cntrb = iV:-1:2;  % for each line of this length, how many smaller lines of length 1 to length l-1 does it add?
                RPstats(iR,iF).Hl(1:iV-1) = RPstats(iR,iF).Hl(1:iV-1) + RPstats(iR,iF).Lhist(iV) .* cntrb;  % multiplied by the number of lines
            end
       end
       % keyboard
       
       % combine to get final histogram
       RPstats(iR,iF).Hl = RPstats(iR,iF).Hl + RPstats(iR,iF).Lhist;
       
       % now estimate correlation entropy (K2)
       % RPstats(iR,iF).cdfL = fliplr(cumsum(fliplr(RPstats(iR,iF).Hl))); % number of lines greater to or equal to L
%        figure
%        loglog(RPstats(iR,iF).Lbins,RPstats(iR,iF).Hl,'.')
       
        

       
       % keyboard
       
%        figure
%        subplot(211),
%        bar(RPstats(iR,iF).Lbins,RPstats(iR,iF).Lhist); hold on
%        line([RPpars.lmin RPpars.lmin],[0 max(RPstats(iR,iF).Lhist)+10],'Color','r')
%        title(['Recording' num2str(iR) ' Prep ' num2str(iF)])
%        xlabel('Diagonal line length (s)') 
       
       % determinism
       RPstats(iR,iF).DetHl = sum(RPstats(iR,iF).Lbins(RPstats(iR,iF).Lbins >= RPpars.lmin) .* RPstats(iR,iF).Hl(RPstats(iR,iF).Lbins >= RPpars.lmin)) / sum(RPstats(iR,iF).Hl.* RPstats(iR,iF).Lbins);
       RPstats(iR,iF).DetL = sum(RPstats(iR,iF).Lbins(RPstats(iR,iF).Lbins >= RPpars.lmin) .* RPstats(iR,iF).Lhist(RPstats(iR,iF).Lbins >= RPpars.lmin)) / sum(RPstats(iR,iF).Lhist.* RPstats(iR,iF).Lbins);
              
       % mean prediction time 
       % RPstats(iR,iF).mPredTime = sum(RPstats(iR,iF).Lbins(RPstats(iR,iF).Lbins >= RPpars.lmin) .* RPstats(iR,iF).Hl(RPstats(iR,iF).Lbins >= RPpars.lmin)) / sum(RPstats(iR,iF).Hl(RPstats(iR,iF).Lbins >= RPpars.lmin)); 
        
       RPstats(iR,iF).mPredTimeL = mean(L(L >= RPpars.lmin)) .* Data(iR).GaussQt;  % mean duration of all diagonal structures (unique set of lines)
       RPstats(iR,iF).mPredTimeHl = Data(iR).GaussQt .* sum(RPstats(iR,iF).Lbins(RPstats(iR,iF).Lbins >= RPpars.lmin) .* RPstats(iR,iF).Hl(RPstats(iR,iF).Lbins >= RPpars.lmin)) / sum(RPstats(iR,iF).Hl(RPstats(iR,iF).Lbins >= RPpars.lmin));       

       % entropy (on beyond tangenital motion) in bits
       Hmin = RPstats(iR,iF).Lhist(RPstats(iR,iF).Lbins >= RPpars.lmin);
       Pmin = Hmin ./ sum(Hmin);
       RPstats(iR,iF).EntL = -sum(Pmin(Pmin > 0) .* log2(Pmin(Pmin > 0))); %P=0 contributes 0
      
       Hmin = RPstats(iR,iF).Hl(RPstats(iR,iF).Lbins >= RPpars.lmin);
       Pmin = Hmin ./ sum(Hmin);
       RPstats(iR,iF).EntHl = -sum(Pmin(Pmin > 0) .* log2(Pmin(Pmin > 0))); %P=0 contributes 0
       
       %% vertical lines
       V = RecurData(iR).RecurPlot(RPpars.ixT).Vlengths; % line lengths in seconds 
       
       % max V;
       RPstats(iR,iF).Vmax = max(V)  .* Data(iR).GaussQt;

       RPstats(iR,iF).Vbins = 1:1:max(V);
       RPstats(iR,iF).Vhist = hist(V,RPstats(iR,iF).Vbins);  % uniqur vertical line lengths
      
       % now make histogram of all occurrences
       RPstats(iR,iF).Hv = zeros(1,max(V));
       for iV = 2:max(V)
            if RPstats(iR,iF).Vhist(iV)
                cntrb = iV:-1:2;  % for each line of this length, how many smaller lines of length 1 to length l-1 does it add?
                RPstats(iR,iF).Hv(1:iV-1) = RPstats(iR,iF).Hv(1:iV-1) + RPstats(iR,iF).Vhist(iV) .* cntrb;  % multiplied by the number of lines
            end
       end
       % combine to get final histogram
       RPstats(iR,iF).Hv = RPstats(iR,iF).Hv + RPstats(iR,iF).Vhist;
       
%        subplot(212),
%        bar(RPstats(iR,iF).Vbins,RPstats(iR,iF).Vhist); hold on
%        line([RPpars.vmin RPpars.vmin],[0 max(RPstats(iR,iF).Vhist)+10],'Color','r')
%        xlabel('Vertical line length (s)') 
       
       % laminarity
       RPstats(iR,iF).LamHv = sum(RPstats(iR,iF).Vbins(RPstats(iR,iF).Vbins >= RPpars.vmin) .* RPstats(iR,iF).Hv(RPstats(iR,iF).Vbins >= RPpars.vmin)) / sum(RPstats(iR,iF).Vbins.* RPstats(iR,iF).Hv);
       RPstats(iR,iF).LamV = sum(RPstats(iR,iF).Vbins(RPstats(iR,iF).Vbins >= RPpars.vmin) .* RPstats(iR,iF).Vhist(RPstats(iR,iF).Vbins >= RPpars.vmin)) / sum(RPstats(iR,iF).Vbins.* RPstats(iR,iF).Hv);
  
       % trapping time
       RPstats(iR,iF).mTrapTimeV = mean(V(V >= RPpars.vmin)) .* Data(iR).GaussQt;
       RPstats(iR,iF).mTrapTimeHv = Data(iR).GaussQt .* sum(RPstats(iR,iF).Vbins(RPstats(iR,iF).Vbins >= RPpars.vmin) .* RPstats(iR,iF).Hv(RPstats(iR,iF).Vbins >= RPpars.vmin)) / sum(RPstats(iR,iF).Hv(RPstats(iR,iF).Vbins >= RPpars.vmin));       

   end
end

%% plots

cmap = cbrewer('qual','Paired',nfiles); 
% summary plot: 
% (1) color code by prep to show consistency within prep
% (2) link within prep to draw eye to groups
% (3) markersize according to stim sequence: no patterns
% (4) markersize according to "sensitised"

Msize = [5 7 9];  % stim order
% Msize = [5 15 10];  % does sensitised prep show clear effect?

hT = figure; hold on
hE = figure; hold on
hDV = figure; hold on 
h3 = figure; 

join ='o';
sym = 'o';
for iR = 1:nfiles
    % ID P10 preps using squares
%     if any(iR == P10preps')
%         sym = 's'; 
%     else
%         sym = 'o';
%     end
    
    % compare diagonal and vertical trapping time
    figure(hDV)
%     subplot(211),plot([RPstats(iR,:).mPredTime],[RPstats(iR,:).mTrapTime],join,...
%         'Color',cmap(iR,:),'MarkerSize',2,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5); hold on
    plot([RPstats(iR,:).Lmax],[RPstats(iR,:).Vmax],join,...   
        'Color',cmap(iR,:),'MarkerSize',2,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5); hold on
    
    % compare max time on attractor with determinism of movement
    figure(hT)
    % plot([RPstats(iR,:).mPredTime],[RPstats(iR,:).Det],join,...
    plot([RPstats(iR,:).Lmax],[RPstats(iR,:).DetL],join,...   
        'Color',cmap(iR,:),'MarkerSize',2,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
    
    figure(hE) 
    % plot([RPstats(iR,:).mPredTime],[RPstats(iR,:).Ent],join,...
    plot([RPstats(iR,:).Lmax],[RPstats(iR,:).EntL],join,...
        'Color',cmap(iR,:),'MarkerSize',2,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
    
    figure(h3)
    plot3([RPstats(iR,:).Lmax],[RPstats(iR,:).DetL],[RPstats(iR,:).EntL],join,...
        'Color',cmap(iR,:),'MarkerSize',8,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5); hold on
   
     % increasing MarkerSize for sequence
     for iF = 1:numel(fnames)
        figure(hDV)
%         subplot(211),plot(RPstats(iR,iF).mPredTime,RPstats(iR,iF).mTrapTime,sym,...
%             'Color',cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',cmap(iR,:)); 
%         if iF == 3
%             text([RPstats(iR,iF).mPredTime],[RPstats(iR,iF).mTrapTime],num2str(iR),'Fontsize',8,'Color',[1 1 1])
%         end

        % subplot(212)
        plot([RPstats(iR,iF).Lmax],[RPstats(iR,iF).Vmax],join,...   
            'Color',cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',cmap(iR,:)); 
        if iF == 3
            text([RPstats(iR,iF).Lmax],[RPstats(iR,iF).Vmax],num2str(iR),'Fontsize',8,'Color',[1 1 1])
        end
        
        figure(hT)
        % plot(RPstats(iR,iF).mPredTime,RPstats(iR,iF).Det,sym,...
        plot([RPstats(iR,iF).Lmax],[RPstats(iR,iF).DetL],sym,...     
            'Color',cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',cmap(iR,:));
        if iF == 3
            %text([RPstats(iR,iF).mPredTime],[RPstats(iR,iF).Det],num2str(iR),'Fontsize',8,'Color',[1 1 1])
            text([RPstats(iR,iF).Lmax],[RPstats(iR,iF).DetL],num2str(iR),'Fontsize',8,'Color',[1 1 1])
        end
        
        figure(hE)
        % plot(RPstats(iR,iF).mPredTime,RPstats(iR,iF).Ent,sym,...
        plot(RPstats(iR,iF).Lmax,RPstats(iR,iF).EntL,sym,...
            'Color',cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',cmap(iR,:));
        if iF == 3
            % text([RPstats(iR,iF).mPredTime],[RPstats(iR,iF).Ent],num2str(iR),'Fontsize',8,'Color',[1 1 1])
            text([RPstats(iR,iF).Lmax],[RPstats(iR,iF).EntL],num2str(iR),'Fontsize',8,'Color',[1 1 1])
        end
        
    end
end
figure(hDV)
% subplot(211), xlabel('Mean prediction time (s)'); ylabel('Mean trapping time (s)')
%subplot(212), 
xlabel('Max prediction time (s)'); ylabel('Max trapping time (s)')
exportPPTfig(gcf,'MaxPredTimevsMaxTrapTime',[10 15 7 7])

figure(hT)
% xlabel('Mean prediction time (s)')
xlabel('Max prediction time (s)')
ylabel('Determinism')
exportPPTfig(gcf,'MaxPredTimevsDET',[10 15 7 7])

figure(hE)
% xlabel('Mean prediction time (s)')
xlabel('Max prediction time (s)')
ylabel('Entropy (bits)')
exportPPTfig(gcf,'MaxPredTimevsEntropy',[10 15 7 7])

figure(h3)
xlabel('Max prediction time (s)')
ylabel('Determinism')
zlabel('Entropy (bits)')
grid on

ds = [RPstats(:,:).maxD];
HD = hist(ds,1:1:10);
figure
bar(1:1:10,HD,'FaceColor',[0 0 0],'LineWidth',0);
xlabel('Embedding dimensions')
ylabel('#Programs')
exportPPTfig(gcf,'NdimensionsVAF80',[10 15 5 5])


save RecurPlotStats RPstats RPpars
