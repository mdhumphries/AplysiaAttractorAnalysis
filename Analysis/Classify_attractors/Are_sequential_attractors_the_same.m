%% are sequential attractors actually the same?
%% could there be more than one stable limit cycle/spiral in system? 
%
% NB this also checking across cases where there is clear transition in the
% system while it's running!

clear all; close all
addpath ../../Functions/

% analysis parameters
spars.stimstart = 30;  % stimulation time
spars.VmaxD = 3;  % number of dimensions to visualise
spars.minDens = 0; % min proportion of recurrence points with period
spars.ixT = 3; % which threshold to use
spars.minRecur = 90; % minimum percentage of recurrence points per time-window
spars.manifoldT = 5; % distance to be on same manifold
spars.divPeriod = 0;  % radius on checking min distance as division of oscillation period (0 = no radius; 1 = period; 2 = half period etc)
spars.nShuffles = 100;  % 100!

fname_Pre = '';
load ../da01_DataProperties_FunctionAndWindowSize FileTable
nfiles = numel(FileTable);
load RecurrenceStats RcrStats
fnames = {'da01','da02','da03'};

ixExample = 1;

% uses Parallel Computing Toolbox do 100 shuffled repeats
feature('numCores');
if isempty(gcp('nocreate'))
    parpool('local',10);
end

%% compare attractors
for iR = 1:nfiles
    All = [];
   tic
   for iF = 1:numel(fnames) 
       
       % load stuff and store across all 3 programs
       load(['../' fname_Pre fnames{iF} '_DataProperties_FunctionAndWindowSize']) % to get PCA
       load(['../Ensembles/' fname_Pre fnames{iF} '_StaticEnsembles.mat'])  % to get Cxy
       load([fname_Pre fnames{iF} '_RecurrencePointsData_FilterPts']); % to get recurrence: density per window, and winmids etc
       All(iF).n = numel(Data(iR).rates); % number of neurons
       
       %  store for all comparisons
       All(iF).Data = Data(iR);
       All(iF).PCAdata = PCAdata(iR);
       All(iF).SDF = SDF(iR);
       All(iF).Cxy = StaticEnsemblesAll(iR).Cxy;
       All(iF).winDT = RcrPtData(iR).winDT; All(iF).winmids = RcrPtData(iR).winmids; All(iF).winstrt = RcrPtData(iR).winstrt;  
       All(iF).PrctWnone = 100*RcrPtData(iR).Orbit(spars.ixT).NoRecurrence.nWinPoints/RcrPtData(iR).winDT;  % density of any recurrence point per window
       All(iF).ixWpts = find(All(iF).PrctWnone <= 100 - spars.minRecur);  % find all windows with high recurrence density
       
       % find set of admissible time-points in each attractor
        ixbreaks = [0 find(diff(All(iF).ixWpts) > 1) numel(All(iF).ixWpts)]; % discontinuties
        ixAdPts = [];
        % get all unbroken sequences of recurrence points
        for iB = 1:numel(ixbreaks)-1
            strts = All(iF).winstrt(All(iF).ixWpts(ixbreaks(iB)+1)); % start of first window in this unbroken sequence
            ends = All(iF).winstrt(All(iF).ixWpts(ixbreaks(iB+1)))+All(iF).winDT; % end of last window in this unbroken sequence
            ixAdPts = [ixAdPts strts:ends];
        end
        All(iF).ixAdPts = ixAdPts;
       % keyboard
       
        % collect stats for comparison
        Prep(iR).Prog(iF).nPCsVAF80 = find(cumsum(PCAdata(iR).eigvalues./sum(PCAdata(iR).eigvalues)) >= 0.8,1,'first');   
        Prep(iR).Prog(iF).nPCsVAF95 = find(cumsum(PCAdata(iR).eigvalues./sum(PCAdata(iR).eigvalues)) >= 0.95,1,'first');   

       % do PCA projection unnormalised - 80% VAF
        All(iF).ixts = find(Data(iR).bins < spars.stimstart);
        All(iF).ixStim = find(Data(iR).bins >= spars.stimstart);
        All(iF).SponPCprojU = zeros(numel(All(iF).ixts),PCAdata(iR).nPCs);  % spon PCs on stim axes
        All(iF).StimPCprojU = zeros(numel(All(iF).ixStim),PCAdata(iR).nPCs);  % stim PCs on stim axes
        for iP = 1:Prep(iR).Prog(iF).nPCsVAF80 
            All(iF).SponPCprojU(:,iP) = sum(bsxfun(@times,SDF(iR).spkfcn(All(iF).ixts,:),PCAdata(iR).coeffs(:,iP)'),2); % unnormalised 
            All(iF).StimPCprojU(:,iP) = sum(bsxfun(@times,SDF(iR).spkfcn(All(iF).ixStim,:),PCAdata(iR).coeffs(:,iP)'),2); % unnormalised 
        end
        
        % project onto da01 axes - 80% VAF
        if iF > 1
            All(iF).SponOnda01 = zeros(numel(All(iF).ixts),All(1).PCAdata.nPCs);  % spon PCs on stim axes
            All(iF).StimOnda01 = zeros(numel(All(iF).ixStim),All(1).PCAdata.nPCs);  % stim PCs on stim axes
            ixR = randperm(All(iF).n);
            for iP = 1:Prep(iR).Prog(1).nPCsVAF80 % onto da01 axes
                All(iF).SponOnda01(:,iP) = sum(bsxfun(@times,SDF(iR).spkfcn(All(iF).ixts,:),All(1).PCAdata.coeffs(:,iP)'),2); % unnormalised 
                All(iF).StimOnda01(:,iP) = sum(bsxfun(@times,SDF(iR).spkfcn(All(iF).ixStim,:),All(1).PCAdata.coeffs(:,iP)'),2); % unnormalised 
                All(iF).ControlSponOnda01(:,iP) = sum(bsxfun(@times,SDF(iR).spkfcn(All(iF).ixts,ixR),All(1).PCAdata.coeffs(:,iP)'),2); % shuffle neuron indices, and project onto da01 axes
                All(iF).ControlStimOnda01(:,iP) = sum(bsxfun(@times,SDF(iR).spkfcn(All(iF).ixStim,ixR),All(1).PCAdata.coeffs(:,iP)'),2); % shuffle neuron indices, and project onto da01 axes
            end 
        else
            All(iF).SponOnda01 = All(iF).SponPCprojU;
            All(iF).StimOnda01 = All(iF).StimPCprojU;
            % control: randomly assigned SDFs
            ixR = randperm(All(iF).n);
            for iP = 1:Prep(iR).Prog(1).nPCsVAF80 % onto da01 axes
                All(iF).ControlSponOnda01(:,iP) = sum(bsxfun(@times,SDF(iR).spkfcn(All(iF).ixts,ixR),All(1).PCAdata.coeffs(:,iP)'),2); % unnormalised 
                All(iF).ControlStimOnda01(:,iP) = sum(bsxfun(@times,SDF(iR).spkfcn(All(iF).ixStim,ixR),All(1).PCAdata.coeffs(:,iP)'),2); % unnormalised 
            end 
        end
        
        % keyboard
   end
   maxD = min([Prep(iR).Prog(:).nPCsVAF80]);  % max comparable dimensions across all three programs
        
   % (1) plot PCA - separate; project onto same axes

    % keyboard
    % projmap = cbrewer('qual','Set2',3);
    projmap = cbrewer('seq','Reds',3);
    hS = figure; 
    hR = figure; 
    for iF = 1:numel(fnames)
        figure(hS)
        plot3(All(iF).SponOnda01(:,1),All(iF).SponOnda01(:,2),All(iF).SponOnda01(:,3),'Color',[0 0 0]); hold on
        plot3(All(iF).SponOnda01(1,1),All(iF).SponOnda01(1,2),All(iF).SponOnda01(1,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);        
        plot3(All(iF).StimOnda01(:,1),All(iF).StimOnda01(:,2),All(iF).StimOnda01(:,3),'Color',projmap(iF,:)); hold on
                
        figure(hR)
        plot3(All(iF).SponPCprojU(:,1),All(iF).SponPCprojU(:,2),All(iF).SponPCprojU(:,3),'Color',[0 0 0]); hold on
        plot3(All(iF).SponPCprojU(1,1),All(iF).SponPCprojU(1,2),All(iF).SponPCprojU(1,3),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);        
        plot3(All(iF).StimPCprojU(:,1),All(iF).StimPCprojU(:,2),All(iF).StimPCprojU(:,3),'Color',projmap(iF,:)); hold on
        
        % build axes vectors
        Prep(iR).Prog(iF).PCvecs80 = All(iF).PCAdata.coeffs(1:All(iF).n*Prep(iR).Prog(iF).nPCsVAF80)';
        Prep(iR).Prog(iF).PCvecs95 = All(iF).PCAdata.coeffs(1:All(iF).n*Prep(iR).Prog(iF).nPCsVAF95)';
       
    end
    figure(hR); grid on; title(['Recording' num2str(iR) 'All on own axes'])
    figure(hS); grid on; title(['Recording' num2str(iR) 'All on da01 axes'])
    plot3(All(iF).ControlSponOnda01(:,1),All(iF).ControlSponOnda01(:,2),All(iF).ControlSponOnda01(:,3),'--','Color',[0.6 0.6 0.6]); hold on
    plot3(All(iF).ControlSponOnda01(1,1),All(iF).ControlSponOnda01(1,2),All(iF).ControlSponOnda01(1,3),'o','MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6]);        
    plot3(All(iF).ControlStimOnda01(:,1),All(iF).ControlStimOnda01(:,2),All(iF).ControlStimOnda01(:,3),'Color',[0.6 0.6 0.6]); hold on
    

    for iF = 1:numel(fnames)
        % store for plotting in example figure
        Viz(iR,iF).ProjSponOnda01 = All(iF).SponOnda01(:,1:3); % for plotting only, so just 3-D!
        Viz(iR,iF).ProjStimOnda01 = All(iF).StimOnda01(:,1:3);
        Viz(iR,iF).ProjControlSponOnda01 = All(iF).ControlSponOnda01(:,1:3);
        Viz(iR,iF).ProjControlStimOnda01 = All(iF).ControlStimOnda01(:,1:3);
    end
    
    % (1b): for projections on to same axis (da01), get min. distance between points on
    % attractor; compare to min distance to shuffle control
    % IF all 3 are actually moving on the attractor sufficiently to
    % estimate this!
    
    if all([RcrStats(iR,:).densAllPeriod] > spars.minDens)
        % initialise storage
        tmax = numel(All(1).Data.stimbins);
        Prep(iR).Dmin = zeros(numel(fnames),numel(fnames),tmax) + nan;  % set as Nans = non admissible points
        DminCtrl = zeros(numel(fnames),numel(fnames),spars.nShuffles,tmax) + nan; 
  
        ctr = 1;
        figure
        subplot(2,2,4); title(['Recording #' num2str(iR)]) 
        for i = 1:numel(fnames)
            % set of admissible recurrent points on target
            targetPCs = All(i).StimOnda01(All(i).ixAdPts,1:Prep(iR).Prog(1).nPCsVAF80);  % target program
            % targetts = All(i).Data.stimbins(All(i).ixAdPts); % corresponding times of admissible points
            
            
            % check against other programs 
            for j = i+1:numel(fnames)   
                
                mPeriod = (RcrStats(iR,i).mRot + RcrStats(iR,j).mRot) / 2; % mean of period
                radius = mPeriod/spars.divPeriod; % radius to check for distance
                
                % set of admissible time-points in source program
                srcPCs = All(j).StimOnda01(All(j).ixAdPts,1:Prep(iR).Prog(1).nPCsVAF80);
                srcts = All(j).Data.stimbins(All(j).ixAdPts); % corresponding times of admissible points
                targetts = All(i).Data.stimbins(All(i).ixAdPts); % corresponding times of admissible points
                
                % find common start point of target and source programs
                cmnStrt = max(srcts(1),targetts(1));
                ixSrc = find(srcts >= cmnStrt,1);
                ixTgt = find(targetts >= cmnStrt,1);
                srcts(srcts < cmnStrt) = []; targetts(targetts < cmnStrt) = []; % remove all time points before joint convergence...
                
                % for every point in source attractor, find min distance to
                % any point in target attractor, within period...
                ix = zeros(numel(srcts),1);
                for iN = 1:numel(srcts)
                    pt0 = srcPCs(iN+ixSrc-1,:);  % current x_i
                    ptts = srcts(iN);  % current time in src program
                    ixts = targetts < ptts + radius & targetts > ptts - radius;  % find all times in target within defined period
                    D = sqrt(sum(bsxfun(@minus,targetPCs(ixts,:),pt0).^2,2)); % Euclidean distance to all points within time radius
                    ix(iN) = find(All(j).Data.stimbins == ptts);  % store index for use in shuffled batch
                    Prep(iR).Dmin(i,j,ix(iN)) = min(D); % above diagonal: source -> target
                    
                end
                
                Dalltemp = zeros(spars.nShuffles,numel(srcts));
                parfor iB = 1:spars.nShuffles
                    %tic
                    % get shuffled target
                    ixR = randperm(All(i).n);
                    ControlStimOnda01 = zeros(numel(All(i).ixStim),Prep(iR).Prog(1).nPCsVAF80);
                    for iP = 1:Prep(iR).Prog(1).nPCsVAF80
                        ControlStimOnda01(:,iP) = sum(bsxfun(@times,All(i).SDF.spkfcn(All(i).ixStim,ixR),All(1).PCAdata.coeffs(:,iP)'),2); % shuffle neuron indices, and project onto da01 axes
                    end
                    ctrltargetPCs = ControlStimOnda01(All(i).ixAdPts,1:Prep(iR).Prog(1).nPCsVAF80); % control source: itself, shuffled
                    Dtemp = zeros(1,numel(srcts));
                    for iN = 1:numel(srcts)
                        pt0 = srcPCs(iN+ixSrc-1,:);  % current x_i
                        ptts = srcts(iN);  % current time in src program
                        ixts = targetts < ptts + radius & targetts > ptts - radius;  % find all times in target within defined period
                        Dctrl = sqrt(sum(bsxfun(@minus,ctrltargetPCs(ixts,:),pt0).^2,2)); % Euclidean distance to all points within time radius
                        % DminCtrl(i,j,iB,ix(iN)) = min(Dctrl); % below diagonal: source -> target*
                        Dtemp(iN) = min(Dctrl); % below diagonal: source -> target*    
                    end
                    Dalltemp(iB,:) = Dtemp;  % dealing with "slicing" issues
                    %toc
                end
                %for iB = 1:spars.nShuffles
                DminCtrl(i,j,:,ix) = Dalltemp;  % map Dtemp to correct indices...
                %end
                
                % and vice versa: from target to source
                ix = zeros(numel(targetts),1);
                for iN = 1:numel(targetts)
                    % current time....
                    pt0 = targetPCs(iN+ixTgt-1,:);  % current x_i
                    ptts = targetts(iN);
                    ixts = srcts < ptts + radius & srcts > ptts - radius;
                    D = sqrt(sum(bsxfun(@minus,srcPCs(ixts,:),pt0).^2,2)); % Euclidean distance to all points within time radius
                    ix(iN) = find(All(i).Data.stimbins == ptts);  % store index for use in shuffled batch
                    Prep(iR).Dmin(j,i,ix(iN)) = min(D); % below diagonal: target -> source
                end
                
                Dalltemp = zeros(spars.nShuffles,numel(targetts));
                parfor iB = 1:spars.nShuffles
                    % get shuffled source
                    ixR = randperm(All(j).n);
                    ControlStimOnda01 = zeros(numel(All(j).ixStim),Prep(iR).Prog(1).nPCsVAF80);                    
                    for iP = 1:Prep(iR).Prog(1).nPCsVAF80
                        ControlStimOnda01(:,iP) = sum(bsxfun(@times,All(j).SDF.spkfcn(All(j).ixStim,ixR),All(1).PCAdata.coeffs(:,iP)'),2); % shuffle neuron indices, and project onto da01 axes
                    end
                    ctrlsrcPCs = ControlStimOnda01(All(j).ixAdPts,1:Prep(iR).Prog(1).nPCsVAF80); % control source: itself, shuffled
                    Dtemp = zeros(1,numel(targetts));
                    for iN = 1:numel(targetts)
                        pt0 = targetPCs(iN+ixTgt-1,:);  % current x_i
                        ptts = targetts(iN);
                        ixts = srcts < ptts + radius & srcts > ptts - radius;
                        Dctrl = sqrt(sum(bsxfun(@minus,ctrlsrcPCs(ixts,:),pt0).^2,2)); % Euclidean distance to all points within time radius
                        % DminCtrl(j,i,iB,ix(iN)) = min(Dctrl); % below diagonal: target -> source*
                        Dtemp(iN) = min(Dctrl); % below diagonal: source -> target* 
                    end
                    Dalltemp(iB,:) = Dtemp;
                end
                DminCtrl(j,i,:,ix) = Dalltemp;  % map Dtemp to correct indices...
          
                % keyboard
                dd1 = squeeze(Prep(iR).Dmin(i,j,:));
                Prep(iR).MedDmin(i,j) = nanmedian(dd1);
                Prep(iR).MaxDmin(i,j) = nanmax(dd1);
                % propotion of distances that are within tolerance of being
                % on same manifold
                Prep(iR).prctDmin(i,j) = sum(dd1 <= spars.manifoldT)./numel(dd1);
               
                dd2 = squeeze(Prep(iR).Dmin(j,i,:));
                Prep(iR).MedDmin(j,i) = nanmedian(dd2);
                Prep(iR).MaxDmin(j,i) = nanmax(dd2);
                % propotion of distances that are within tolerance of being
                % on same manifold
                Prep(iR).prctDmin(j,i) = sum(dd2 <= spars.manifoldT)./numel(dd2);
             
                
                Prep(iR).HaussD(i,j) = max(nanmax(dd1),nanmax(dd2)); % maximum of minium distances between sets
                
                
                dctrl1 = squeeze(DminCtrl(i,j,:,:));
                Prep(iR).MedDminCtrl(i,j,:) = nanmedian(dctrl1,2);
                Prep(iR).MaxDminCtrl(i,j,:) = nanmax(dctrl1,[],2);
                Prep(iR).prctDminCtrl(i,j,:) = sum(dctrl1 <= spars.manifoldT,2)./numel(dctrl1);
                
                dctrl2 = squeeze(DminCtrl(j,i,:,:));
                Prep(iR).MedDminCtrl(j,i,:) = nanmedian(dctrl2,2);
                Prep(iR).MaxDminCtrl(j,i,:) = nanmax(dctrl2,[],2);
                Prep(iR).prctDminCtrl(j,i,:) = sum(dctrl2 <= spars.manifoldT,2)./numel(dctrl2);
                
                Prep(iR).HaussDCtrl(i,j,:) = max([nanmax(dctrl1,[],2),nanmax(dctrl2,[],2)],[],2); % maximum of minium distances between sets

                % plot of whole distribution of min distance
                [f,x] = ecdf(dd1);
                [f2,x2] = ecdf(dd2);
                [fctrl,xctrl] = ecdf(nanmean(dctrl1));   % mean over shuffles
                [fctrl2,xctrl2] = ecdf(nanmean(dctrl2));
               
                subplot(2,2,ctr),stairs(x,f,'r'); hold on
                stairs(x2,f2,'k')
                stairs(xctrl,fctrl,'r--')
                stairs(xctrl2,fctrl2,'k--')
                ctr = ctr+1;
                
               % keyboard
                
            end
        end
    end


    %% (2) correlate PC axes: separately, and as group (up to common number between 3 for 80% VAF) 
    Prep(iR).D = []; Prep(iR).L = [];
    for iF = 1:numel(fnames)
        % over all 3. build design matrix for correlation
        Prep(iR).D = [Prep(iR).D Prep(iR).Prog(iF).PCvecs80(1:All(iF).n*maxD)];
        Lt = [];
        for iD = 1:maxD 
            % PC loadings: divided by sqrt(Ev)
            Lt = [Lt; Prep(iR).Prog(iF).PCvecs80(1+All(iF).n*(iD-1):All(iF).n*iD)./  sqrt(All(iF).PCAdata.eigvalues(iD))];
        end
        Prep(iR).L = [Prep(iR).L Lt];
    end
%     % correlate stuff - no point unless we correct for sign flips
%     [rD,pD] = corr(D,'type','Spearman');
%     [rL,pL] = corr(L,'type','Spearman');
% 
    % each separate PC's correlation 
    for i = 1:numel(fnames)
        A = reshape(Prep(iR).Prog(i).PCvecs80(1:All(iF).n*maxD),All(iF).n,maxD);
        % AL = ...
        for j = 1:numel(fnames)
            B = reshape(Prep(iR).Prog(j).PCvecs80(1:All(iF).n*maxD),All(iF).n,maxD);
            [Prep(iR).rAll(i,j).rho,Prep(iR).pAll(i,j).prho] = corr(A,B,'type','Spearman');  % all possible corrs between PC weights in pair of programs
            Prep(iR).rhoAllD(:,i,j) = diag(Prep(iR).rAll(i,j).rho);  % same-same corrs on diagonal
            [Prep(iR).rAll(i,j).r,Prep(iR).pAll(i,j).p] = corr(A,B,'type','Pearson');  % all possible corrs between PC weights in pair of programs
            Prep(iR).rAllD(:,i,j) = diag(Prep(iR).rAll(i,j).r);  % same-same corrs on diagonal
        end
        
    end

    
%     %% visualise only up to maxD...
%     % Dmap = cbrewer('seq','Purples',pars.VmaxD);
%     Dmap = cbrewer('qual','Set2',spars.VmaxD);
%     matR = zeros(3,spars.VmaxD);
%     hscat = figure;
%     hpar = figure;
%     % plot each component separately
%     for iD = 1:spars.VmaxD 
%         figure(hpar)
%         r = squeeze(Prep(iR).rhoAllD(iD,:,:));
%         matR(:,iD) = [r(1,2),r(1,3),r(2,3)]';
%         plot(1:3,abs(matR(:,iD)),'.-','Color',Dmap(iD,:)); hold on
%         
%         figure(hscat)
%         subplot(2,2,1),plot(Prep(iR).D(1+All(iF).n*(iD-1):All(iF).n*iD,1),Prep(iR).D(1+All(iF).n*(iD-1):All(iF).n*iD,2),'.','Color',Dmap(iD,:)); hold on
%         xlabel('da01'); ylabel('da02'); 
%         
%         subplot(2,2,2),plot(Prep(iR).D(1+All(iF).n*(iD-1):All(iF).n*iD,1),Prep(iR).D(1+All(iF).n*(iD-1):All(iF).n*iD,3),'.','Color',Dmap(iD,:)); hold on 
%         xlabel('da01'); ylabel('da03')
%         
%         subplot(2,2,4),plot(Prep(iR).D(1+All(iF).n*(iD-1):All(iF).n*iD,2),Prep(iR).D(1+All(iF).n*(iD-1):All(iF).n*iD,3),'.','Color',Dmap(iD,:)); hold on
%         xlabel('da02'); ylabel('da03')
%     end
%     subplot(2,2,3), title(['Recording' num2str(iR) 'PCs1-3 weight correlations'])
%     
% %     figure(hscat)
% %     title(['Recording' num2str(iR) 'scatter correlation'])
%     
%     figure(hpar)
%     set(gca,'XTick',[1,2,3],'XTickLabel',{'da01 vs da02','da01 vs da03','da02 vs da03'})
%     ylabel('PC weight correlation')
%     title(['Recording #' num2str(iR) ' Parallel Coord correlation'])
    
    % legend('PC1','PC2','PC3')
    
%     figure
%     h = parallelcoords(abs(matR),'labels',{'da01 vs da02','da01 vs da03','da02 vs da03'},...
%             'Marker','.');

    % keyboard
    
%     figure
%     subplot(2,2,1),plot(L(:,1),L(:,2),'k.'); xlabel('da01'); ylabel('da02')
%     subplot(2,2,2),plot(L(:,1),L(:,3),'k.'); xlabel('da01'); ylabel('da03')
%     subplot(2,2,4),plot(L(:,2),L(:,3),'k.'); xlabel('da02'); ylabel('da03')
%     subplot(2,2,3), title(['Recording' num2str(iR) 'PCs1-3 loading correlations'])

    % (3) correlate Cxy
    Prep(iR).rhoCxy = triu(corr([All(1).Cxy(:), All(2).Cxy(:), All(3).Cxy(:)],'type','Spearman'),1);
    Prep(iR).rCxy = triu(corr([All(1).Cxy(:), All(2).Cxy(:), All(3).Cxy(:)],'type','Pearson'),1);
    
    for iF = 1:numel(fnames) % null model for correlation
        Sxy = All(iF).Cxy;
        Sxy(Sxy < 0) = 0; Sxy(eye(size(Sxy))==1) = 0;
        All(iF).Sxy = Sxy;
        All(iF).P = expectedA(Sxy);
        Prep(iR).rhoPxy(iF) = corr(All(iF).Sxy(:),All(iF).P(:),'type','Spearman');
        Prep(iR).rPxy(iF) = corr(All(iF).Sxy(:),All(iF).P(:),'type','Pearson');
        Viz(iR,iF).Sxy = Sxy;
        Viz(iR,iF).Pxy = All(iF).P;
    end
    
    Prep(iR).rhoSxy = triu(corr([All(1).Sxy(:), All(2).Sxy(:), All(3).Sxy(:)],'type','Spearman'),1);
    Prep(iR).rSxy = triu(corr([All(1).Sxy(:), All(2).Sxy(:), All(3).Sxy(:)],'type','Pearson'),1);
    
    toc
end

%% save stuff:
save([fname_Pre 'RecurrenceManifoldStats_AllPreps'],'Prep','Viz','spars') 


%% SUMMARISE

% cmap = cbrewer('qual','Paired',nfiles); 
cmap = brewermap(10,'Paired'); 
% summary plot: 
% (1) color code by prep to show consistency within prep
% (2) link within prep to draw eye to groups
% (3) markersize according to stim sequence: no patterns
% (4) markersize according to "sensitised"

Msize = [5 7 9];  % stim order
% Msize = [5 15 10];  % does sensitised prep show clear effect?


hPC= figure; hold on
hDs = figure; hold on
hCxy = figure; hold on
hDmin = figure; hold on

sym = 'o';

for iR = 1:nfiles
    % ID P10 preps using squares
%     if any(iR == P10preps')
%         sym = 's'; 
%     else
%         sym = 'o';
%     end
    % (1) plot all r for PC correlations (for each PC, and together)

    figure(hPC)
    for iD = 1:spars.VmaxD
        subplot(spars.VmaxD,1,iD), hold on
        rs = abs(squeeze(Prep(iR).rAllD(iD,:,:))); rhos = abs(squeeze(Prep(iR).rhoAllD(iD,:,:)));
        plot(rs(rs>0),rhos(rhos>0),'o-',...
            'Color',cmap(iR,:),'MarkerSize',10,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
        text(rs(1,2),rhos(1,2),num2str(iR),'Fontsize',8,'Color',[1 1 1])
        
%         % increasing MarkerSize for sequence
%         for iF = 1:numel(fnames)
%             plot(rs(rs>0),rhos(rhos>0),'o-',...
%                 'Color',cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',cmap(iR,:));
%             if iF == 3
%                 text([RcrStats(iR,iF).mRot],[RcrStats(iR,iF).mEigR],num2str(iR),'Fontsize',8,'Color',[1 1 1])
%             end
%         end
    end
    
    
    figure(hDs)
    rs1 = abs(squeeze(Prep(iR).rAllD(1,:,:))); rhos1 = abs(squeeze(Prep(iR).rhoAllD(1,:,:)));
    rs2 = abs(squeeze(Prep(iR).rAllD(2,:,:))); rhos2 = abs(squeeze(Prep(iR).rhoAllD(2,:,:)));
    plot(rs1(rs1>0),rs2(rs2>0),'o-',...
        'Color',cmap(iR,:),'MarkerSize',10,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
    text(rs1(1,2),rs2(1,2),num2str(iR),'Fontsize',8,'Color',[1 1 1])

    
    % (2) plot all r for Sxy and control model
    figure(hCxy)
%     plot(Prep(iR).rSxy(Prep(iR).rSxy>0),Prep(iR).rhoSxy(Prep(iR).rhoSxy > 0),'o-',...
%         'Color',cmap(iR,:),'MarkerSize',10,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
%     plot(Prep(iR).rPxy,Prep(iR).rhoPxy,'o-',...
%         'Color',cmap(iR,:),'MarkerSize',10,'MarkerFaceColor',[0.6 0.6 0.6],'Linewidth',0.5);

%     plot(Prep(iR).rPxy,Prep(iR).rSxy(Prep(iR).rSxy>0),'o-',...
%         'Color',cmap(iR,:),'MarkerSize',10,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
%    text(Prep(iR).rPxy(1),Prep(iR).rSxy(1,2),num2str(iR),'Fontsize',8,'Color',[1 1 1])
     ctrl =   Prep(iR).rhoPxy; data =  Prep(iR).rhoSxy(Prep(iR).rhoSxy>0);
     plot(ctrl,data,'o-',...
        'Color',cmap(iR,:),'MarkerSize',8,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
     text(Prep(iR).rhoPxy(3),Prep(iR).rhoSxy(2,3),num2str(iR),'Fontsize',8,'Color',[1 1 1])
%       for iF = 1:numel(fnames)
%             plot(ctrl(iF),data(iF),'o',...
%                 'Color',cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',cmap(iR,:));
%             if iF == 3
%                 text([RcrStats(iR,iF).mRot],[RcrStats(iR,iF).mEigR],num2str(iR),'Fontsize',8,'Color',[1 1 1])
%             end
%         end
%      % plot min distance
%      if ~isempty(Prep(iR).prctDmin)
%          prctDdata = 100*Prep(iR).prctDmin; prctDcontrol = 100*Prep(iR).prctDminCtrl;
%          figure(hDmin)
%          plot(prctDcontrol(prctDcontrol > 0),prctDdata(prctDdata > 0),'o-',...
%             'Color',cmap(iR,:),'MarkerSize',8,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
%          text(prctDcontrol(2,3),prctDdata(2,3),num2str(iR),'Fontsize',8,'Color',[1 1 1])
%      end

     % plot Hausdorff distance
     if ~isempty(Prep(iR).HaussD)
         ctrl = mean(Prep(iR).HaussDCtrl,3); ctrl = ctrl(ctrl > 0); 
         sem = std(Prep(iR).HaussDCtrl,0,3) / sqrt(spars.nShuffles); sem = sem(sem > 0);
         data = Prep(iR).HaussD(Prep(iR).HaussD > 0);
         figure(hDmin)
         plot(ctrl,data,'o',...
            'Color',cmap(iR,:),'MarkerSize',8,'MarkerFaceColor',cmap(iR,:),'Linewidth',0.5);
         line([ctrl-2*sem ctrl+2*sem]',[data data]','Linewidth',0.5,'Color',cmap(iR,:))
         % text(prctDcontrol(2,3),prctDdata(2,3),num2str(iR),'Fontsize',8,'Color',[1 1 1])
        % increasing MarkerSize for sequence
%         for iF = 1:numel(fnames)
%             plot(ctrl(iF),data(iF),'o',...
%                 'Color',cmap(iR,:),'MarkerSize',Msize(iF),'MarkerFaceColor',cmap(iR,:));
%             if iF == 3
%                 text([RcrStats(iR,iF).mRot],[RcrStats(iR,iF).mEigR],num2str(iR),'Fontsize',8,'Color',[1 1 1])
%             end
%         end
     end
end

figure(hDs)
xlabel('PC1 correlation')
ylabel('PC2 correlation')

figure(hCxy)
ylabel('Data correlation')
xlabel('Control correlation')
line([0 1],[0 1],'Color',[0.6 0.6 0.6])
axis([0 1 0 1])
% exportPPTfig(gcf,'Similarity_between_attractors',[10 15 7 7])

figure(hDmin)
title('Haussdorf distance')
ylabel('Data (spikes/s)')
xlabel('Control (spikes/s)')
axis square
axis([0 30 0 30])
line([0 30],[0 30],'Color',[0.6 0.6 0.6])

% title('Percentage points on same manifold')
% ylabel('Data (%)')
% xlabel('Control (%)')
% line([0 100],[0 100],'Color',[0.6 0.6 0.6])
% axis([0 100 0 100])


