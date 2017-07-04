%% largest contributors to joint PCA
clear all; close all

addpath ../../Functions/

fname_Pre = '';
load RecurrenceManifoldStats_CommonAxes % joint PCA projection
fnames = {'da01','da02','da03'};

spars.VarExplained = 0.8;  % just get it to work....

nfiles = numel(Prep);

%% analyse
hStrong = figure; hWeak = figure; hMaxVar = figure; hMinVar = figure;
for iR = 1:nfiles
    iR
    
% Analysis only of projection onto the single concatenated axis...
    figure
    Coeffs(iR).sumW = sum(abs(Prep(iR).PCA.coeffs(:,1:Prep(iR).PCA.nPCs)),2);
    Coeffs(iR).sumWAll = sum(abs(Prep(iR).PCA.coeffs),2);
    
%     for iP = 1:3
%        subplot(5,1,iP),
%        bar(Prep(iR).PCA.coeffs(:,iP))
%     end
    [srt,Coeffs(iR).Rank] = sort(Coeffs(iR).sumW,'descend');
    subplot(3,1,1),
    bar(srt)
    ylabel('Participation: sum W'); xlabel('Ranked neurons');
    subplot(3,1,2)
    ecdf(Coeffs(iR).sumW)
    xlabel('Participation: sum W'); ylabel('P(participation)');

    subplot(3,1,3)
    plot(Coeffs(iR).sumW, Coeffs(iR).sumWAll,'k.');
    xlabel('Participation: sum W (VAF)'); ylabel('Participation: sum W (all PCs)');
    
    
    % hist(Coeffs(iR).sumW,20);
    conSDF = []; pre = 0;
    for iF = 1:numel(fnames)
       % load stuff and store across all 3 programs
       load(['../' fname_Pre fnames{iF} '_DataProperties_FunctionAndWindowSize']) % to get density functions and basic data
       All(iF).ixts = find(Data(iR).bins < spars.stimstart);
       All(iF).ixStim = find(Data(iR).bins >= spars.stimstart);
       All(iF).n = numel(Data(iR).rates);
       
       conSDF = [conSDF; SDF(iR).spkfcn(All(iF).ixStim,:)];
       % endProg(iF) = Prep(iR).PIdx(iF).ts(end);   % concatenated PCA
       endProg(iF) = pre + numel(All(iF).ixStim); pre = endProg(iF);
       
       All(iF).PCAdata = PCAdata(iR); 
       
       
       %% vital choice here?
       % retain VAF% dimensions in each
       % Coeffs(iR).nPCs(iF) = find(cumsum(PCAdata(iR).eigvalues./sum(PCAdata(iR).eigvalues)) >= spars.VarExplained,1,'first');
       
       % or just use the same number of dimension for all 3 programs...
       % (i.e. the number from the joint projection....
       Coeffs(iR).nPCs(iF) = Prep(iR).PCA.nPCs;
       

       % keyboard
       % compute weights = L1-norm
       % Coeffs(iR).sumWProg(iF,:) = sum(abs(All(iF).PCAdata.coeffs(:,1:Coeffs(iR).nPCs(iF))),2);   % linear summation
       matEigs = repmat(PCAdata(iR).eigvalues(1:Coeffs(iR).nPCs(iF))',All(iF).n,1);  % repmat to same size as below
       Coeffs(iR).sumWProg(iF,:) = sum(abs(matEigs .* All(iF).PCAdata.coeffs(:,1:Coeffs(iR).nPCs(iF))),2);   % weighted by eigenvalue
       
       % rank order
       [~,ix] = sort(Coeffs(iR).sumWProg(iF,:),'descend'); % get order of neurons
       Coeffs(iR).RankProg(iF,ix) = 1:All(iF).n; % assign them their ranks
       % normalised versions
       Coeffs(iR).normWMax(iF,:) = Coeffs(iR).sumWProg(iF,:) ./ max(Coeffs(iR).sumWProg(iF,:));  % normalised to maximum
       Coeffs(iR).normWSum(iF,:) = Coeffs(iR).sumWProg(iF,:) ./ sum(Coeffs(iR).sumWProg(iF,:));  % normalised to total
    
    end
    
    
%     % plot by global ranking - concatenated PCA only
%     mS = round(max(max(conSDF)))+1;
%     figure(hStrong)
%     subplot(nfiles,1,iR),plot(conSDF(:,Coeffs(iR).Rank(1:2))); hold on
%     line([endProg; endProg],[0 0 0; mS mS mS],'Color',[0 0 0])
%     if iR == 1 title('Maximum participation'); end
% 
%     % subplot(nfiles,2,iR+nfiles),plot(conSDF(:,ix(3:10)))
%     
%     figure(hWeak)
%     subplot(nfiles,1,iR),plot(conSDF(:,Coeffs(iR).Rank(end-1:end))); hold on
%     line([endProg; endProg],[0 0 0; mS mS mS],'Color',[0 0 0])
%     if iR == 1 title('Minimum participation'); end

   
    %% variable participation between programs? 
    
    % keyboard
%     [~,rnkix] = sort(Coeffs(iR).RankProg(1,:),'ascend');
%     figure
%     plot(Coeffs(iR).RankProg(:,rnkix),'k.-')
%     figure
%     imagesc(Coeffs(iR).RankProg(:,rnkix)')
    % variable participation by change in rank
    Coeffs(iR).dRank = [(Coeffs(iR).RankProg(1,:) - Coeffs(iR).RankProg(2,:))' (Coeffs(iR).RankProg(2,:) - Coeffs(iR).RankProg(3,:))']; 
    Coeffs(iR).maxDrank = max(abs(Coeffs(iR).dRank),[],2);  % largest jump of each neuron
    
    Coeffs(iR).minRank = min(Coeffs(iR).RankProg);  % highest rank obtained of each neuron
    
%     figure
%     subplot(211),plot(Coeffs(iR).minRank,Coeffs(iR).maxDrank,'k.')
%     xlabel('Highest Rank')
%     ylabel('Maximum change')
%     title(['Variable participation plot for Program ' num2str(iR)])
%     subplot(212),plot(Coeffs(iR).sumW,Coeffs(iR).maxDrank,'k.')
%     xlabel('Total contribution to joint PCA')
%     ylabel('Maximum change')
%     title(['Variable participation plot for Program ' num2str(iR)])
    
    % variable participation by change in relative proportion of maximum
    Coeffs(iR).dMaxP = [(Coeffs(iR).normWMax(1,:) - Coeffs(iR).normWMax(2,:))' (Coeffs(iR).normWMax(2,:) - Coeffs(iR).normWMax(3,:))' (Coeffs(iR).normWMax(1,:) - Coeffs(iR).normWMax(3,:))']; 
    Coeffs(iR).maxDMaxP = max(abs(Coeffs(iR).dMaxP(:,1:2)),[],2);  % largest jump of each neuron between consecutive programs
    
    Coeffs(iR).HighMaxP = max(Coeffs(iR).normWMax)';  % highest P(max) of each neuron
    Coeffs(iR).NmaxDmaxP = Coeffs(iR).maxDMaxP ./ Coeffs(iR).HighMaxP; % scale of change compared to maximum...
    
    % fit a model to quantify variability of program
    tbl = table(Coeffs(iR).HighMaxP,Coeffs(iR).maxDMaxP ,'VariableNames',{'Max_Pi','Max_DPi'});
    lm=fitlm(tbl,'linear');
    Coeffs(iR).Variability = lm.Coefficients.Estimate(2); 
    Coeffs(iR).VariabilityR2 = lm.Rsquared.Ordinary;
    lm=fitlm(tbl,'linear','RobustOpts','on');   % default: bisquare
    Coeffs(iR).VariabilityRobust = lm.Coefficients.Estimate(2); 
    Coeffs(iR).VariabilityR2Robust = lm.Rsquared.Ordinary;
 
    figure
    subplot(211),
    plot(Coeffs(iR).HighMaxP,Coeffs(iR).maxDMaxP,'k.')
    xlabel('Highest Proportion of Maximum')
    ylabel('Maximum change in proportion')
    title(['Variable participation plot for Program ' num2str(iR)])
    axis([0 1.05 0 1])
    subplot(212),
    plot(Coeffs(iR).HighMaxP,Coeffs(iR).NmaxDmaxP,'k.')
    xlabel('Highest Proportion of Maximum')
    ylabel('Scale of max change in proportion')
    title(['Variable participation plot for Program ' num2str(iR)])
    axis([0 1.05 0 1])
    
    % keyboard
    
    % plot most and least variable neurons in this prep...
    [Varsrt,Varrnk] = sort(Coeffs(iR).maxDMaxP,'descend'); % from most to least variable
    mS = round(max(max(conSDF)))+1;
    figure(hMaxVar)
    subplot(nfiles,1,iR),plot(conSDF(:,Varrnk(1:2))); hold on
    line([endProg; endProg],[0 0 0; mS mS mS],'Color',[0 0 0])
    if iR == 1 title('Maximum variation'); end
    figure(hMinVar)
    subplot(nfiles,1,iR),plot(conSDF(:,Varrnk(end-1:end))); hold on
    line([endProg; endProg],[0 0 0; mS mS mS],'Color',[0 0 0])
    if iR == 1 title('Minimum variation'); end

%     figure
%     plot(conSDF(:,20)); hold on
%     line([endProg; endProg],[0 0 0; mS mS mS],'Color',[0 0 0])
 
%     %%  variable phase between programs? WRT to all others 
%     for iF = 1:numel(fnames)  
%         for iM = 1:All(1).n
%             %H1 = hilbert(conSDF(Prep(iR).PIdx(iF).ts,Coeffs(iR).Rank(iM)) - mean(conSDF(Prep(iR).PIdx(iF).ts,Coeffs(iR).Rank(iM))));
%             H1 = hilbert(conSDF(Prep(iR).PIdx(iF).ts,iM) - mean(conSDF(Prep(iR).PIdx(iF).ts,iM)));
%             phase1 = angle(H1); % phase on [-pi,pi]
% 
%             for iN = iM+1:All(1).n
%                 %H2 = hilbert(conSDF(Prep(iR).PIdx(iF).ts,Coeffs(iR).Rank(iN)) - mean(conSDF(Prep(iR).PIdx(iF).ts,Coeffs(iR).Rank(iN))));
%                 H2 = hilbert(conSDF(Prep(iR).PIdx(iF).ts,iN) - mean(conSDF(Prep(iR).PIdx(iF).ts,iN)));
%                 phase2 = angle(H2); % phase on [-pi,pi]
% 
%                 Pdiffs = phase1 -  phase2;  % instantaneous phase difference
%                 Project = exp(1i*(Pdiffs));  % wrapped to unit circle
%                 Rz = mean(Project);  % mean of the vectors
%                 Coeffs(iR).Pdiff(iM,iN,iF) = angle(Rz); % mean phase angle between signals: delay from rs to cs (+ve = rs ahead of cs)
%                 Coeffs(iR).PL(iM,iN,iF) = sqrt(real(Rz).^2+imag(Rz).^2); % phase synch
% 
%             end
%         end
%     end
%     
%     % summarise over all pairs of programs: change in phase locking and
%     % phase offset between pairs
%     for i = 1:numel(fnames)
%         PLvec1 = nonzeros(squeeze(Coeffs(iR).PL(:,:,i))); Phvec1 = nonzeros(squeeze(Coeffs(iR).Pdiff(:,:,i))); 
%         for j = i+1:numel(fnames) 
%             PLvec2 = nonzeros(squeeze(Coeffs(iR).PL(:,:,j))); Phvec2 = nonzeros(squeeze(Coeffs(iR).Pdiff(:,:,j))); 
%             Coeffs(iR).dPL(i,j,:) = abs(PLvec1 - PLvec2); % change in phase-locking
%             Coeffs(iR).mPL(i,j,:) = (PLvec1 + PLvec2)/2;  % mean phase-locking
%             Coeffs(iR).maxPL(i,j,:) = max(PLvec1,PLvec2);  % max phase-locking
%             dPhase = abs(Phvec1 - Phvec2);  % get all changes in phase-offsets; this will be on [0, 2pi];
%             dPhase(dPhase > pi) = abs(dPhase(dPhase > pi) - 2*pi); % wrap so that phase difference expressed as [0,pi]: pi is maximum possible phase-flip
%             Coeffs(iR).dPhase(i,j,:) = dPhase;  % difference in phase-offset
%             % keyboard
%             figure
%             subplot(211),plot(squeeze(Coeffs(iR).maxPL(i,j,:)),squeeze(Coeffs(iR).dPL(i,j,:)),'k.')
%             xlabel('Max PL'); ylabel('\Delta PL')
%             title(['Programs ' num2str(i) ' vs ' num2str(j) ': jumps in phase-locking?']);
% 
%             subplot(212),plot(squeeze(Coeffs(iR).mPL(i,j,:)),squeeze(Coeffs(iR).dPhase(i,j,:))/pi,'k.')
%             xlabel('Mean PL'); ylabel('Delta Phase (% of pi)')
%             title(['Programs ' num2str(i) ' vs ' num2str(j) ': jumps in phase-offset?']);
% 
% %             subplot(313),plot(squeeze(Coeffs(iR).dPL(i,j,:)),squeeze(Coeffs(iR).dPhase(i,j,:))/pi,'k.')
% %             xlabel('\Delta PL'); ylabel('Delta Phase (% of pi)')
%         end
%     end
    % keyboard
end

%% pool over all preps to get overall variability view...
hP = figure; allPmax = []; allDPmax = []; allN_DPmax = [];
hPhi = figure;

% plot population scatters, and collate for distributions
for iR = 1:nfiles
    figure(hP)
    subplot(211), hold on
    plot(Coeffs(iR).HighMaxP,Coeffs(iR).maxDMaxP,'k.')
    xlabel('Highest Proportion of Maximum')
    ylabel('Maximum change in proportion')
    title(['Variable participation plot for all programs'])
    axis([0 1.05 0 1])
    subplot(212), hold on
    plot(Coeffs(iR).HighMaxP,Coeffs(iR).NmaxDmaxP,'k.')
    xlabel('Highest Proportion of Maximum')
    ylabel('Scale of max change in proportion')
    title(['Variable participation plot for all programs'])
    axis([0 1.05 0 1])

    % find outliers of change in proportion, fitting Normal distribution
    data = Coeffs(iR).dMaxP(:,1:2);  % changes between consecutive programs
    [Out,Coeffs(iR).Pars,LogL,Lcheck] = iterateoutliers(data(:),'Normal');
    Coeffs(iR).MinVarPart = max(abs(Coeffs(iR).Pars.mu-3*Coeffs(iR).Pars.std),abs(Coeffs(iR).Pars.mu+3*Coeffs(iR).Pars.std)); % threshold based on extremes of fitted model
    Coeffs(iR).VarPart = find(Coeffs(iR).maxDMaxP >= Coeffs(iR).MinVarPart);  % all neurons who exceeded the threshold at least once
    
    % collate for population distributions
    allPmax = [allPmax; Coeffs(iR).normWMax(:)];
    allDPmax = [allDPmax; Coeffs(iR).dMaxP(:)];
    Normed = Coeffs(iR).dMaxP(:,1:2) ./  Coeffs(iR).normWMax(1:2,:)';
    allN_DPmax = [allN_DPmax; Normed(:)]; 
    
    % and for phase
    figure(hP)
    subplot(211), hold on
    plot(Coeffs(iR).HighMaxP,Coeffs(iR).maxDMaxP,'k.')
    xlabel('Highest Proportion of Maximum')
    ylabel('Maximum change in proportion')
    title(['Variable participation plot for all programs'])
    axis([0 1.05 0 1])
    subplot(212), hold on
    plot(Coeffs(iR).HighMaxP,Coeffs(iR).NmaxDmaxP,'k.')
    xlabel('Highest Proportion of Maximum')
    ylabel('Scale of max change in proportion')
    title(['Variable participation plot for all programs'])
    axis([0 1.05 0 1])
end

% find outliers of change in proportion
[Out,Pars,LogL,Lcheck] = iterateoutliers(allDPmax,'Normal');
% [OutLog,ParsLog,LogLLog,LcheckLog] = iterateoutliers(allDPmax,'Logistic');


x = -1:0.01:1;
y = pdf(Pars,x);
yc = cdf(Pars,x);

[hDPmax,bins] = hist(allDPmax,50);
figure
subplot(311),hist(allPmax,50);
xlabel('Proportion of maximum contribution')
ylabel('Count')
subplot(312),bar(bins,hDPmax/sum(hDPmax)); hold on
plot(x,y/sum(y),'r')
xlabel('Change in proportion of maximum contribution')
ylabel('Count')

subplot(313),hist(allN_DPmax,50);
xlabel('Normed Delta proportion of maximum contribution')
ylabel('Count')

% cumulative view
[f,fx] = ecdf(allDPmax);
figure
stairs(fx,f); hold on
stairs(x,yc,'r')



%% plot per program scatters...

% %% % what are the outliers if we apply population level model?
% vars = abs(allDPmax(Out));  % all outlier values: not in Gaussian
% MinVarPart = min(vars);  % threshold for identiyfing outliers
% MinVarPart = max(abs(Pars.mu-3*Pars.std),abs(Pars.mu+3*Pars.std)); % threshold based on extremes of fitted model
% 
% % now find all neurons with max change exceeding this minium
% % for iR = 1:nfiles
% %     Coeffs(iR).VarPart = find(Coeffs(iR).maxDMaxP >= MinVarPart);
% % end


%% save data
save([fname_Pre 'ParticipationChanges_SamePCs'],'Coeffs');











