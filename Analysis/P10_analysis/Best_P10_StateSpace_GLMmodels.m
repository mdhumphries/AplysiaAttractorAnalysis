%% for each P10 GLM model from state-space dynamics, get:
% (1) the best model (CV by forecasting)
% (1b) it's time-series of filters: consistency? 
% (2) Fit to first 50%; 
% (3) do Likelihood Ratio on coefficients
% (4) Forecast next 50% (do with and without LR coefficients)
% (5) Examine filters after LR etc key time-points and key dimensions in
% attractor (need N>2?)
%
% Mark Humphries
% 24/6/2016: added summary stats for training fits

clear all; close all

if ispc
    datapath = 'C:\Users\mqbssmhg.DS\Dropbox\paper2\P10\';
else
    datapath = '/Users/mqbssmhg/Dropbox/paper2/P10/';
end
    % spikepath = '/Users/mqbssmhg/Dropbox/paper2/datasets/deconcatenated_spks/experimental/';    
addpath ../Functions/

fname = {'da01','da02','da03'}; 
prefix = [];
VAFsuffix = []; % '_VAF_90'; % [];  % make sure this matches Bpars.VarIx, below...
Bpars.VarIx = 2;  % 2= 80%, 3 = 90%

% viz parameters
Divs = 4;  % division of fitting time-series into groups (4 = quartiles etc)
alpha = 0.001;  % signifiance level for likelihood-ratio
        
% get models and fitting parameters
load([fname{1} '_StateSpace_P10_GLMmodel' VAFsuffix]);
npreps = numel(P10data);

% % % get basic data
% load(['../' prefix fname '_DataProperties_FunctionAndWindowSize'],'FileTable','Data')
% get PCA projections used
load(['../Classify attractors/Projects_For_VAF.mat'])

color = 'rgb';
format = 'png';
dpi = 800;

ixexprt = [inf, inf]; % [2,1; 2,2; 3,3];  % examples to export (P10 prep, P10 stim]

%% get best models: for each prep, loop over programs evoked

for iF = 1:npreps
    for iP = 1:numel(fname)
        load([fname{iP} '_StateSpace_P10_GLMmodel' VAFsuffix],'P10data','pars'); %
        load(['../' prefix fname{iP} '_DataProperties_FunctionAndWindowSize'],'Data')
        % get recurrence point analysis, including state-space (currPCs)

        D = size(P10data(iF).Hist(1).filters,1); % dimensionality of the current model
        mid = floor(D/2);
        
        %% plot the forecast scores, and choose a model
        % combination of min MAE and max R (how to combine)
        % min MAE forecast
        % figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[20 1 20 30]); 
        figure
        subplot(6,D,1:mid)
        [ax,h1,h2] = plotyy(pars.Hmax,[P10data(iF).Hist(:).meanForecastMAE],pars.Hmax,[P10data(iF).Hist(:).meanForecastR]);
        xlabel(ax(1),'History (s)')
        ylabel(ax(1),'mean MAE')
        title('Forecast error')
        ylabel(ax(2),'mean R')

        % mean L of forecast?
        
        % best history model is lowest MAE fit
        ixbest = find([P10data(iF).Hist(:).meanForecastMAE] == min([P10data(iF).Hist(:).meanForecastMAE]));
        
        % time-series of fitting for best model
        T = Data(P10data(iF).iR).stimbins(end)  - pars.strtMatch;
        Traints = round(T * pars.prctTrain);    % length of test set  
        tsTrainStart = pars.strtMatch + (0:P10data(iF).Ntests-1) * (pars.Wstep*Data(P10data(iF).iR).GaussQt);
        tsForecastStart = tsTrainStart + Traints;
        subplot(6,D,mid+1:2*mid)
        [ax,h1,h2] = plotyy(tsForecastStart,P10data(iF).Hist(ixbest).mae_forecast,tsForecastStart,P10data(iF).Hist(ixbest).R_forecast);
        axis(ax,'tight')
        xlabel(ax(1),'forecast start (s)')
        ylabel(ax(1),'MAE')
        ylabel(ax(2),'R')
        title('Best model: time-series of forecast fits')       
        
        %% plot filters of best, and compute change over time (R of filter per training window)
        % should suggest dimensions contributing noise...
        Nhist = round(pars.Hmax(ixbest) * pars.Hstep);   % number of history bins
        n = floor(P10data(iF).Ntests / Divs);
        mD = zeros(Divs,Nhist);
        cmapAll = flipud(repmat([1:P10data(iF).Ntests]',1,3) * 1/P10data(iF).Ntests);
        cmapDivs = cbrewer('seq','Reds',Divs);  % flipud([[1:Divs]' zeros(Divs,1) [1:Divs]'] * 1/Divs);
        cmapCorr = flipud(cbrewer('div','RdBu',10)); % flip so red=1 blue = -1  
        for iD = 1:D
            for loop = 1:Divs    
                mD(loop,:) = mean(squeeze(P10data(iF).Hist(ixbest).filters(iD,1+n*(loop-1):loop*n,:)));
            end
            % every filter
            subplot(6,D,D+iD)
            set(gca, 'ColorOrder', cmapAll, 'NextPlot', 'replacechildren');
            plot(1:Nhist,squeeze(P10data(iF).Hist(ixbest).filters(iD,:,:)));
            title(['Filter#' num2str(iD) ': all filters'])
            
            % mean filters per quartile
            subplot(6,D,2*D+iD)
            set(gca, 'ColorOrder', cmapDivs, 'NextPlot', 'replacechildren');
            plot(1:Nhist,mD);
            xlabel('Time-step in past')    
            ylabel('filter weight')
            title(['Filter#' num2str(iD) ': quartile mean'])
            
            % correlation in time of filters
            subplot(6,D,3*D+iD)
            imagesc(tsTrainStart,tsTrainStart,P10data(iF).Hist(ixbest).filterCorr(iD).r); 
            colormap(cmapCorr)
            xlabel('Training start(s)'); ylabel('Training start(s)');
            title(['Filter#' num2str(iD) ' correlation'])
        end
        
  
        %% do first training window fit again
        % fit model to first X% of time-series, then test forecast on remaining set...
        ixStim = str2double(fname{iP}(end)); % column index
        SSts = Projects(P10data(iF).iR,ixStim,Bpars.VarIx).currPCs;
        SSDims =  Projects(P10data(iF).iR,ixStim,Bpars.VarIx).maxD;
        SSts(Data(P10data(iF).iR).stimbins < pars.strtMatch,:) = [];  % remove all potential stimulation artefacts; aligns both P10 and state-space time-series

        TrainStrt = pars.strtMatch;  
        TrainEnd = pars.strtMatch + Traints; % in seconds
        ixTrainStrt = find(P10data(iF).bins <= TrainStrt,1,'last');  % find indices of training set
        ixTrainEnd = find(P10data(iF).bins <= TrainEnd,1,'last');  % find indices of training set
        
        Tstep = round((pars.Hmax(ixbest)/Nhist) / Data(P10data(iF).iR).GaussQt);  % step between history points, in entries of time-series
        ixH =  Tstep:Tstep:Tstep*Nhist;  % history bin indices, back in time from target bin

        % build design matrix...
        Nt = ((ixTrainStrt-1) + (Tstep+pars.Hmax(ixbest)/Data(P10data(iF).iR).GaussQt)):ixTrainEnd; % indices of time points to use for design matrix
        Xflat = zeros(numel(Nt),D*Nhist); 
        % matrix of indices into SSts
        ixD = repmat(Nt',1,Nhist) - repmat(ixH,numel(Nt),1);  % index of target time-series, minus steps in past

        for iD = 1:numel(Nt)
            Xflat(iD,:) = reshape(SSts(ixD(iD,:),:),1,SSDims*Nhist);
        end

        % fit using general GLM model object, to get log-likelihood
        mdl = GeneralizedLinearModel.fit(Xflat, P10data(iF).spkfcn(Nt),'linear',...
                    'Distribution','poisson','Link','log');
        thetaF = mdl.Coefficients.Estimate;
        modelP10F = predict(mdl,Xflat);
        HmatF = reshape(thetaF(2:end),Nhist,D); % get matrix of filters
        
        %% do LR elimination of small coefficients
       % likelihood ratio test...
        Lmodel = sum(log(poisspdf(round(P10data(iF).spkfcn(Nt)),modelP10F)));  % likelihood of best fit model
        theta = thetaF(2:end);  % don't include bias term - only filters 
        [sTheta,I] = sort(abs(theta)); 
        ntheta = theta; 
        df = 0;
    
        for iC = 1:numel(ntheta)
            ntheta(I(iC)) = 0;  % test next weakest coefficient
            % [f,L0] = evalP10_fit(P10rate(Thist:end),X,Bsize,ntheta);
            modelP10new = glmval([thetaF(1); ntheta],Xflat,'log','constant','on'); 
            L0 = sum(log(poisspdf(round(P10data(iF).spkfcn(Nt)),modelP10new)));
            LR = -2* (L0 - Lmodel);  % likelihood ratio test statistic
            df = df + 1;  % difference in non-zero parameters between original and reduced model
            % critical value of Chi-squared: P = (1-alpha) of getting value greater than z by chance
            z = chi2inv(1-alpha,df); 

            if LR > z
                % then parameter is significantly contributing to model fit:
                % reinstate!
                ntheta(I(iC)) = theta(I(iC)); 
                df = df - 1;
                % AND STOP? (as all are sorted...)
            end
        end
        modelP10_final = glmval([thetaF(1); ntheta],Xflat,'log','constant','on'); 
        Rfinal = corr(P10data(iF).spkfcn(Nt),modelP10_final);        
        
        %% do forecast: with and without pruned filters
        
        if ~all(theta == ntheta)
            % some coefficients eliminated
            warning('Some coefficients were removed!')
        end
        
        % original model's forecast
        forecastNt = ixTrainEnd+1:numel(P10data(iF).bins); % all of the remaining P10 time-series
        Xfore = zeros(numel(forecastNt),D*Nhist); 
        ixD = repmat(forecastNt',1,Nhist) - repmat(ixH,numel(forecastNt),1);  % index of target time-series, minus steps in past

        for iD=1:numel(forecastNt)
            Xfore(iD,:) = reshape(SSts(ixD(iD,:),:),1,D*Nhist);
        end
        modelP10forecast = glmval(thetaF,[Xflat; Xfore],'log','constant','on');  % evaluate model forecast
     
        %% plot fit, and forecast
        subplot(6,D,4*D+(1:mid*2))
        plot(P10data(iF).bins,P10data(iF).spkfcn,'k'); hold on
        plot(P10data(iF).bins([Nt forecastNt]),modelP10forecast,'r')
        plot(P10data(iF).bins(Nt),modelP10F,'b')
        xlabel('Time (s)')
        ylabel('P10 firing rate (spikes/s)')
        title(['Fit and forecast of P10. Recording ' num2str(P10data(iF).iR) '; Program ' fname{iP}]);
        legend('P10',' Testing forecast','Training fit')
        
        %% plot filters of this specific model
        for iD = 1:D
            subplot(6,D,5*D+iD)
            line([0 Nhist],[0 0],'Color',[0.5 0.5 0.5]); hold on
            plot(1:Nhist,HmatF(:,iD),'k')
            title(['Filter#' num2str(iD) ': first filter'])
        end
        
        % exportfig(gcf,['P10_summary_' num2str(P10data(iF).iR) '_' fname{iP}],'Color',color,'Format',format,'Resolution',dpi)
        
        %% data for forecast summary figure
        SummaryStats(iF,iP).D = D; % dimensionality of state-space
        SummaryStats(iF,iP).meanMAE = mean(P10data(iF).Hist(ixbest).mae_forecast);
        SummaryStats(iF,iP).medianMAE = median(P10data(iF).Hist(ixbest).mae_forecast);
        SummaryStats(iF,iP).stdMAE = std(P10data(iF).Hist(ixbest).mae_forecast);
        SummaryStats(iF,iP).iqrMAE = prctile(P10data(iF).Hist(ixbest).mae_forecast,[25 75]);     
        
        SummaryStats(iF,iP).meanR = mean(P10data(iF).Hist(ixbest).R_forecast);
        SummaryStats(iF,iP).medianR = median(P10data(iF).Hist(ixbest).R_forecast);
        SummaryStats(iF,iP).stdR = std(P10data(iF).Hist(ixbest).R_forecast);
        SummaryStats(iF,iP).iqrR = prctile(P10data(iF).Hist(ixbest).R_forecast,[25 75]);     

        SummaryStats(iF,iP).meanR2 = mean(P10data(iF).Hist(ixbest).R_forecast.^2);
        SummaryStats(iF,iP).medianR2 = median(P10data(iF).Hist(ixbest).R_forecast.^2);
        SummaryStats(iF,iP).stdR2 = std(P10data(iF).Hist(ixbest).R_forecast.^2);
        SummaryStats(iF,iP).iqrR2 = prctile(P10data(iF).Hist(ixbest).R_forecast.^2,[25 75]);     
        
        %% data for fit summary figure
        SummaryStats(iF,iP).Train.meanMAE = mean(P10data(iF).Hist(ixbest).mae_train);
        SummaryStats(iF,iP).Train.medianMAE = median(P10data(iF).Hist(ixbest).mae_train);
        SummaryStats(iF,iP).Train.stdMAE = std(P10data(iF).Hist(ixbest).mae_train);
        SummaryStats(iF,iP).Train.iqrMAE = prctile(P10data(iF).Hist(ixbest).mae_train,[25 75]);     
        
        SummaryStats(iF,iP).Train.meanR = mean(P10data(iF).Hist(ixbest).R_train);
        SummaryStats(iF,iP).Train.medianR = median(P10data(iF).Hist(ixbest).R_train);
        SummaryStats(iF,iP).Train.stdR = std(P10data(iF).Hist(ixbest).R_train);
        SummaryStats(iF,iP).Train.iqrR = prctile(P10data(iF).Hist(ixbest).R_train,[25 75]);     

        SummaryStats(iF,iP).Train.meanR2 = mean(P10data(iF).Hist(ixbest).R_train.^2);
        SummaryStats(iF,iP).Train.medianR2 = median(P10data(iF).Hist(ixbest).R_train.^2);
        SummaryStats(iF,iP).Train.stdR2 = std(P10data(iF).Hist(ixbest).R_train.^2);
        SummaryStats(iF,iP).Train.iqrR2 = prctile(P10data(iF).Hist(ixbest).R_train.^2,[25 75]);     
   
        
        %% export examples
        if any(find(iF == ixexprt(:,1)) == find(iP == ixexprt(:,2)))
            figure
            plot(P10data(iF).bins,P10data(iF).spkfcn,'k'); hold on
            plot(P10data(iF).bins([Nt forecastNt]),modelP10forecast,'r')
            plot(P10data(iF).bins(Nt),modelP10F,'b')
            xlabel('Time (s)')
            ylabel('P10 firing rate (spikes/s)')
            exportPPTfig(gcf,['Prep' num2str(iF) '_Program' num2str(iP) ... 
                    '_full_forecast_R' num2str(SummaryStats(iF,iP).medianR) '_MAE' num2str(SummaryStats(iF,iP).medianMAE) '.png'],[10 15 6 4])
            
            close
            
            figure
            plot(P10data(iF).bins,P10data(iF).spkfcn,'k'); hold on
            plot(P10data(iF).bins([forecastNt(1:1000)]),...
                modelP10forecast(numel(modelP10F)+1:numel(modelP10F)+1000),'r')
            plot(P10data(iF).bins(Nt),modelP10F,'b')
            xlabel('Time (s)')
            ylabel('P10 firing rate (spikes/s)')
            exportPPTfig(gcf,['Prep' num2str(iF) '_Program' num2str(iP) ... 
                    '_example_forecast_R' num2str(SummaryStats(iF,iP).medianR) '_MAE' num2str(SummaryStats(iF,iP).medianMAE) '.png'],[10 15 6 4])
            close

            
            figure
            [ax,h1,h2] = plotyy(tsForecastStart,P10data(iF).Hist(ixbest).mae_forecast,tsForecastStart,P10data(iF).Hist(ixbest).R_forecast);
            axis(ax,'tight')
            xlabel(ax(1),'forecast start (s)')
%             ylabel(ax(1),'MAE')
%             ylabel(ax(2),'R')
            exportPPTfig(gcf,['Prep' num2str(iF) '_Program' num2str(iP) ... 
                    '_forecastfits_R' num2str(SummaryStats(iF,iP).medianR) '_MAE' num2str(SummaryStats(iF,iP).medianMAE) '.png'],[10 15 5 5])
            close
            
            %% plot filters of this specific model
            figure
            thist = pars.Hmax(ixbest) / Nhist:pars.Hmax(ixbest) / Nhist:pars.Hmax(ixbest);
            for iD = 1:D
                subplot(1,D,iD)
                line([0 thist(end)],[0 0],'Color',[0.5 0.5 0.5]); hold on
                plot(thist,HmatF(:,iD),'k')
                xlabel('History length (s)')
                ylabel('Weight')
                axis tight
            end
            exportPPTfig(gcf,['Prep' num2str(iF) '_Program' num2str(iP) ... 
                    '_examplefilters_R' num2str(SummaryStats(iF,iP).medianR) '_MAE' num2str(SummaryStats(iF,iP).medianMAE) '.png'],[10 15 15 4])
            close
        end
        
        %% store stuff for replotting above visualisations
        VizP10(iF,iP).Nt = Nt;
        VizP10(iF,iP).forecastNt = forecastNt;
        VizP10(iF,iP).tsP10Model = modelP10forecast;
        
        VizP10(iF,iP).tsForecastStart = tsForecastStart;
        VizP10(iF,iP).ixbest = ixbest;  % chosen filter set for visualisation: min of MAE
        VizP10(iF,iP).tsMAEforecast = P10data(iF).Hist(ixbest).mae_forecast;
        VizP10(iF,iP).tsRforecast = P10data(iF).Hist(ixbest).R_forecast;

        
    end
end

% do summary figure over all....
iqrMAE = reshape([SummaryStats(:).iqrMAE],2,npreps*numel(fname));
iqrR = reshape([SummaryStats(:).iqrR],2,npreps*numel(fname));

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 7]); hold on;
line([0,11],[0,0],'Color',[0.5 0.5 0.5]); hold on
%errorbar(1:npreps*numel(fname),[SummaryStats(:).meanMAE],[SummaryStats(:).stdMAE],'k.'); hold on
errorbar(0.2+(1:npreps*numel(fname)),[SummaryStats(:).medianMAE],...
                [SummaryStats(:).medianMAE]- iqrMAE(1,:),iqrMAE(2,:)-[SummaryStats(:).medianMAE],'k.'); 

view(90,90)
%pos = get(gca,'Position');
axis([0.5 9.5 0 50])
set(gca,'Position',[0.25 0.2 0.63 0.74])
ylabel('Median Absolute Error');
set(gca,'XTick',[])
% exportPPTfig(gcf,['Summary_MAE_median_IQR'],[10 15 4 7])

% xlabel('Session','FontSize',fontsize); 
% ylabel('Convergence (%)','FontSize',fontsize); 
% set(gca,'FontName','Helvetica','FontSize',fontsize);
% set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 4 7]); hold on;
line([0,11],[0,0],'Color',[0.5 0.5 0.5]); hold on
% errorbar(1:npreps*numel(fname),[SummaryStats(:).meanR],[SummaryStats(:).stdR],'k.'); hold on
errorbar(0.2+(1:npreps*numel(fname)),[SummaryStats(:).medianR],...
                [SummaryStats(:).medianR] - iqrR(1,:),iqrR(2,:)-[SummaryStats(:).medianR],'k.'); 
view(90,90)
axis([0.5 9.5 -0.4 1])
ylabel('R');
set(gca,'XTick',[])
set(gca,'Position',[0.25 0.2 0.63 0.74])
% exportPPTfig(gcf,['Summary_R_median_IQR'],[10 15 4 7])


figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 5]); hold on;
errorbarxy([SummaryStats(:).meanMAE],[SummaryStats(:).meanR],[SummaryStats(:).stdMAE],[SummaryStats(:).stdR],{'ko','r','r'}); hold on
% plot([SummaryStats(:).medianMAE],[SummaryStats(:).medianR],'k.')  %
% WAITING FOR IQR POINTS...

axis([0 35 -0.25 1])
xlabel('MAE'); ylabel('R')

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 5]); hold on;
% plot([SummaryStats(:).meanMAE],[SummaryStats(:).meanR],'r.'); hold on
plot([SummaryStats(:).medianMAE],[SummaryStats(:).medianR],'k.');
axis([0 35 0 1])
xlabel('Median Absolute Error'); ylabel('R')
% exportPPTfig(gcf,['Summary_R_vs_MAE'],[10 15 5 5])

save(['P10_statespace_Stats' VAFsuffix],'SummaryStats','VizP10')