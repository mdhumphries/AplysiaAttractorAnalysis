%%% script to decode P10 from state-space using GLM

clear all; close all


datapath = '../../Data/P10/';
addpath ../../Functions/

fname = 'da03';
prefix = [];

% get PCA projections
load(['../' prefix fname '_DataProperties_FunctionAndWindowSize'],'FileTable','Data')  % filenames and 

% get PCA projections - has every recording in its matrix (10x3)
load(['../Classify attractors/Projects_For_VAF.mat'])
pars.VarIx = 3;  % 2= 80%, 3 = 90%
pars.Stim = str2double(fname(end)); % column index

% analysis parameters
pars.stimstart = 30;
pars.endwin = 123;  % allow for overrun of recording
pars.strtMatch = 33;

Bsize = 0.02;  % binsize in s  (for histogram - not used at present)

pars.Hstep = 100; % number of bins of spike history per second
pars.Hmax = 0.05:0.05:0.25;   % length of history (in s) 
% pars.Hstep = 5; % number of bins of spike history per second
% pars.Hmax = 1:1:4;   % length of history (in s) 

pars.Wstep = 20;  % number of bins to advance the forecasting window (1=10ms)

% GaussQt = 0.01;     % time-step for building P10 SDF
% GaussSD = 0.1;

pars.prctTrain = 0.5;  % min. proportion of time-series on which to train model
pars.kstep = 1:1:10; % forecast this many seconds ahead

%% get state-space

files = dir(datapath); files(1:2) = []; 

ixP10 = find(arrayfun(@(x) any(findstr(x.name,fname)),files)); % index of all fname files in P10 directory
tic
for iF = 1:numel(ixP10)
    iF
    % filename of P10 file
    P10data(iF).name =  files(ixP10(iF)).name;
        
    % get P10 spikes
    load([datapath P10data(iF).name]);
    P10spks = spks(:,2);
    % stimspks = P10spks(P10spks > pars.stimstart);  % predict stim-spikes
    
    % find matching dataset
    P10data(iF).iR = find(cellfun(@(x) strncmp(P10data(iF).name,x,7),FileTable));  % first 7 characters; index of matching dataset
    
    % retained time period of P10 spike train
    T = Data(P10data(iF).iR).stimbins(end)  - pars.strtMatch;
    % T = pars.endwin - pars.strtMatch;
    
%     % make full histogram of P10....
%     bins = pars.strtMatch:Bsize:pars.endwin;   % spike-counting bins, after stim offset
%     P10rate = histc(P10spks,bins);  % target rate of P10 to fit
    
    % convolve P10 with Gaussian 
    sig = Data(P10data(iF).iR).GaussSD / Data(P10data(iF).iR).GaussQt; % SD in time-steps
    x = [-5*sig:1:5*sig]';  % x-axis values of the discretely sampled Gaussian, out to 5xSD
    h = (1/(sqrt(2*pi*sig^2)))*exp(-((x.^2*(1/(2*sig^2))))); % y-axis values of the Gaussian
    h = h ./sum(h); % make sure kernel has unit area, then can use for rate functions
    shiftbase = floor(length(x)/2); 

    P10data(iF).bins = 0:Data(P10data(iF).iR).GaussQt:Data(P10data(iF).iR).stimbins(end);  % all time-points to evaluate
    
    % compute spike-train convolution functions
    [P10data(iF).spkfcn,idxs] = convolve_spiketrains([ones(numel(P10spks),1) P10spks],h,shiftbase,1,P10data(iF).bins,Data(P10data(iF).iR).GaussQt,[0 Data(P10data(iF).iR).stimbins(end)],'Gaussian');
    
    % constrict to after stimulation period
    P10data(iF).spkfcn= P10data(iF).spkfcn(P10data(iF).bins >= pars.strtMatch);
    P10data(iF).bins = P10data(iF).bins(P10data(iF).bins >= pars.strtMatch);    
    figure
    plot(P10data(iF).bins,P10data(iF).spkfcn)
    
    
    % fit model to first X% of time-series, then test forecast on remaining set...
    SSts = Projects(P10data(iF).iR,pars.Stim,pars.VarIx).currPCs;
    SSDims =  Projects(P10data(iF).iR,pars.Stim,pars.VarIx).maxD;
    % SSts = RcrPtData(P10data(iF).iR).currPCs;  % time-series in state-space
    SSts(Data(P10data(iF).iR).stimbins < pars.strtMatch,:) = [];  % remove all potential stimulation artefacts; aligns both P10 and state-space time-series
    
    stepsFore = round(pars.kstep ./ Data(P10data(iF).iR).GaussQt); % bins-in-the-future for each forecast time
    
    % number of forecast tests to do: starting from pars.strtMatch, and
    % advancing by Wstep bins each time
    Traints = round(T * pars.prctTrain);    % length of test set  
    P10data(iF).Ntests = floor((T - Traints - max(pars.kstep))/Data(P10data(iF).iR).GaussQt / pars.Wstep);  
    % P10data(iF).Ntests = 1; % check history effects on training fits
    
    for iH = 1:numel(pars.Hmax)
        iH
        Nhist = round(pars.Hmax(iH) * pars.Hstep);   % number of history bins
        Tstep = round((pars.Hmax(iH)/Nhist) / Data(P10data(iF).iR).GaussQt);  % step between history points, in entries of time-series
        ixH =  Tstep:Tstep:Tstep*Nhist;  % history bin indices, back in time from target bin

        P10data(iF).Hist(iH).filters = zeros(SSDims,P10data(iF).Ntests,Nhist);

        for iN = 1:P10data(iF).Ntests
            % iN
            % get window of training data
            TrainStrt = pars.strtMatch + (iN-1) * (pars.Wstep*Data(P10data(iF).iR).GaussQt);  % in seconds - advance by one history step
            TrainEnd = pars.strtMatch + Traints + (iN-1) * (pars.Wstep*Data(P10data(iF).iR).GaussQt); % in seconds
            ixTrainStrt = find(P10data(iF).bins <= TrainStrt,1,'last');  % find indices of training set
            ixTrainEnd = find(P10data(iF).bins <= TrainEnd,1,'last');  % find indices of training set

            % build design matrix...
            Nt = ((ixTrainStrt-1) + (Tstep+pars.Hmax(iH)/Data(P10data(iF).iR).GaussQt)):ixTrainEnd; % indices of time points to use for design matrix
            Xflat = zeros(numel(Nt),SSDims*Nhist); 
            % matrix of indices into SSts
            ixD = repmat(Nt',1,Nhist) - repmat(ixH,numel(Nt),1);  % index of target time-series, minus steps in past

            for iD = 1:numel(Nt)
                Xflat(iD,:) = reshape(SSts(ixD(iD,:),:),1,SSDims*Nhist);
            end

            % fit using general GLM model object, to get log-likelihood
            mdl = GeneralizedLinearModel.fit(Xflat, P10data(iF).spkfcn(Nt),'linear',...
                        'Distribution','poisson','Link','log');
%             P10data(iF).Hist(iH).LL_train2(iN) = mdl.LogLikelihood;
            thetaF = mdl.Coefficients.Estimate;
            modelP10F = predict(mdl,Xflat);
            
%             [thetaF,dev,stats] = glmfit(Xflat,P10data(iF).spkfcn(Nt),'poisson','link','log','constant','off');
%             modelP10F = glmval(thetaF,Xflat,'log','constant','off');
%             [thetaF,dev,stats] = glmfit(Xflat,P10data(iF).spkfcn(Nt),'poisson','link','log','constant','on'); % 'estdisp','on');
%             modelP10F = glmval(thetaF,Xflat,'log','constant','on');  % evaluate model fit
            

            
            P10data(iF).Hist(iH).bias(iN) = thetaF(1);
            HmatF = reshape(thetaF(2:end),Nhist,SSDims); % get matrix of filters
            for iD = 1:SSDims
                P10data(iF).Hist(iH).filters(iD,iN,:) = HmatF(:,iD);
            end
            
            % within-training correlation
            P10data(iF).Hist(iH).R_train(iN) = corr(P10data(iF).spkfcn(Nt),modelP10F);
            P10data(iF).Hist(iH).RMSE_train(iN) = sqrt(mean((P10data(iF).spkfcn(Nt) - modelP10F).^2));
            P10data(iF).Hist(iH).mae_train(iN) = median(abs(P10data(iF).spkfcn(Nt) - modelP10F));

                        
            % Use GLMVAL to compute fitted values for each observation. 
            % Then use the appropriate PDF function to compute individual likelihoods, 
            % then log and sum. So for a binomial GLM (e.g. logistic regression),
            P10data(iF).Hist(iH).LL_train(iN) = sum(log(poisspdf(round(P10data(iF).spkfcn(Nt)),modelP10F)));
           
            % k-step forecast
            forecastNt = ixTrainEnd+1:ixTrainEnd + max(pars.kstep)/Data(P10data(iF).iR).GaussQt;
            Xfore = zeros(numel(forecastNt),SSDims*Nhist); 
            ixD = repmat(forecastNt',1,Nhist) - repmat(ixH,numel(forecastNt),1);  % index of target time-series, minus steps in past

            for iD=1:numel(forecastNt)
                Xfore(iD,:) = reshape(SSts(ixD(iD,:),:),1,SSDims*Nhist);
            end
%             modelP10forecast = glmval(thetaF,[Xflat; Xfore],'log','constant','off');  % evaluate model forecast
%             modelP10forecast = glmval(thetaF,[Xflat; Xfore],'log','constant','on');  % evaluate model forecast
            
            modelP10forecast = predict(mdl,[Xflat; Xfore]);

            msteps = numel(Nt)+stepsFore;  % model forecast entries to compare
            P10data(iF).Hist(iH).Error_test(iN,:) = P10data(iF).spkfcn(forecastNt(stepsFore)) - modelP10forecast(msteps);
            P10data(iF).Hist(iH).R_forecast(iN) = corr(P10data(iF).spkfcn(forecastNt),modelP10forecast(numel(Nt)+1:end));
            P10data(iF).Hist(iH).mae_forecast(iN) = median(abs(P10data(iF).spkfcn(forecastNt) - modelP10forecast(numel(Nt)+1:end)));

            % log likelihood forecast
            P10data(iF).Hist(iH).LL_forecast(iN) = sum(log(poisspdf(round(P10data(iF).spkfcn(forecastNt)),modelP10forecast(numel(Nt)+1:end))));

%             % view the training and forecast fits
%             figure
%             plot(P10data(iF).bins,P10data(iF).spkfcn,'k'); hold on
%             plot(P10data(iF).bins([Nt forecastNt]),modelP10forecast,'r')
%             plot(P10data(iF).bins(Nt),modelP10F,'b')
            
        end
         
        
        % summary statistics of training and test fits
        P10data(iF).Hist(iH).rmseForecast_Points = sqrt(mean(P10data(iF).Hist(iH).Error_test.^2));
        P10data(iF).Hist(iH).maeForecast_Points = median(abs(P10data(iF).Hist(iH).Error_test)); % mean over each g
        P10data(iF).Hist(iH).meanForecastMAE = mean(P10data(iF).Hist(iH).mae_forecast);  % mean over all forecasts
        P10data(iF).Hist(iH).meanForecastR = mean(P10data(iF).Hist(iH).R_forecast); % mean over all forecasts
        P10data(iF).Hist(iH).meanTrainR = mean(P10data(iF).Hist(iH).R_train);  % mean over all training windows
        P10data(iF).Hist(iH).meanTrainMAE = mean(P10data(iF).Hist(iH).mae_train); % mean over all training windows
        
        % stability of filters over training
        for iD = 1:SSDims
            allfilters = squeeze(P10data(iF).Hist(iH).filters(iD,:,:));
            r = corr(allfilters'); % correlation of every filter with all others
            r(eye(P10data(iF).Ntests)==1) = 0;
            P10data(iF).Hist(iH).filterCorr(iD).r = r;
%             figure
%             imagesc(r); colormap(hot)
            P10data(iF).Hist(iH).filterStability(iD) = iqr(squareform(r)); % one option: others include?
            
        end
        % summarise: sum stability over all dimensions....
        P10data(iF).Hist(iH).TotalfilterStability = sum(P10data(iF).Hist(iH).filterStability); % one option: others include?

        
%         figure
%         plot(HmatF)
        
%         % time-series of training fits
%         tsTrainStart = pars.strtMatch + (0:P10data(iF).Ntests-1) * (pars.Wstep*Data(P10data(iF).iR).GaussQt);
%         figure
%         subplot(311),plot(tsTrainStart,P10data(iF).Hist(iH).R_train)
%         xlabel('training window start (s)')
%         ylabel('R of training fit to data')
%         subplot(312),plot(tsTrainStart,P10data(iF).Hist(iH).mae_train)
%         xlabel('training window start (s)')
%         ylabel('MAE of training fit to data')
%         subplot(313),plot(tsTrainStart,P10data(iF).Hist(iH).LL_train)
%         xlabel('training window start (s)')    
%         ylabel('log-L of training fit to data')
% 
%         % time-series of forecast correlation
%         tsForecastStart = tsTrainStart + Traints;
%         figure
%         subplot(311),plot(tsForecastStart,P10data(iF).Hist(iH).R_forecast)
%         xlabel('forecast start (s)')
%         ylabel('R of forecast fit to data')
%         subplot(312),plot(tsForecastStart,P10data(iF).Hist(iH).mae_forecast)
%         xlabel('forecast start (s)')
%         ylabel('MAE of forecast fit to data')        
%         subplot(313),plot(tsForecastStart,P10data(iF).Hist(iH).LL_forecast)
%         xlabel('forecast start (s)')    
%         ylabel('log-L of forecast fit to data')
% %         
%         % change in history kernels
%         Divs = 4; n = floor(P10data(iF).Ntests / Divs);
%         mD = zeros(Divs,Nhist);
%         cmapAll = flipud(repmat([1:P10data(iF).Ntests]',1,3) * 1/P10data(iF).Ntests);
%         cmapDivs = cbrewer('seq','Reds',Divs);  % flipud([[1:Divs]' zeros(Divs,1) [1:Divs]'] * 1/Divs);
%         for iD = 1:SSDims
%             for loop = 1:Divs    
%                 mD(loop,:) = mean(squeeze(P10data(iF).Hist(iH).filters(iD,1+n*(loop-1):loop*n,:)));
%             end
%             figure
%             subplot(211),
%             set(gca, 'ColorOrder', cmapAll, 'NextPlot', 'replacechildren');
%             plot(1:Nhist,squeeze(P10data(iF).Hist(iH).filters(iD,:,:)));
%             subplot(212)
%             set(gca, 'ColorOrder', cmapDivs, 'NextPlot', 'replacechildren');
%             plot(1:Nhist,mD);
%             xlabel('Time-step in past')    
%             ylabel('filter weight')
%         end
%         % bias
%         figure
%         plot(tsTrainStart,P10data(iF).Hist(iH).bias,'k')
%         xlabel('training window start (s)')
%         ylabel('bias term')
        
    end
    
    maeForecast = []; 
    for iH = 1:numel(pars.Hmax)
        maeForecast = [maeForecast; P10data(iF).Hist(iH).maeForecast_Points];
    end
    
    cmapH = cbrewer('seq','Reds',numel(pars.Hmax)); 

    figure
    subplot(321),
    set(gca, 'ColorOrder', cmapH, 'NextPlot', 'replacechildren');
    % plot(pars.Hmax,maeForecast,'Color',[0.5 0.5 0.5]); hold on
    plot(pars.Hmax,maeForecast); hold on
    plot(pars.Hmax,mean(maeForecast,2),'Color',[0 0 0])
    
    xlabel('History (s)')
    ylabel('Forecast error - points (MAE)')
    
    subplot(323),
    plot(pars.Hmax,[P10data(iF).Hist(:).meanForecastMAE])
    xlabel('History (s)')
    ylabel('Forecast error - all (mean MAE)')

    subplot(325),
    plot(pars.Hmax,[P10data(iF).Hist(:).meanTrainMAE])
    xlabel('History (s)')
    ylabel('Train error (MAE)')


    
    
    subplot(324),
    plot(pars.Hmax,[P10data(iF).Hist(:).meanForecastR])
    xlabel('History (s)')
    ylabel('Forecast fit (mean R)')
   
    subplot(326),
    %plot(pars.Hmax*Nhist*SSDims,[P10data(iF).Hist(:).meanTrainR2])
    % xlabel('No. coefficients')
    plot(pars.Hmax,[P10data(iF).Hist(:).meanTrainR])    
    xlabel('History (s)')
    ylabel('Train fit (mean R)')

    % keyboard
end
toc

strVar = num2str(round(100*pars.VarExplained(pars.VarIx)));
save([fname '_StateSpace_P10_GLMmodel_VAF_' strVar],'P10data','pars'); % ,'-v7.3');
