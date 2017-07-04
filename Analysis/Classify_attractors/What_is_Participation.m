%% contributions of synchrony to participation
clear all; close all

fnames = {'da01','da02','da03'};
fname_Pre = [];

load([fname_Pre 'ParticipationChanges_SamePCs'],'Coeffs');  % load Pariticpation scores

npreps = size(Coeffs,2);
nprogs = numel(fnames);

%% is participation just firing rates?

% correlation between participation and firing rate
 for iProg = 1:numel(fnames)
    for iPrep = 1:npreps
        load(['../' fname_Pre fnames{iProg} '_DataProperties_FunctionAndWindowSize'],'Data') % to get density functions and basic data
        PartContribution(iPrep).rates(iProg,:) = Data(iPrep).rates;
        [PartContribution(iPrep).rRatePart(iProg),PartContribution(iPrep).pRatePart(iProg)] = corr(Coeffs(iPrep).normWMax(iProg,:)',PartContribution(iPrep).rates(iProg,:)'); 
    end
 end

 % correlation between change in participation and change in firing rate
 for iPrep = 1:npreps
    PartContribution(iPrep).dRates = [PartContribution(iPrep).rates(1,:)' - PartContribution(iPrep).rates(2,:)' PartContribution(iPrep).rates(1,:)' - PartContribution(iPrep).rates(3,:)' PartContribution(iPrep).rates(2,:)' - PartContribution(iPrep).rates(3,:)'];
    for iChange = 1:numel(fnames)  % for each change, correlate the change in rate and the change in participation (normalised)
         [PartContribution(iPrep).rDrateDpart(iChange),PartContribution(iPrep).pDrateDpart(iChange)] = corr(PartContribution(iPrep).dRates(:,iChange),Coeffs(iPrep).dMaxP(:,iChange));
    end
    % correlate maximum
    PartContribution(iPrep).maxDrates = max(abs(PartContribution(iPrep).dRates(:,1:2)),[],2);
    [PartContribution(iPrep).rDrateDpartMax,PartContribution(iPrep).pDrateDpartMax] = corr(PartContribution(iPrep).maxDrates, Coeffs(iPrep).maxDMaxP);
 end



 
%% correlation between participation and total synchrony
 for iProg = 1:numel(fnames)
    for iPrep = 1:npreps
        load(['../Ensembles/' fname_Pre fnames{iProg} '_StaticEnsembles']) % to get Cxy matrices
        SxyP = StaticEnsemblesAll(iPrep).Cxy; SxyN = SxyP;
        SxyP(StaticEnsemblesAll(iPrep).Cxy < 0) = 0;
        SxyN(StaticEnsemblesAll(iPrep).Cxy > 0) = 0;
        PartContribution(iPrep).TSxy(iProg,:) = max([abs(sum(SxyN));sum(SxyP) - 1]);  % total synchrony as max total; except self (on diagonal)
        PartContribution(iPrep).TSxy(iProg,:) = abs(sum(SxyN)) + sum(SxyP) - 1;  % total synchrony, except self (on diagonal)

        PartContribution(iPrep).NormTSxy(iProg,:) = PartContribution(iPrep).TSxy(iProg,:) ./ max(PartContribution(iPrep).TSxy(iProg,:));
        
        [PartContribution(iPrep).rSynchPart(iProg),PartContribution(iPrep).pSynchPart(iProg)] = corr(Coeffs(iPrep).normWMax(iProg,:)',PartContribution(iPrep).TSxy(iProg,:)'); 
      
    end
 end

 
 
%% correlation between change in participation and change in total synchrony
for iPrep = 1:npreps
    PartContribution(iPrep).dSynch = [PartContribution(iPrep).TSxy(1,:)' - PartContribution(iPrep).TSxy(2,:)' PartContribution(iPrep).TSxy(1,:)' - PartContribution(iPrep).TSxy(3,:)' PartContribution(iPrep).TSxy(2,:)' - PartContribution(iPrep).TSxy(3,:)'];
    
    % add relative change here...
    
    for iChange = 1:numel(fnames)  % for each change, correlate the change in rate and the change in participation (normalised)
         [PartContribution(iPrep).rDsynchDpart(iChange),PartContribution(iPrep).pDsynchDpart(iChange)] = corr(PartContribution(iPrep).dSynch(:,iChange),Coeffs(iPrep).dMaxP(:,iChange));
    end
end



%% save

save([fname_Pre 'WhatIsParticipation_SamePCs'],'PartContribution');