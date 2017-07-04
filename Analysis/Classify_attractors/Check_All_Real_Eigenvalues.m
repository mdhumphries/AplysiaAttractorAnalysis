% script to visualise other eigenvalues....
clear all; close all;


addpath ../../Functions/

% pick which threshold to look at
Eigpars.ixT = 3; % threshold of 10%
Eigpars.Same = 1e-7;

load da01_RecurrencePointsData_FilterPts RcrPtData
fname_Pre = '';
fnames = {'da01','da02','da03'};

nPreps = numel(RcrPtData);

cmap = brewermap(10,'Paired'); 
Msize = [3 9 5];   % for 5HT & washout: pick out the 5HT prep
Msize = [4 4 4];


for iPrep = 1:nPreps
   for iStim = 1:numel(fnames) 
       % load 
       load([fname_Pre fnames{iStim} '_RecurrencePointsData_FilterPts'],'RcrPtData','pars');
       load(['../' fname_Pre  fnames{iStim} '_DataProperties_FunctionAndWindowSize'],'Data','FileTable')

       %% get all eigenvalues
       nPoints = [RcrPtData(iPrep).Orbit(Eigpars.ixT).Period(:).nPoints];
       ixPrd = find(nPoints == max(nPoints));  
       
       % use only all eigenvalues from accurately estimated models
       usedPts = RcrPtData(iPrep).Orbit(Eigpars.ixT).Period(ixPrd).Npts >= pars.minNeighbours;
       AllEigs(iPrep,iStim).nUsedPts = numel(usedPts); % points used to compute each mean etc!
       allegs = RcrPtData(iPrep).Orbit(Eigpars.ixT).Period(ixPrd).allEgs(usedPts,:);  % the set of permissible eigenvalues
       allRealEgs = mean(real(allegs)); % take all means of the real parts

       % how do we determine which dimensions to skip? 
       % because means of real parts will be the same? 
       % No: because occasionally, there will not be an imaginary part
       % (during e.g. a perturbation without rotation) - so the eigenvalues
       % all change... Hence why we have stuck to the domiant eigenvalue
%        dReal = diff(allRealEgs); 
%        ix = dReal > -Eigpars.Same && dReal < Eigpars.Same;  % find all indices consisent with 0, given rounding error 
       
       % therefore Simple Version: just use all of them....              
       AllEigs(iPrep,iStim).mAllEigR = allRealEgs;
      
       %% Full Version: check all complex and not complex
       % keyboard
       blnReal = abs(imag(allegs))< eps('double').*abs(real(allegs)); % set of time-points that are only real in each dimension 
        
       AllEigs(iPrep,iStim).Pcomplex1 = sum(~blnReal(:,1)) ./ AllEigs(iPrep,iStim).nUsedPts; % proportion of dominant is complex
       AllEigs(iPrep,iStim).Pcontract1 = sum(real(allegs(:,1)) < 0) ./ AllEigs(iPrep,iStim).nUsedPts;  % proportion of time contracting
 
       AllEigs(iPrep,iStim).Pcomplex12 = sum(~blnReal(:,1) | ~blnReal(:,2)) ./ AllEigs(iPrep,iStim).nUsedPts; % proportion of top two dimensions in which at least one is complex
       AllEigs(iPrep,iStim).Pcontract12 = sum(real(allegs(:,1)) < 0 | real(allegs(:,2)) < 0) ./ AllEigs(iPrep,iStim).nUsedPts; % proportion of top two dimensions in which at least one is negative

       %% create list of unique magnitudes of dominant eigenvalues
       nDims = size(allegs,2);
       maxDomEgs = max(sum(~blnReal,2))/2 + rem(nDims,2); % maximum number of magnitude eigenvalues: +1 if number of dimensions is odd
       AllEigs(iPrep,iStim).UniqueMagnitudeEgs = zeros(AllEigs(iPrep,iStim).nUsedPts,maxDomEgs);
       for iDom = 1:AllEigs(iPrep,iStim).nUsedPts
            % for each recurrent point
            rCtr = 1; % counter for which recurrent points to check next
            for iDim = 1:maxDomEgs
                % assign values to each of the dominant eigenvalue columns in turn
                AllEigs(iPrep,iStim).UniqueMagnitudeEgs(iDom,iDim) = real(allegs(iDom,rCtr));  % magnitude of this entry
                if blnReal(iDom,rCtr)
                    rCtr = rCtr + 1;    % current entry is real, so check the next number
                else
                    rCtr = rCtr + 2;    % complex numbers in conjugate pairs: skip two
                end
            end
       end
       AllEigs(iPrep,iStim).mUniqueEigMag = mean(AllEigs(iPrep,iStim).UniqueMagnitudeEgs);
   end
end

%% plot
hPar = figure; hold on;
for iPrep = 1:nPreps
   for iStim = 1:numel(fnames) 
        % or store for parallel plot... 1->2->3->N 
        figure(hPar)
        plot(1:numel(AllEigs(iPrep,iStim).mAllEigR),AllEigs(iPrep,iStim).mAllEigR,'o-',...
            'Color',cmap(iPrep,:),'MarkerSize',Msize(iStim),'MarkerFaceColor',cmap(iPrep,:));
        line([0.5 8.5],[0 0],'Color',[0 0 0],'Linewidth',1)
   end
end
xlabel('Dimensions')
ylabel('Mean magnitude of eigenvalue')

hParUnique = figure; hold on;
for iPrep = 1:nPreps
   for iStim = 1:numel(fnames) 
        % or store for parallel plot... 1->2->3->N 
        figure(hParUnique)
        plot(1:numel(AllEigs(iPrep,iStim).mUniqueEigMag),AllEigs(iPrep,iStim).mUniqueEigMag,'o-',...
            'Color',cmap(iPrep,:),'MarkerSize',Msize(iStim),'MarkerFaceColor',cmap(iPrep,:));
        line([0.5 4.5],[0 0],'Color',[0 0 0],'Linewidth',1)
   end
end
xlabel('Unique magnitude eigenvalues')
ylabel('Mean magnitude of eigenvalue')


save CheckAllEigs AllEigs Eigpars