% script to look at dimensionality at different VAF%

clear all; close all;

fnames = {'da01','da02','da03'};
prefix = [];

load(['../' prefix fnames{1} '_DataProperties_FunctionAndWindowSize'],'Data')

nPreps = numel(Data);

pars.VarExplained = [0.7 0.8 0.9 0.95];

%% get for variance explained
for iStim = 1:numel(fnames)
    load(['../' prefix fnames{iStim} '_DataProperties_FunctionAndWindowSize'],'PCAdata')

    for iPrep = 1:nPreps
        for iV = 1:numel(pars.VarExplained)
            % pick the subset of dimensions for whole thing 
            Projects(iPrep,iStim,iV).maxD = sum(PCAdata(iPrep).prctVar <= pars.VarExplained(iV));
            Projects(iPrep,iStim,iV).currPCs = PCAdata(iPrep).PCproj(:,1:Projects(iPrep,iStim,iV).maxD);  % from stimulation
        end
    end
end


for iV = 1:numel(pars.VarExplained)
    plotD(:,iV) = [Projects(:,:,iV).maxD];
end
    
figure
plot(pars.VarExplained*100,plotD,'.-')
xlabel('VAF(%)')
ylabel('Number of dimensions')

save Projects_For_VAF Projects pars