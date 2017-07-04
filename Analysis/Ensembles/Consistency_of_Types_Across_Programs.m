%% compare each neuron's classification across programs within the same prep
%
% Types:
% 1: tonic
% 2: oscillator
% 3: burster
% 4: pauser
%
% Mark Humphries 24/5/2017

clear all; close all

fname_Pre = []; % [];
stimsets = {'da01','da02','da03'};  % 'da01': first; 'da02': second; and 'da03' : third ("rest" / control)

addpath ../../Functions/  % add list of local functions to path

load(['../' fname_Pre stimsets{1} '_DataProperties_FunctionAndWindowSize'],'DataTable');
nPreps = size(DataTable,1);

for iPrep = 1:nPreps
    % set up storage
    TypeConsistency(iPrep).NeuronType = zeros(DataTable(iPrep,1),numel(stimsets));

    for iProg = 1:numel(stimsets)
        % load ensemble statistics
        % groupdata: N-length struct, where N is the total number of
        % ensembles in the stimset(iProg} across all 10 preps
        % GroupList = [PrepID EnsembleID] for all N ensembles - gives index
        % into groupdata
        load([fname_Pre stimsets{iProg} '_Analyses_Neurons_and_Groups'],'groupdata','GroupList');

        load([fname_Pre stimsets{iProg} '_Ensemble_Types'],'AcorrTypes');
        
       
        % get this set of ensembles
        ixGrps = find(GroupList(:,1) == iPrep); % index into group data
        
        % look up their classification
        grpTypes = AcorrTypes(ixGrps);
        
        % assign each of their member neurons that classification
        for iG = 1:numel(ixGrps)
            TypeConsistency(iPrep).NeuronType(groupdata(ixGrps(iG)).IDs,iProg) = grpTypes(iG);
        end     
        
    end
    
    % summarise
    TypeConsistency(iPrep).FirstTransition = zeros(4); 
    TypeConsistency(iPrep).SecondTransition = zeros(4); 
    TypeConsistency(iPrep).BothTransitions = zeros(4); 
    for iType = 1:4
        ixFirst = find(TypeConsistency(iPrep).NeuronType(:,1) == iType);  % in first program, find neurons of this type
        ixSecond = find(TypeConsistency(iPrep).NeuronType(:,2) == iType);

        for iChange = 1:4  % for each type, find out if each neuron transtioned to it in the next program
            TypeConsistency(iPrep).FirstTransition(iType,iChange) =  sum(TypeConsistency(iPrep).NeuronType(ixFirst,2) == iChange);
            TypeConsistency(iPrep).SecondTransition(iType,iChange) =  sum(TypeConsistency(iPrep).NeuronType(ixSecond,3) == iChange);
        end
        
        TypeConsistency(iPrep).BothTransitions = TypeConsistency(iPrep).FirstTransition + TypeConsistency(iPrep).SecondTransition;
        
        % isolate consistently labelled neurons
        ixThird = find(TypeConsistency(iPrep).NeuronType(:,3) == iType);

        TypeConsistency(iPrep).MaxPropType(iType) = max([numel(ixFirst) numel(ixSecond) numel(ixThird)]) ./DataTable(iPrep,1);
        TypeConsistency(iPrep).PersistentTypeNeuronIDs{iType} = intersect(intersect(ixFirst,ixSecond),ixThird);
        
        % proportion of all neurons that were consistently identified as this Type
        TypeConsistency(iPrep).AllPropPersistent(iType) = numel(TypeConsistency(iPrep).PersistentTypeNeuronIDs{iType}) ./ DataTable(iPrep,1);

        % proportion of max identified Type neurons that were consistently this type
        TypeConsistency(iPrep).TypePropPersistent(iType) = numel(TypeConsistency(iPrep).PersistentTypeNeuronIDs{iType}) ./ TypeConsistency(iPrep).MaxPropType(iType);

    end
end

save EnsembleTypeConsistency TypeConsistency


