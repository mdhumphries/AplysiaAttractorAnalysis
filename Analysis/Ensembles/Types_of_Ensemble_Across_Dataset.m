%%% script to detect the types of ensembles in the data-set
%
% Key outputs:
% dataset_spikes: the matrix defining the fit-space - one row per ensemble,
%                 first set of columns the ISI P(model) vector, 
%                 the second set of columns the CV2 P(model) vector
% SxyData.Spikes: the similarity matrix for all pairs of ensembles in the
%                 fit-space
% Ccon.Spikes:    the ensemble-types identified by the consensus community
%                 detection algorithm (format: [EnsembleType#]
% Cmax.Spikes:    the ensemble-types identified by the maximum modularity (Q) clustering from the community
%                 detection algorithm (format: [EnsembleType#]
% AcorrTypes:     the ensemble-types idenitified by significant peaks or
%                 troughs in their autocorrelogram
% 
% Mark Humphries 16/6/2014
clear all; close all

fname_Pre = []; % [];
stimset = 'da01';  % 'da01': first; 'da02': second; and 'da03' : third ("rest" / control)
stimstart = 30; % 30 s into recording

addpath ../../Functions/

% load ensemble statistics
fname = [fname_Pre stimset '_Analyses_Neurons_and_Groups'];
load(fname)


M = 5; % marker size
flag = '3'; % plotting flag

dists = {'sqEuclidean'};  % options for all-eigenvector consensus
rpts = 100;

%% classify by auto-correlogram

AcorrTypes = zeros(length(groupdata),1);

% four possible outcomes: 
% (1) no significant peaks or troughs = "non-oscillatory"
AcorrTypes([groupdata.thisP] == 0 & [groupdata.thisN] == 0) = 1;  

% (2) significant peaks and troughs = "oscillator"
AcorrTypes([groupdata.thisP] == 1 & [groupdata.thisN] == 1) = 2;  

% (3) significant peaks but no significant troughs = "burster"
AcorrTypes([groupdata.thisP] == 1 & [groupdata.thisN] == 0) = 3;  

% (4) no significant peaks but significant troughs = "pauser"
AcorrTypes([groupdata.thisP] == 0 & [groupdata.thisN] == 1) = 4;  


%% spike-train metrics: fits to distributions

% concatenate vectors of model-fits: make fit-space
Spikes.dataset= [[groupdata.pAICs_isi]' [groupdata.pAICs_cv2s]'];

Y = pdist(Spikes.dataset,'euclidean'); % find distance between each pair of ensembles in fit-space
DxyData = squareform(Y);
nnodes = size(DxyData,1)

% modularity clustering: using similarity
Spikes.SxyData = exp(-DxyData.^2);  % exponential conversion to similarity
Spikes.SxyData(eye(nnodes)==1) = 0; % zeros on diagonal

% cluster using consensus community detection

tic
[Spikes.Cmax,Spikes.Qmax,Spikes.Ccon,Spikes.Qc,N,Qvec] = allevsplitConTransitive(Spikes.SxyData,dists,rpts);
toc

% Spikes.Ccon IDs ensemble-type of each ensemble across recordings

%% save stuff
save([fname_Pre stimset '_Ensemble_Types'],'Spikes','dists','rpts','AcorrTypes')


