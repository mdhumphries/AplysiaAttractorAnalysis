%% colormap and linewidths for every figure
% default parametersL overwrite locally if needed

format = 'png'; %'eps' is possible, PNG for files EPS struggles with
color = 'rgb'; % for PNG only
dpi = 1200; % for PNG only

fontsize = 6;
plotlinewidth = 1;
pclinewidth = 0.75;
errorlinewidth = 0.75;
axlinewidth = 0.5;
M = 2;  % default marker size

% color code of stimulus bars etc
stimcolor = [0.8 0.3 0.3];

% colour code of PC projections
sponcolor = [0.5 0.5 0.5];

stim_set = [0.3 0.3 0.3;... 
            0.6 0.6 0.6;...
            0.8 0.8 0.8];

        % pertubations:
jumpcolor = [0.8 0.5 0];
endcolor = [0.2 0.2 1];

% color-code of every preparation
prep_Cmap = brewermap(10,'Paired'); 

% colour-code of ensemble types
Typeclrs = {'Purples','Greys','Blues','Reds'};
