% script to count diagonal and vertical lines in all recurrence plots

clear all; close all

addpath ../../Functions/

fname = 'da03';
prefix = [];

% get PCA projections
load(['../' prefix fname '_DataProperties_FunctionAndWindowSize'],'PCAdata','Data','SDF')
nfiles = numel(PCAdata);

% parameters
stimstart = 30; % (s)
T0 = 35;  % start counting after stim (s)
maxWindow = 10;  % stop before end of time-series (s)
VarExplained = 0.80;    % how many PCs to keep

prctTheta = [2,5,10]; %  % of distance distribution to use as threshold; max distance to check (spikes/s)
% stepTheta = 0.5;
% Theta = stepTheta:stepTheta:maxTheta;

bstep = 1;  % s
bins = 0:bstep:50;  % recurrence time bins (s);

tmin = 5;  % s; minimum recurrence time to include in diagonal line counting

NeighbourD = 6;  % multiplier of threshold to use as neighbourhood around (m,e) point
minNeighbours = 50; % minimum number of points to regression

minbin = bins(find(bins <= tmin,1,'last')); % lowest retained bin in histogram

%% find attractor points
for iR = 1:nfiles
    iR

    % pick the subset of dimensions for whole thing 
    RecurData(iR).maxD = sum(PCAdata(iR).prctVar <= VarExplained);
    currPCs = PCAdata(iR).PCproj(:,1:RecurData(iR).maxD);  % from stimulation
    
    % ALSO NEED TO DO THIS INCLUDING SPONTANEOUS STATE.... 
    % make recurrence plot and related
    dPs = pdist(currPCs,'euclidean');
    RecurData(iR).Theta = prctile(dPs,prctTheta);
    
    % keyboard
    
    N = size(currPCs,1);    % number of time-points
    k = round(tmin / Data(iR).GaussQt); % number of bins with local dynamics only

    for iT=1:numel(RecurData(iR).Theta)
        iT
        % find recurrence points
        RecurData(iR).RecurPlot(iT).Rp = sparse(squareform(dPs <= RecurData(iR).Theta(iT)));  % automatically puts 0s on diagonal, so no LOI

        % recurrence rate 
        RecurData(iR).RecurPlot(iT).RR = sum(sum(RecurData(iR).RecurPlot(iT).Rp)) / N^2;
        
        %% get all diagonal lines above minimum recurrence time k
        
        % debug: get small submatrix 
%         N = 1000;
%         submatrix = RecurPlot(iT).Rp(2000:2000+N-1,4500:4500+N-1); % make submatrix; fix look up below when getting L...
%         subs = find(triu(submatrix));  % all subscripts of non-zero entries; plot is symmetric, so just need upper tri
%         checkmatrix = triu(submatrix);
        
        % fast counting using checkmatrix
        RecurData(iR).RecurPlot(iT).dlengths = [];
        checkmatrix = triu(RecurData(iR).RecurPlot(iT).Rp,k);
        %tic
        left = sum(sum(checkmatrix));
        Cctr = 1;  % column ctr
        while left > 0
            S = sum(checkmatrix(:,Cctr)); 
            if S
                seed = find(checkmatrix(:,Cctr)==1); % all first points of diagonals remaining in this column of checkmatrix
                maxDiag = N - Cctr  + 1;

                for iS = 1:numel(seed)  % loop over all possible i's
                    % get diagonal line
                    allD = [seed(iS)+(1:maxDiag)'-1   Cctr+(1:maxDiag)'-1];  % all +k points
                    L = find(checkmatrix(sub2ind([N,N],allD(:,1),allD(:,2)))==0,1,'first')-1; %^debugging
                    if isempty(L) L = maxDiag; end  % if empty, then diagonal is to the edge of the plot
                    RecurData(iR).RecurPlot(iT).dlengths = [RecurData(iR).RecurPlot(iT).dlengths; L];  % store line length
                 
                    % set checkmatrix entries to zeros
                    del = sub2ind([N,N],allD(1:L,1),allD(1:L,2));  % set of indices to delete         
                    checkmatrix(del) = 0;
                    left = left - L;  % remove from "left"
                    
                end
                    
            end
            Cctr = Cctr+1;  % check next column
        end
        %toc
        
%         tic
%        subs = find(triu(RecurPlot(iT).Rp,k));  % all subscripts of non-zero entries; plot is symmetric, so just need upper tri
%
%         RecurPlot(iT).dlengths = [];
%         while ~isempty(subs)
%             % convert first in queue to subscripts
%             [i,j] = ind2sub([N,N],subs(1)); 
% 
%             maxDiag = N - j + 1;  % max length of diagonal from the current starting position
% %             blnLine = 1; L = 1;
% %             while blnLine
% %                 ix = intersect(find(i(1)+L == i), find(j(1)+L == j)); 
% %                 if ix
% %                     i(ix) = []; % remove next point on the line from queue
% %                     j(ix) = [];
% %                     L = L+1
% %                 else 
% %                     blnLine = 0;
% %                 end                
% %             end  % finish searching for line
% %             i(1) = []; j(1) = [];  % remove start point
% %             RecurPlot(iT).dlengths = [RecurPlot(iT).dlengths; L];  % store line length
%            
%             allD = [i+(1:maxDiag)'-1   j+(1:maxDiag)'-1];  % all +k points
%             % L = find(RecurPlot(iT).Rp(sub2ind([N,N],allD(:,1),allD(:,2)))==0,1,'first')-1;  % length is entries until first 0
%             L = find(submatrix(sub2ind([N,N],allD(:,1),allD(:,2)))==0,1,'first')-1; %^debugging
%             if isempty(L) L = maxDiag; end  % if empty, then diagonal is to the edge of the plot
%             RecurPlot(iT).dlengths = [RecurPlot(iT).dlengths; L];  % store line length
%             % indices to delete must be the first entry, and then
%             
%             del = sub2ind([N,N],allD(1:L,1),allD(1:L,2));  % set of indices to delete         
%             ixdel = arrayfun(@(x) find(subs==x),del);
%             subs(ixdel) = []; % slow for huge set of indices
%         end
%         toc
        
        %% get all vertical lines 
        sumC = sum(RecurData(iR).RecurPlot(iT).Rp);
        % sumC = sum(submatrix);
        RecurData(iR).RecurPlot(iT).Vlengths = [];
        for iC = 1:N  % each column
            if sumC(iC)
                %V = find(submatrix(:,iC)==1); 
                V = find(RecurData(iR).RecurPlot(iT).Rp(:,iC)==1); % find indices of all possible vertical lines entries
                dV = [diff(V); inf];   % get change in index
                ix = [0; find(dV > 1)]; % line ends are discontinuous changes in index
                Vs = ix(2:end) - ix(1:end-1); % line lengths
                RecurData(iR).RecurPlot(iT).Vlengths = [RecurData(iR).RecurPlot(iT).Vlengths; Vs];  % store line length
            end
        end

        
    end
%     figure
%     imagesc(RecurData(1).RecurPlot(3).Rp); % (2000:2000+N,4500:4500+N));
%     colormap(gray)
%     % keyboard
   
   

end  % over all files

save([fname '_RecurrencePlotData'],'RecurData'); % ,'-v7.3');
