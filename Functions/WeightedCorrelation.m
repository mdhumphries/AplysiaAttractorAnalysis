function [r,p] = weightedcorrelation(X,Y,W,type)

% WEIGHTEDCORRELATION compute weighted correlation coefficients
% [R,P] = WEIGHTEDCORRELATION(X,Y,W,TYPE) computes the weighted correlation
% coefficient between the data in arrays X and Y, using weights W (all
% N-element arrays). Choose TYPE = {'Pearson','Spearman'}
%
% Returns:
%       R: the weighted correlation coefficient (of TYPE)
%       P: approximate P-value for R (permutation test: limit is p=0.0001: 1e-4)
%
% To Do:
% (i) generalise to matrices (X,Y,W) using repmat and bsxfun
% (ii) set permutation test options
%
% Mark Humphries 6/7/2016

nP = 1e5;

% nP = 10;

% strip NaNs
ixnan = isnan(X) | isnan(Y);
X = X(~ixnan);
Y = Y(~ixnan);
W = W(~ixnan);

% keyboard

n = numel(X);

% convert to ranks?
ranks = 1:n;
Xrank = zeros(n,1); Yrank = Xrank;
if findstr(type,'Spearman')
    % convert to ranks
    [~,I] = sort(X,'ascend');
    Xrank(I) = ranks;  
    [~,I] = sort(Y,'ascend');
    Yrank(I) = ranks;
    r = computeR(Xrank,Yrank,W);  % data correlation
else
    r = computeR(X,Y,W);  % data correlation

end


% get p-value from permutation test
rperm = zeros(nP,1);
pool = [X; Y];
for iP = 1:nP
    % permute
    ixR = randperm(2*n);
    Xp = pool(ixR(1:n)); 
    Yp = pool(ixR(n+1:end));
    if findstr(type,'Spearman')
        % convert to ranks
        [~,I] = sort(Xp,'ascend');
        Xrank(I) = ranks;  
        [~,I] = sort(Yp,'ascend');
        Yrank(I) = ranks;
        rperm(iP) = computeR(Xrank,Yrank,W);  % correlation
    else
        rperm(iP) = computeR(Xp,Yp,W);  %correlation

    end
    
end

if r > 0
    p = sum(rperm > r) ./ nP;
else
    p = sum(rperm < r) ./ nP;
end

% function to compute r
function r = computeR(X,Y,W)
sW = sum(W);

mX = sum(W.*X)/sW;
mY = sum(W.*Y)/sW;

covX = sum(W.*(X-mX).*(X-mX)) / sW;
covY = sum(W.*(Y-mY).*(Y-mY)) / sW;
covXY = sum(W.*(X-mX).*(Y-mY)) / sW;

r = covXY ./ sqrt(covX * covY);

