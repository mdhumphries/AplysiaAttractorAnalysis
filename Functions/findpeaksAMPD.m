function p = findpeaksAMPD(x,varargin)

% FINDPEAKSAMPD uses the AMPD algorithm to find all peaks in a time-series
% P = findspeaksAMPD(T) given the single time-series in array T, returns all identified
% peaks in array P. Assumes that T is uniformly sampled.
%
% ... = findpeaksAMPD(..,'r') will detrend the time-series using linear
% regression before applying AMPD
%
% NOTES:
% The automatic multiscale-based peak detection (AMPD) algorithm is parameter-free. 
% 
% References:
% 
% Scholkmann, F.; Boss, J. & Wolf, M. "An Efficient Algorithm for Automatic 
% Peak Detection in Noisy Periodic and Quasi-Periodic Signals". Algorithms, 2012, 5, 588-603 
%
% Mark Humphries 24/10/2014

% offset constant
alpha = 1;

blnR = 0;
if nargin >= 2
    if findstr(varargin{1},'r')
        blnR = 1;
    end
end

% detrend x
if size(x,2) > 1
    x = x';  % make into column vector
end

N = numel(x);
if blnR
    y = [1:N]';
    [b,bint] = regress(x,[ones(N,1) y]);
    x = x - (b(1) + b(2)*y);
end

L = ceil(N/2) - 1; % maximum window size

k = 1:L;
wsize = 2 * k;

m = zeros(numel(k),N);

% for each scale, check for local peaks
for iK = 1:numel(k)
    m(iK,1:k(iK)) = 1; %rand(1,k(iK)) + alpha;
    m(iK,N-k(iK)+1:N) = 1; %rand(1,k(iK)) + alpha;
    for i = k(iK)+1:N-k(iK)
        % check if current value (x (i)) is greater than values k steps
        % away ahead/behind
        if ~(x(i) > x(i-k(iK)) & x(i) > x(i+k(iK)))  % if not, then set random number
            m(iK,i) = 1; % rand + alpha; 
        end
    end
end

gamma = sum(m,2);  % row-wise sum of elements
ixMin = find(gamma == min(gamma)); % minimum = scale with greatest number of local peaks

mR = m(1:ixMin,:); % remove all rows below it
sumMR = sum(mR); % column-wise sum: if 0 then agrees at all scales

p = find(sumMR == 0);  % indices of peaks


