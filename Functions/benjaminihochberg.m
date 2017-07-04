function [H,T] = benjaminihochberg(P,alpha)

% BENJAMANIHOCHBERG do Benjamini-Hochberg procedure for multiple comparisons
%   [H,T] = BENJAMANIHOCHBERG(P,ALPHA) applies the Benjamini-Hochberg procedure for multiple
%   comparison tests to the p-values in array P, given an overall required
%   significance level of ALPHA. Returns  an array of hypothesis test 
%   values H (entry = 0  if null hypothesis not rejected, otherwise entry =
%   1); and also returns the rejection threshold T.
%
%   Notes: the Benjamini-Hochberg procedure controls the False Discovery
%   Rate, not the Family-Wise Error rate; it is thus less conservative and
%   rejects more null hypotheses (cf Shaffer, 1995). FWE: minimise
%   rejection of true null hypotheses; FDR: minimise expected proportion of
%   all erroneous rejections amongst all rejections.
%
%   TO DO:
%   (i) Check is meaningful to derive adjusted p-values
%   (ii) implement False Discovery Rate control for dependent p-values from Benjamini & Yekuteili
%    (2001) Annals of Stats: very similar to original Benjamini-Hochberg, but more
%     robust...
%
%   References:
%   Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery rate
%   in multiple testing under dependency The Annals of Statistics, 29, 1165-1188
%
%   Shaffer, J. P. (1995) "Multiple hypothesis testing" Annu. Rev.
%   Psychol., 46, 561-584.
%
%   Wasserman,L. (2002) All of Statistics. Springer: Berlin.
%
%   Wright, S. P. (1992) "Adjusted p-values for simultaneous inference",
%   Biometrics, 48, 1005-1013.
%   
%   Mark Humphries 1/4/2012

[r c] = size(P);
if c == 1
    P = P';
end

num_pairs = numel(P);

[Ps,Pi] = sort(P);

% 1. Benjamini-Hochberg: from Wasserman
% assuming p-values are independent
L = alpha.*(1:num_pairs) ./ num_pairs;
R = find(Ps < L,1,'last');    % find first P-value to exceed corrected FDR threshold
reject = 1:R;

if isempty(reject)
    % no null hypothesis rejected, so no threshold P
    T = nan;
else
    T = Ps(R);
end

H = zeros(num_pairs,1);
H(Pi(reject)) = 1;

%keyboard

% [temp,idx] = sort(Pi);
% Padj = Padj(idx)';