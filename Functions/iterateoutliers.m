function [Out,Pars,LogL,varargout] = iterateoutliers(data,model)

% ITERATEOUTLIERS find outliers by maximising likelihood fit of model
% [O,P,L] = ITERATEOUTLIERS(D,M) given the vector of data D, finds the 
%  potential outliers according to model M: all unimodal models available from
%  function 'FITDIST' 
%
% Returns: O, the array of indices into D of identified outliers; P, the
% final fitted model as a probability distribution object; and L the log likelihood as a
% function of removed data points
%
% Algorithm:
% (1) Fit unimodal model M to D using max likelihood: L_D
% (2) Remove furtherest data-point from estimated location parameter, giving D*
% (3) Fit new model M to D*, giving likelihood L_D*
% (4) Repeat from (2) until maximised L_D*
%
% Mark Humphries 25/11/2015

% switch model
%     case 'Normal'
%     case 'Lognormal'
% end

checkstop = inf;

pd1 = fitdist(data,model);
ll(1) = -pd1.negloglik;
ixD = ones(numel(data),1);
ixSeq = 1:numel(data);

% keyboard

blnStop = 0; ctr = 1; checkctr = 0;
while ~blnStop | checkctr <= checkstop
    dist = abs(data(ixD==1)-pd1.mu);   % distances of remaining data-points
    ixOmit = find(dist == max(dist));  % find max distance of remaining data-points 
    ixD(ixSeq(ixOmit)) = 0;  % now omit this data-point from overall index
    ixSeq(ixOmit) = []; % running vector of remaining indexes into data array
    try
        pd2 = fitdist(data(ixD==1),model);
    catch
        keyboard
    end
    ll(ctr+1) = -pd2.negloglik;
    if ll(ctr+1) < ll(ctr) 
        if blnStop == 0
            ixD(ixSeq(ixOmit)) = 1; % put data-point back in....
            Out = find(ixD == 0); % save result as this point
            Pars = pd1;
            LogL = ll(1:ctr+1);
            blnStop = 1;
            checkctr = 1;
            checkstop = min((numel(data) - numel(Out)) / 2,50);  % at most check half of remaining data or 50 more points
        else
            checkctr = checkctr+1;
            % run another few to check this is true maximum
        end
    end
    pd1 = pd2; % iterate models
    ctr=ctr+1;

end

varargout{1} = ll;
