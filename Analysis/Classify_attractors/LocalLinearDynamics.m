function [Emax,Egs,RotPeriod,Npts,Cond,R2] = LocalLinearDynamics(IDs,statespace,ix0,Theta,NeighbourD,DT)

% LOCALLINEARDYNAMICS compute local dynamics around point in state space
% [Emax,R,N,E,C,R2] = LOCALLINEARDYNAMICS(ID,S,I,T,D,DT) computes the linear system in 
% the neighbourhood of every point in ID within the N-dimensional state-space S 
% (a MxN matrix, of M time-points). The index I is the offset of the IDs
% from index 1 of the state-space (i.e. IDs=1 is 1+I in the state-space).
% The linear model is fit to all state-space points in the neighbourhood of
% +/- T*D, where T is the threshold used to find recurrence points, and D
% is a scalar. DT is the time-step of the state-space (in seconds)
%
% Returns: Emax, the array of maximum eigenvalues; R, the array of rotation periods; 
% N, array of the number of points used for the regression to fit the linear model;
% and E, all eigenvalues fo every tested point in a MxN matrix. Extra
% information returned for every points are: C, the condition number for the matrix of
% coefficients; R2, the linear correlation between the model and data
% trajectories (computed from first data-point; using concatenated
% time-series).
%
% Reference:
% Lathrop, D. P. & Kostelich, E. J. Characterization of an experimental 
% strange attractor by periodic orbits. Phys Rev A, 1989, 40, 4028-4031
% 
% ChangeLog:
% 24/3/17: ensure all eigenvalues are sorted in magnitude order (high to
% low)
% 
% Mark Humphries 16/11/2015

dims = size(statespace,2); 

Egs = zeros(numel(IDs),dims);

% feature('numCores');
if isempty(gcp('nocreate'))
    parpool('local',10);
end

% parfor iN = 1:numel(IDs)  % for each recurrence point....
for iN = 1:numel(IDs)
   % get points within orbital neighbourhood
    ixRecr = IDs(iN);
    pt0 = statespace(ix0+ixRecr,:);  % x_i of recurrence point (ix0 is first sampled point in state-space)
    D = sqrt(sum(bsxfun(@minus,statespace,pt0).^2,2)); % Euclidean distance to all points
    currixs = find(D <= NeighbourD * Theta);  % all distances within Theta neighbourhood
    
    back = 0; fwd=numel(currixs+1);
    
    % reduce to set of points contiguous with the chosen point
    ixPt = find(currixs == ix0+ixRecr);  % recurrence point's location within the neighbourhood
    ixshift = abs(currixs - ixPt);  
    dshift = diff(ixshift); % difference between adjacent indices in neighbourhood
    back = find(dshift(1:ixPt-1) > 1,1,'last'); if isempty(back) back = 0; end
    fwd = find(dshift(ixPt:end) > 1,1,'first')+ixPt; if isempty(fwd) fwd = numel(currixs)+1; end
    
    
    % keyboard
    Npts(iN) = numel(currixs(back+1:fwd-1));
    ixPtshft = ixPt - back+1; % index of current point within retained contiguous portion
    ptfit = statespace(currixs(back+1:fwd-1),:);

   % fit linear model
   try
   coeffs = ptfit(1:end-1,:)\diff(ptfit); % on dx/dt  = Ax(t)
   catch
       keyboard
   end
   % coeffs = ptfit(1:end-1,:)\ptfit(2:end,:); % on map xt+1 = Axt
   % coeffs = mvregress([ones(Npts(iN)-1,1),ptfit(1:end-1,:)],ptfit(2:end,:)); % include Ax+b model
   % coeffs = mvregress(ptfit(1:end-1,:),ptfit(2:end,:));

   % get eigenvalues
   Egs(iN,:) = sort(eig(coeffs),'descend');
   Emax(iN) = max(Egs(iN,:));  % keep full complex number -take ABS later
   Cond(iN) = cond(coeffs);  % condition number of this matrix
   
   
   % local linear model prediction
   modelpred = zeros(size(ptfit));
   
   % forward
   modelpred(ixPtshft,:) = ptfit(ixPtshft,:);
   for i=ixPtshft:Npts(iN)-1
        modelpred(i+1,:) = [modelpred(i,:)' + coeffs * modelpred(i,:)']'; % dx/dt
        % modelpred(i+1,:) = [coeffs * modelpred(i,:)']'; % map
   end
   
   % backward
   for i=ixPtshft-1:-1:1
        modelpred(i,:) = [modelpred(i+1,:)' + coeffs * modelpred(i+1,:)']'; % dx/dt
        % modelpred(i+1,:) = [coeffs * modelpred(i,:)']'; % map
   end   

   r = corr(modelpred(:),ptfit(:)); 
   R2(iN) = r^2; 

%    
%    figure
%    plot(ptfit(:,1),ptfit(:,2)); hold on
%    plot(modelpred(:,1),modelpred(:,2),'r')

end

% for this period,                
% compute orbital frequency
RotPeriod = 2*pi ./ imag(Emax) * DT;  % Strogatz, pg 134: local estimate assuming no contraction or expansion
RotPeriod(isinf(RotPeriod)) = 0;  % inf = stable (no rotation locally) 
