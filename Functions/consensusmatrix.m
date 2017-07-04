function Sc = consensusmatrix(Grps)
    % Pass: NxC matrix of C clusterings of N objects  
    [nIDs nreps] = size(Grps);
    Sc = zeros(nIDs);
    
    % pair-wise similarity matrix
    for nr = 1:nIDs
        for nC = nr:nIDs
            if nr ~= nC
                Gi = Grps(nr,:); Gj = Grps(nC,:);
                Sc(nr,nC) = sum(Gi == Gj) / nreps;
                Sc(nC,nr) = Sc(nr,nC);
            end
        end
    end