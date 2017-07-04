function [spkfcn,idxs] = convolve_spiketrains(spkdata,h,shiftbase,Didxs,bins,bin,T,window)
    Nidxs = numel(Didxs);
    
    %% go round and compute spike-train binless functions
    spkfcn = zeros(numel(bins),Nidxs);
    nspikes = NaN; % just in case there are no spikes....
    
    for j = 1:Nidxs
        currix = find(spkdata(:,1) == Didxs(j));
        nspikes(j) = numel(currix);
        [spk,bts] = spike_train_from_times(spkdata(currix,2),bin,T);
        if nspikes(j) > 0   % only bother doing convolution if there's something to convolve!!
            switch window
                case 'Gaussian'
                    try
                        y = conv(h,spk);
                        [r c] = size(y); if c>1 y = y'; end  % for reasons best known to Matlab, certain convolutions will return this as a row vector rather than a column vector
                        shifty = y(shiftbase+1:end-shiftbase);   % shift convolved signal to line up with spike-times
                        if numel(shifty) < numel(bins)
                            % pad with zeros
                            diffbins = numel(bins) - numel(shifty);
                            shifty = [zeros(diffbins,1); shifty]; 
                        end % can occasionally happen with width pars that are not integer multiples of step-size 
                        spkfcn(:,j) = shifty;
                    catch
                        disp('I ran into a problem convolving the Gaussian')
                        keyboard
                    end
                    
                case 'exponential'
                    y = conv(h,spk);
                    [r c] = size(y); if c>1 y = y'; end  % for reasons best known to Matlab, certain convolutions will return this as a row vector rather than a column vector
                    % keyboard
                    spkfcn(:,j) = y(1:numel(bins));   % truncate convolved signal to line up with recording time   
            end
        end
        % keyboard
    end
   
    % keyboard
    
    try
    % if not firing, then strip out cells
    spkfcn(:,nspikes==0) = []; 
    idxs = Didxs; idxs(nspikes==0) = [];
    catch
        disp('I ran into a problem removing non-firing cells')
        keyboard
    end
    
    % convert to firing rate (spikes/s)
    spkfcn = spkfcn ./ bin; % sums to 1 if spike in every bin
