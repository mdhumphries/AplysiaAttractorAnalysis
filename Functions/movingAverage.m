function y = movingAverage(x, w)
   % moving average using fixed window, applied to every point in signal 
   % x = data; w = window-size (must be odd)
   % y is size of original signal; plot from (w+1)/2 for first full window
   % on signal
   if ~rem(w,2) error('window size must be odd'); end
   k = ones(1, w) / w;    
   y = conv(x, k);
   trim = (w-1)/2;
   y = y(trim+1:end-trim);
end