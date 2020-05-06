function [out, delta] = huber_loss(x, y, beta, delta)
% Loss function used in robust regression.
% source: http://docs.h2o.ai/h2o/latest-stable/h2o-docs/data-science/algo-params/huber_alpha.html
% The Huber loss function is a combination of the squared-error loss function and absolute-error loss function. 
% It applies the squared-error loss for small deviations from the actual response value and the absolute-error loss for large 
% deviations from the actual respone value. The alpha parameter dictates the threshold between quadratic and linear loss 
% (i.e. the top percentile of error that should be considered as outliers). 
% This value must be between 0 and 1 and defaults to 0.9.
%
% Function:
%
% err = y - beta x
% L_d (err) = 1/2 err^2          if |err| < d
%             d(|err| - 1/2 d)   else
%
% Source for function: Wikipedia Huber Loss

  
  alpha = 0.9 ;
  err = abs( y - ( beta .* x ) ) ;

  % unless delta is specified, default to 90th percentile of residuals
  if ~exist('delta', 'var')
    delta = prctile(err, 100*alpha) ;
  end
  
  small = abs(err) < delta ;
  big = ~small ;
  err(small) = 0.5 .* err(small) .^2 ;
  err(big) = delta .* ( err(big) - (0.5 .* delta) ) ;
  out = sum(err) ;

end
