function [c] = coef(T, scale)
% Calculates coefficient c at different resolutions
% Divide scale by 10 to get degrees
% Estimates coefficient using OLS estimator, 
% then sets cutoff as 90th percentile of OLS residuals
% then re-estimates c by minimizing Huber's loss

    tbar = block_avg_scale(T, scale) ; 
    sig_TT  = block_single_sd_scale( T, scale) ;
    lgt = len_grad(tbar) ;
    
    
    nonan = ~isnan(sig_TT(:)+lgt(:)) ;
    x = lgt(nonan) ;
    y = sig_TT(nonan) ;

    % OLS: forcing intercept to be zero
    b = (x'* x)\ (x'*y) ;
    err = abs( y - ( b.* x ) ) ;

    % set threshold between quadratic and linear loss to be
    % 90th percentile of OLS errors
    alpha = 0.9 ;
    delta = prctile(err, 100*alpha) ;

    betas = 0.08:.001:0.28 ;
    huber = zeros(size(betas)) ;
    for i = 1:length(betas)
        huber(i) = huber_loss(x, y, betas(i), delta) ;
    end

    [d, Iall] = min(huber) ; 
    c = betas(Iall) ;

end