function r2 = r_squared(x, y)
% Calculates R^2 where y is a vector (or matrix) of observations and x is the predicted values
% R^2 = 1 - total sum of squares / residual sum of squares

    these = ~isnan(x+y) ;
    s_tot = sum( ( y(these) - mean(y(these)) ) .^2 ) ; % total sum of squares
    s_resid = sum( (y(these) - x(these)).^2 ) ; % residual sum of squares
    r2 = 1 - s_resid./s_tot ;

end