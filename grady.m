function grady = grady(field, dy)
% computes gradient using Euler method
  
  nz = size(field, 3) ;

  temp_up = field ;
  temp_down = circshift(field, [-1 0 0]) ;
  dy_recip = repmat( 1./dy, 1, 1, nz) ;

  grady = dy_recip .* (temp_up - temp_down) ;
  % TO DO: set last row to nan

end

