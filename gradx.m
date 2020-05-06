function gradx = gradx(field, dx)
% computes gradient using Euler method
% field is periodic in x-direction
  
  nz = size(field, 3) ;

  temp_left = field ;
  temp_right = circshift(field, [0, -1, 0]) ;
  dx_recip = repmat( 1./dx, 1, 1, nz) ;

  gradx = dx_recip .* (temp_right - temp_left) ;

end

