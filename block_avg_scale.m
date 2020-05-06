function avrg = block_avg_scale( mat , scale)
  % input is a 3d matrix mat to be smoothed and a scalar scale. In mat land has value NaN.
  % output is a 3d matrix smoothed horizontally in scale x scale blocks

  nrow = scale * floor( size(mat, 1) / scale )  ;
  ncol = scale * floor( size(mat, 2) / scale ) ;
  ndeep = size(mat, 3) ;

  nrow_avg = floor(nrow / scale) ;
  ncol_avg = floor(ncol / scale) ;
  avrg = zeros(nrow_avg, ncol_avg, ndeep) ;

  i = 0 ;
  for start_row = 1:scale:nrow
    i = i+1 ; 
    end_row = start_row + scale - 1 ; 
    j = 0 ;
    for start_col = 1:scale:ncol
      j = j+1 ;
      end_col = start_col + scale - 1 ;
      ocean_mass = sum( sum( mat(start_row:end_row, start_col:end_col, :) ) ) ;
      avrg(i, j, :) = ocean_mass ./ ( scale^2 ) ;
    end
  end

end
