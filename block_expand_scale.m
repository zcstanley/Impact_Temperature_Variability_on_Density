function out_mat = block_expand_scale( mat , scale)

  nrow = size(mat, 1) ;
  ncol = size(mat, 2) ;
  ndeep = size(mat, 3) ;
  out_mat = zeros( scale*nrow, scale*ncol, ndeep ) ;

  for i = 1:nrow
    start_row = scale * (i - 1) +  1 ;
    end_row =  scale * i ;
    for j = 1:ncol
      start_col = scale * (j - 1) + 1 ;
      end_col = scale * j ;
      out_mat( start_row:end_row, start_col:end_col, : ) = repmat(mat(i, j, :), [scale, scale, 1]) ;
    end
  end

end
