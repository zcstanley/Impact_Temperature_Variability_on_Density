function [sd11 sd12 sd22 ] = block_sd_scale( mat1, mat2, scale)

  nrow = scale * floor( size(mat1, 1) / scale ) ;
  ncol = scale * floor( size(mat1, 2) / scale ) ;
  ndeep = size(mat1, 3) ;

  avrg1 = block_avg_scale(mat1, scale) ;
  big_avrg1 = block_expand_scale(avrg1, scale) ;
  anomaly1 = (mat1 - big_avrg1 ) ;  

  avrg2 = block_avg_scale(mat2, scale) ;
  big_avrg2 = block_expand_scale(avrg2, scale) ;
  anomaly2 = (mat2 - big_avrg2 ) ;

  sd11 = block_avg_scale(anomaly1.^2, scale) ;
  sd22 = block_avg_scale(anomaly2.^2, scale) ;
  sd12 = block_avg_scale(anomaly1.*anomaly2, scale) ;
 
end
