function sd11  = block_single_sd_scale( mat1, scale )
% input: mat1 is high resolution
% output: sd11 is low resolution
% average over blocks of size scale x scale

  avrg1 = block_avg_scale(mat1, scale) ;
  big_avrg1 = block_expand_scale(avrg1, scale) ;
  anomaly1 = (mat1 - big_avrg1) ;  

  sd11 = block_avg_scale(anomaly1.^2, scale) ;

end
