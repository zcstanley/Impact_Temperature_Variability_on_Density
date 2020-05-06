function [rho_TT rho_ST rho_SS] = dens_wright_sot(T, S, p)

 % Second order correction to Wright equation of state for sea water
 % Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
 % Units: T[degC], S[PSU], p[Pa]
 % Returns density per temperature^2 [kg m^-3 degC^-2]

  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7 ;
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4; b3 = 2.084372e2;  b4 = 5.944068e5; b5 = -9.643486e3;
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422;   c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464;

 % Zero-th order terms
  al0 = a0 + a1.*T + a2.*S ;
  p0  = b0 + b4.*S + T .* (b1 + T.*(b2 + b3.*T) + b5.*S) ;
  l = c0 + c4.*S + T .* (c1 + T.*(c2 + c3.*T) + c5.*S) ;
  d0 = l + al0.*(p + p0) ;
  d02 = d0.^2 ;
  d03 = d0 .^3 ;
   
 % First order terms
  p0_T = b1 + b5.*S + T .* (2*b2 + 3*b3.*T) ;
  p0_S = b4 + b5.*T ;
  l_T = c1 + c5.*S + T .* (2*c2 + 3*c3.*T) ;
  l_S = c4 + c5.*T ;
  d0_T = l_T + a1.*(p + p0) + al0.*p0_T ;
  d0_S = l_S + a2.*(p + p0) + al0.*p0_S ;

 % Second order terms
  p0_TT = 2*b2 + 6*b3.*T ;
  l_TT = 2*c2 + 6*c3.*T ;
  d0_TT = l_TT + 2*a1.*p0_T + al0.*p0_TT ;
  d0_ST = c5 + a1.*p0_S + a2.*p0_T + b5.*al0 ;
  d0_SS = 2*a2.*p0_S ;

 % Density corrections
 rho_TT = ( 2.*(p+p0).*d0_T.^2 - 2.*d0.*d0_T.*p0_T - (p+p0).*d0.*d0_TT + d02.*p0_TT ) ./ d03 ;
 rho_ST = ( 2.*(p+p0).*d0_T.*d0_S - d0.*d0_T.*p0_S - d0.*d0_S.*p0_T - (p+p0).*d0.*d0_ST + b5.*d02 ) ./ d03 ;
 rho_SS = ( 2.*(p+p0).*d0_S.^2 - 2.*d0.*d0_S.*p0_S - (p+p0).*d0.*d0_SS ) ./ d03 ;

end

