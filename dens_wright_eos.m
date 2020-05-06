function rho = dens_wright_eos(T, S, p)

 % Equation of state for sea water given by Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
 % Units: T[degC],S[PSU],p[Pa]
 % Returns density [kg m-3]

  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7 ;
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4; b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3;
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422;   c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464;

  al0 = a0 + a1.*T + a2.*S;
  p0  = b0 + b4.*S + T .* (b1 + T.*(b2 + b3.*T) + b5.*S);
  l = c0 + c4.*S + T .* (c1 + T.*(c2 + c3.*T) + c5.*S);
  rho = (p + p0) ./ (l + al0.*(p+p0));
   
end
