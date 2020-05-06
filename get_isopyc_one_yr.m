% save delta rho and R^2 for temperature term as a model for delta rho
% data is isopycnal

delta_rho = zeros(240, 360, 15, 73);
r2_temp_v_delta_rho = zeros(73, 1);

% scale
scale = 10 ;

% load grid
load('/glade/work/zofias/Brankart/Data/dxt_dyt.mat')
% dxt: x-distance between midpoints of T cells
% dyt: y-distance between midpoints of T cells

% coarse grid
dxt(dxt<0) = NaN ;
dyt(dyt<0) = NaN ;
cdxt = scale .* block_avg_scale(dxt, scale) ; 
cdyt = scale .* block_avg_scale(dyt, scale) ;

% loop over one year
years = ["0013"];
months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"];
days =  ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", ...
	     "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31" ];
jj = 0 ;
for Y = years
for M = months
for D = days

  date_str = strcat(Y, '-', M, '-', D) ;
  filename = strcat('/glade/scratch/zofias/isopyc_datag.e01.GIAF.T62_t12.003.pop.iso0.', date_str, '.nc') ;

  if isfile(filename)
  jj = jj + 1;

  fprintf('Iteration %g. Date: %s\n', jj, date_str)

  % open netcdf file
  ncid = netcdf.open(filename, 'NC_NOWRITE');

  % load variables from netcdf
  ZSIG = netcdf.getVar(ncid, 2);
  TEMP = netcdf.getVar(ncid, 4);
  SALT = netcdf.getVar(ncid, 5);

  % rotate variables for ease of plotting
  Z = zeros(2400, 3600, 15);
  T = Z ;
  S = Z ;
  for ii = 1:15
    Z(:,:,ii) = rot90(ZSIG(:,:,ii)) ;
    T(:,:,ii) = rot90(TEMP(:,:,ii)) ;
    S(:,:,ii) = rot90(SALT(:,:,ii)) ;
  end

  % change from cm to m
  Z = Z./100 ;

  % set land and missing values to NaN
  Z(Z>1e20) = NaN ;
  T(T>1e20) = NaN ;
  S(S>1e20) = NaN ;

  % clear netcdf variables
  clear ZSIG TEMP SALT

  % pressure
  rho_0 = 1020 ; % reference density (kg/m^3)
  g = 9.81 ;     % gravity (m/s^2)
  P = rho_0 .* g.* Z ;

  % density
  density = dens_wright_eos(T, S, P) ;

  % horizontal block average to one degree grid from 1/10 degree grid
  tbar = block_avg_scale(T, scale) ; % temperature
  sbar = block_avg_scale(S, scale) ; % salinity
  pbar = block_avg_scale(P, scale) ; % pressure

  % calculate delta rho
  rho_bar = block_avg_scale(density, scale) ; % true average density
  rho_mdl = dens_wright_eos(tbar, sbar, pbar) ; % what the model calcuates
  this_delta_rho = rho_bar - rho_mdl ; % error in density calcualtion

  % save
  delta_rho(:,:,:,jj) = this_delta_rho ;

  % second order tempature derivative
  [rho_TT, ~, ~] = dens_wright_sot(tbar, sbar, pbar) ;

  % subgrid scale temperature variance
  sig_TT = block_single_sd_scale(T, scale) ;

  % temperature term
  term_TT = 0.5 .* sig_TT .* rho_TT ;

  % save R^2
  r2 = r_squared(term_TT(:), this_delta_rho(:)) ;
  r2_temp_v_delta_rho(jj) = r2;

end % end isfile
end % end days
end % end months
end % end years

save('isopyc_one_yr_avg.mat', 'delta_rho', 'r2_temp_v_delta_rho')
