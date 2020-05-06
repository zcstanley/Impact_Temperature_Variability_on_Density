load('/glade/work/zofias/Brankart/Data/dxt_dyt.mat')
% dxt: x-distance between midpoints of T cells
% dyt: y-distance between midpoints of T cells

load('/glade/work/zofias/Brankart/Data/dz_zt_mask.mat', 'z_t', 'dz') 
% z_t depth from surface to midpoint of layer (m)
% dz depth of layer (m)

how_deep = 62 ;

% coefficient c caluclated separately with
% get_coef_one_yr_avg

% half degree
scale = 5;
c = 0.14;

% one degree
%scale = 10 ;
%c = 0.2 ; 

% two degrees
%scale = 20;
%c=0.24;

rho_0 = 1020 ; % reference density (kg/m^3)
g = 9.81 ;     % gravity (m/s^2)
pressure = rho_0 .* g.* z_t(1:how_deep) ; % pressure at midpoint of layer (Pa)
p2 = reshape(pressure, 1, 1, how_deep) ;
P = repmat(p2, [2400, 3600, 1]) ;
pbar = repmat(p2, [2400/scale, 3600/scale, 1]) ;

dxt(dxt<0) = NaN ;
dyt(dyt<0) = NaN ;
cdxt = scale * block_avg_scale(dxt, scale) ;
cdyt = scale * block_avg_scale(dyt, scale) ;

date_strs = ["01-05", "01-10", "01-15", "01-20", "01-25", "01-30", ...
                "02-04", "02-09", "02-14", "02-19", "02-24", ...
                "03-01", "03-06", "03-11", "03-16", "03-21", "03-26", "03-31", ...
                "04-05", "04-10", "04-15", "04-20", "04-25", "04-30", ...
                "05-05", "05-10", "05-15", "05-20", "05-25", "05-30", ...
                "06-04", "06-09", "06-14", "06-19", "06-24", "06-29", ...
                "07-04", "07-09", "07-14", "07-19", "07-24", "07-29", ...
                "08-03", "08-08", "08-13", "08-18", "08-23", "08-28", ...
                "09-02", "09-07", "09-12", "09-17", "09-22", "09-27", ...
                "10-02", "10-07", "10-12", "10-17", "10-22", "10-27", ...
                "11-01", "11-06", "11-11", "11-16", "11-21", "11-26", ...
                "12-01", "12-06", "12-11", "12-16", "12-21", "12-26", "12-31" ] ;

filenames = strcat("/glade/work/zofias/Brankart/Data/T_S/proj_field_0013_", date_strs, ".mat") ;

% calculate over a year
delrho_pos = zeros(73, 1) ;
r2_delrho_rhostar = zeros(73, 1) ;
mean_TT = zeros(73, 1) ;
mean_TS = zeros(73, 1) ;
mean_SS = zeros(73, 1) ;
r2_delrho_termTT = zeros(73, 1) ;
%hub_slope = zeros(73, 1) ;
r2_sigTT_lgt = zeros(73, 1) ;
r2_delrho_mdltermTT = zeros(73, 1) ;
r2_dx_delrho_mdltermTT = zeros(73, 1) ;
r2_dy_delrho_mdltermTT = zeros(73, 1) ;

for ii = 1:73 ;
    
    fprintf('Iteration %d out of 73\n', ii)
    
    % load data
    load(filenames(ii), 'T', 'S')
    
    % calculate density
    T(S<0) = NaN ; % salinity should not be negative
    S(S<0) = NaN ; % salinity should not be negative
    density = dens_wright_eos(T, S, P) ;
    
    % calculate spatial (block) averages
    tbar = block_avg_scale(T, scale) ;
    sbar = block_avg_scale(S, scale) ;
    rho_bar = block_avg_scale(density, scale) ;
    
    % error in density as a result of averaging
    rho_mdl = dens_wright_eos(tbar, sbar, pbar) ;
    delta_rho = rho_bar - rho_mdl ;
    
    % second order derivatives
    [rho_TT rho_TS rho_SS] = dens_wright_sot(tbar, sbar, pbar) ;
    
    % (co)variances 
    [sig_TT sig_TS sig_SS] = block_sd_scale( T, S, scale) ;
    
    % three terms
    term_TT = 0.5 .* sig_TT .* rho_TT ;
    term_TS = sig_TS .* rho_TS ;
    term_SS = 0.5 .* sig_SS .* rho_SS ;
    
    % dx * gradient of T
    lgt = len_grad(tbar) ;
    
    % estimate slope from minimizing hubers loss
    nonan = ~isnan(sig_TT(:)+lgt(:)) ;
    x = lgt(nonan) ;
    y = sig_TT(nonan) ;
    %betas = 0.1:.01:0.25 ;
    %huber = zeros(size(betas)) ;
    %for j = 1:length(betas)
    %    huber(j) = huber_loss(x, y, betas(j)) ;
    %end
    %[C, I] = min(huber) ;
    %c = betas(I)  ;
    
    % our full model of the error in density
    mdl_term_TT = 0.5 .* c .* lgt .* rho_TT ;
    
    % gradients
    [alpha_x, beta_x] = dens_wright_alpha_beta(0.5*(tbar + circshift(tbar, [0 -1 0])), 0.5 .* (sbar + circshift(sbar, [0 -1 0])), pbar) ;
    [alpha_y, beta_y] = dens_wright_alpha_beta(0.5*(tbar + circshift(tbar, [-1 0 0])), 0.5 .* (sbar + circshift(sbar, [-1 0 0])), pbar) ;
    dx_tbar = gradx(tbar, cdxt) ;
    dx_sbar = gradx(sbar, cdxt) ;
    dy_tbar = grady(tbar, cdyt) ;
    dy_sbar = grady(sbar, cdyt) ;
    dx_mdl_term_TT = gradx(mdl_term_TT, cdxt) ;
    delta_dx_rho = gradx(rho_bar, cdxt)  - (alpha_x .* dx_tbar + beta_x .* dx_sbar)  ;
    dy_mdl_term_TT = grady(mdl_term_TT, cdyt) ;
    delta_dy_rho = grady(rho_bar, cdyt)  - (alpha_y .* dy_tbar + beta_y .* dy_sbar)  ;

    % save huber's loss slope
    %hub_slope(ii) = betas(I) ;
    
    % save means
    mean_TT(ii) = nanmean(abs(term_TT(:))) ;
    mean_TS(ii) = nanmean(abs(term_TS(:))) ;
    mean_SS(ii) = nanmean(abs(term_SS(:))) ;
    
    % save percent of time delta rho is positive
    delrho_pos(ii) = nanmean( delta_rho(:) > 0 ) ;
    
    % save r^2 values
    r2_delrho_rhostar(ii) = r_squared(term_TT+term_TS+term_SS, delta_rho) ;
    r2_delrho_termTT(ii) = r_squared( term_TT, delta_rho) ;
    r2_sigTT_lgt(ii) = r_squared(c.*lgt, sig_TT) ;
    r2_delrho_mdltermTT(ii) = r_squared(mdl_term_TT, delta_rho) ;
    r2_dx_delrho_mdltermTT(ii) = r_squared(dx_mdl_term_TT, delta_dx_rho) ;
    r2_dy_delrho_mdltermTT(ii) = r_squared(dy_mdl_term_TT, delta_dy_rho) ;
    
end

save('one_yr_avg_r2_scale_5.mat', 'delrho_pos', 'r2_delrho_rhostar',... 
    'mean_TT', 'mean_TS', 'mean_SS', ...
    'r2_delrho_termTT', 'r2_sigTT_lgt',... 
    'r2_delrho_mdltermTT', 'r2_dx_delrho_mdltermTT', 'r2_dy_delrho_mdltermTT' )
