load('../Data/stoch_proj_fields.mat', 'proj_fields')
load('../Data/ULAT_ULONG.mat', 'lat', 'lon')

log_fields = log(proj_fields);

mask = ~isnan(log_fields(:,:,1)) & abs(lat)<60 ;

% Estimate AR 1 coefficients
mdl = arima(1,0,0) ;
phis = zeros(240, 360) ;
for i = 1:240
    fprintf('Iteration %d out of 240.\n', i)
    for j = 1:360
        if mask(i, j)
            this = log_fields(i, j, :) ;
            est_mdl = estimate(mdl, this(:), 'Display', 'off') ;
            phis(i, j) = est_mdl.AR{1} ;
        end
    end
end

% save('AR_Phi_Final.mat', 'phis')
load('AR_Phi_Final.mat')

% Change to decorrelation time
phis(phis<=0) = NaN;
tau = -5 ./ log(phis) ;

% Estimate model for decorrelation time
load('../Data/surface_energy_0013.mat')
en_time = 100./yr_avg_surface_energy ; % seconds, 100 is correction for cm->m
en_time = en_time ./ (60*60*24) ; % days

new_mask = mask & ~isnan(tau+en_time) ;

b = mean( log(tau(new_mask)) - log(en_time(new_mask)) );


%% Make figures for paper
scale = 10;

load('../Data/ULAT_ULONG.mat')
ulat2 = circshift(lat, [0, 70]);
ulon2 = circshift(lon, [0, 70]);

tau_2 = circshift(log10(tau), [0, 70]);
tau_2(isnan(tau_2)) = 0;
en_time_2 = circshift(log10(exp(b+log(en_time))), [0, 70]);
en_time_2(isnan(en_time_2)) = 0;

% estimated decorrelation time
axesm('eckert4', 'MapLatLim', [-60, 60]); 
framem; gridm; tightmap
t = geoshow(ulat2, ulon2, tau_2, 'DisplayType', 'texturemap') ;
colormap default
set(gca, 'clim', [log10(4) log10(400)])
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5] );
t.AlphaDataMapping = 'none'; % interpet alpha values as transparency values
t.FaceAlpha = 'texturemap'; % Indicate that the transparency can be different each pixel
alpha(t,double(circshift(new_mask, [0,70]))); 

% colorbar
cvals = [5, 10, 30, 100, 300] ;
c = colorbar('Ticks', log10(cvals), 'TickLabels', cvals);
c.FontSize = 14;

%colormap pink
%old_cmap = colormap ;
%colormap ( flipud(old_cmap) ) ;

%grid
framem; tightmap; gridm;

%save
saveas(gcf, 'decorrelation_time_diagnosed.png')

% modeled decorrelation time
axesm('eckert4', 'MapLatLim', [-60, 60]); 
framem; gridm; tightmap
t = geoshow(ulat2, ulon2, en_time_2, 'DisplayType', 'texturemap') ;
colormap default
set(gca, 'clim', [log10(4) log10(400)])
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5] );
t.AlphaDataMapping = 'none'; % interpet alpha values as transparency values
t.FaceAlpha = 'texturemap'; % Indicate that the transparency can be different each pixel
alpha(t,double(circshift(new_mask, [0,70]))); 

% colorbar
cvals = [5, 10, 30, 100, 300] ;
c = colorbar('Ticks', log10(cvals), 'TickLabels', cvals);
c.FontSize = 14;

%grid
framem; tightmap; gridm;

%save
saveas(gcf, 'decorrelation_time_modeled.png')


% histogram 
[N,edges] = histcounts(log(proj_fields(45:205, :, :)), 'Normalization','pdf');
edges = edges(2:end) - (edges(2)-edges(1))/2;
de = diff(edges);
N2 = N ./ sum(N(2:end).*de) ;
plot(edges, N2, 'LineWidth', 2);
xlim([-5, 5]) ;
set(gca, 'FontSize', 14)
ylabel('Probability density');
xlabel('Log of projection')

%saveas(gcf, '../log_proj_pdf_diagnosed.png')

hold on
x = linspace(-5, 5, length(edges));
y = normpdf(x, 0, sqrt(0.39));
plot(x, y, '--', 'LineWidth', 2);
legend('Diagnosed', 'Modeled')
hold off

saveas(gcf, '../log_proj_pdf_diagnosed_modeled.png')

% QQ plot
this = log(proj_fields(45:205, :, :));
qqplot(this(:))
set(gca, 'FontSize', 14)
%title('Q-Q Plot')
title('')
ylabel('Quantiles of log of projection')
xlabel('Standard normal quantiles')
saveas(gcf, '../log_proj_qqplot.png')


% relative error
load('0013_01_05_1degree.mat', 'tbar', 'sgs_t_var', 'dz', 'ulat')
how_deep = 46;
c = 0.2;
scale=10;
tbar = tbar(:,:,1:how_deep) ; % average temperature, top 2km
sgs_t_var = sgs_t_var(:,:,1:how_deep); % subgrid scale temperature variance, top 2km
sgs_t_var_mdl = len_grad(tbar) ; % when multiplied by c this is our model for sgs temp variance 

dz = dz(1:how_deep);
dz2 = reshape(dz, 1, 1, how_deep);
DZ = repmat(dz2, [2400/scale, 3600/scale, 1]) ; % depth matrix

err = sgs_t_var - proj_fields(:,:,1) .* c.* sgs_t_var_mdl ;
norm_err = sum(err.^2 .* DZ, 3)./ sum(DZ,3);
norm_true = sum(sgs_t_var.^2 .* DZ, 3) ./ sum(DZ, 3);
rel_err = sqrt(norm_err)./sqrt(norm_true);

norm_mask = norm_true > 1e-4  & abs(ulat) < 60 & ~isnan(norm_err);
perc_var_exp = 100 .* (1-rel_err(norm_mask).^2);

% plot percent variance explained
[N3,edges] = histcounts(perc_var_exp, 'Normalization','pdf');
edges = edges(2:end) - (edges(2)-edges(1))/2;
de3 = diff(edges);
N4 = N3 ./ sum(N3(2:end).*de3) ;
plot(edges, N4, 'LineWidth', 2);
xlim([0, 100]) ;
set(gca, 'FontSize', 14)
ylabel('Probability density');
xlabel('Percent variance explained')
saveas(gcf, '../perc_var_exp_pdf.png')

%% East West Correlations

EarthRad = 6.371e6 ;

log_fields = log(proj_fields);
mask = ~isnan(log_fields(:,:,1)) & abs(lat)<60 ;


ln_mu_field = mean(log_fields, 3) ;
ln_sd_field = std(log_fields, 0, 3) ;

sd = mean(ln_sd_field(mask)) ;
mu = mean(ln_mu_field(mask)) ;

tau_mdl = exp(1.3) .* en_time ;
phi_mdl = exp(-5./tau_mdl);
phi_mdl(phi_mdl==0)=NaN;

innov = log_fields(:,:,2:73) - phi_mdl .* log_fields(:,:,1:72) ;
theta = innov ./ (sd.*sqrt(1-phi_mdl.^2)) ;
center = theta ;

nits = 20 ;
wests = zeros(240, 360, nits) ;
dists = zeros(240, 360, nits) ;
for i=1:nits
    w = circshift(theta, [0, i] ) ;
    d = ( (EarthRad .* pi) ./ 180) .* distance(lat, lon, circshift(lat, [0, i] ), circshift(lon, [0, i] )) ;
    w = mean( (w - mu) .* (center - mu), 3) ;
    w(~mask) = NaN ;
    d(~mask) = NaN ;
    wests(:,:, i) = w ;
    dists(:,:,i) = d ;
end

% plot east-west correlations

histogram2(dists/1e5, wests, 'DisplayStyle', 'tile') ;
ylim([-1 1]) ;
xlabel('Distance (100 km)') ;
ylabel('Sample correlation') ;
title('East-West Correlation') ;
xlim([0, 20])
set(gca,'fontsize', 14)
saveas(gcf, '../east_west_correlation.png')

%% North-South correlation

nits = 20 ;
norths = zeros(240, 360, nits) ;
dists_ns = zeros(240, 360, nits) ;
for i=1:nits
    n = circshift(theta, [i, 0] ) ;
    d = ( (EarthRad .* pi) ./ 180) .* distance(lat, lon, circshift(lat, [i, 0] ), circshift(lon, [i, 0] )) ;
    n = mean( (n - mu) .* (center - mu), 3) ;
    n(~mask) = NaN ;
    d(~mask) = NaN ;
    norths(:,:, i) = n ;
    dists_ns(:,:,i) = d ;
    dists_ns(:,:,i) = d ;
end


histogram2(dists_ns/1e5, norths, 'DisplayStyle', 'tile') ;
ylim([-1 1]) ;
xlabel('Distance (100 km)') ;
ylabel('Sample correlation') ;
title('North-South Correlation') ;
set(gca,'fontsize', 14)
xlim([0,20])
saveas(gcf, '../north_south_correlation.png')
