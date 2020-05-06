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

how_deep = 46 ;
scale = 10;
c = 0.20;

fprintf('Getting started on calculating the projection field for use in Stochastic Brankart. Scale = %g.\n', scale)

% load grid

% vertical distances
load('/glade/work/zofias/Brankart/Data/dz_zt_mask.mat', 'dz') % dz depth of layer (m)
dz = double(dz);
dz = dz(1:how_deep);
dz2 = reshape(dz, 1, 1, how_deep);
DZ = repmat(dz2, [2400/scale, 3600/scale, 1]) ;

% store projection fields over one year
proj_fields = zeros(2400/scale, 3600/scale, 73);
var_exp = zeros(2400/scale, 3600/scale, 73);

for ii = 1:73 ;
    fprintf('Iteration %d out of 73.\n', ii)
    load(filenames(ii), 'T')
    
    % only use top 2km
    T = T(:,:,1:how_deep);
    
    % horizontal block average to one degree grid from 1/10 degree grid
    tbar = block_avg_scale(T, scale) ;
    sig_TT = block_single_sd_scale(T, scale) ;
    
    % subgrid scale temperature variance
    sig_TT = block_single_sd_scale(T, scale) ;
    
    % mean model
    lgt = len_grad(tbar) ;
    
    % calculate projection
    inner_prod = sum(sig_TT .* lgt .* DZ, 3);
    norm_mean = c .* sum(lgt.^2 .* DZ, 3);
    proj_field = inner_prod./norm_mean;
    
    % save projection
    proj_fields(:,:,ii) = proj_field ;
    
    % relative error
    err = sig_TT - proj_field .* c.* lgt ;
    norm_err = sum(err.^2 .* DZ, 3);
    norm_true = sum(sig_TT.^2 .* DZ, 3);
    rel_err = sqrt(norm_err)./sqrt(norm_true);
    perc_var_exp = 1-rel_err.^2 ;
    
    % save fraction variance explained
    var_exp(:,:,ii) = perc_var_exp ;
    
end

save('stoch_proj_fields.mat', 'proj_fields', 'var_exp')
