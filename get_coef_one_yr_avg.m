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

% resolution
% sample scales:  5 - half degree, 10 - one degree, 20 - two degrees
% scale = 2, 3, 4, 5, 8, 10, 12, 15, 16, 20 

scale = 2;
filename = sprintf('one_yr_avg_coef_scale_%g.mat', scale) ;

% calculate over a year
hub_slope = zeros(73, 1) ;

for ii = 1:73 ;
    
fprintf('Iteration %d out of 73. Scale %g.\n', ii, scale)
    
    % load data
    load(filenames(ii), 'T')
    
    % save hubers loss slope
    hub_slope(ii) = coef(T, scale) ;
    
end

save(filename, 'hub_slope')
