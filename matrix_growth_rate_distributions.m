% MATRIX_GROWTH_RATE_DISTRIBUTIONS builds and saves the hourly
% volume and Qc distributions of Pro and Syn SeaFlow data in preparation
% for the size-based matrix model.  Also grabs PAR from the sfl table in
% the SQL database.
%
% References:
% 
%   Sosik, et al, 2003.  Limnol. Oceanogr. 48:1756-1765.
%   Hunter-Cevera, et al, 2014.  PNAS. 111:9852-9857.
%   GitHub:  https://github.com/khuntercevera/phyto-division-rate-model
%
% Other files used:
%   *.db
%   suplabel.m
%   viridis.m
%
% Started:  08/Feb/2019 Annette Hynes, UW
% Modified: 11/Mar/2019 Fillgaps for hourly PAR
%           25/Mar/2019 Use shell to connect to db instead of 'sqlite'

close all; clearvars

%ds = datastore('smb://d-128-208-239-126.dhcp4.washington.edu/SeaFlow-OPP/latest');

%db_dir = '~/Documents/SeaFlow/Cruise_databases/';  % Directory containing SeaFlow db files
%db_list = dir([big_dir, '*.db']);                  % Structure containing all the .db files    
cruise_dir = '~/Documents/SeaFlow/Cruise_data/';    % Local directory containing cruise data

OPP_dir = '/Volumes/SeaFlow-OPP/latest/';           % Francois's Mac mini (smb://Francois Mac._smb._tcp.local/SeaFlow-OPP/latest)
this_cruise = 'SCOPE_16';                                 

%% Get required data from SQLite database
% Establishing a connection through Matlab's 'sqlite' command no longer
% works for me.  Send the SQL command to a shell via 'system' instead.

%dbfile = fullfile(db_dir, [this_cruise, '.db']);
%dbfile = fullfile([cruise_dir, this_cruise, '/'], [this_cruise, '.db']);
dbfile = fullfile([OPP_dir, this_cruise, '/'], [this_cruise, '.db']);

%conn = sqlite(dbfile);         % This function doesn't function 

%sqlquery = 'SELECT * FROM sqlite_master WHERE type="table"';        % list of tables
%results = fetch(conn, sqlquery);

%sqlquery = 'SELECT * FROM sfl';                                     % Get all data in specified table
%data = fetch(conn, sqlquery);

%sqlquery = 'SELECT sql FROM sqlite_master WHERE tbl_name = "sfl" AND type = "table"';   % Get fields in specified table
%extracted_data = fetch(conn, sqlquery);

%sqlquery = 'SELECT ifnull(conductivity,0) FROM sfl';
%data2 = fetch(conn, sqlquery);

%sqlquery = 'SELECT file, date, par FROM sfl';         % Take only the variables I need
%par_data = fetch(conn, sqlquery);

%sqlquery = 'SELECT file, date FROM sfl';         % Take only the variables I need
%par_data1 = fetch(conn, sqlquery);

%sqlquery = 'SELECT ifnull(par,-99) FROM sfl';         % Replace NULLs with a flag
%par_data2 = fetch(conn, sqlquery);
%ind = find([par_data2{:}] == -99);
%par_data2(cell2mat(cellfun(@(elem) elem == -99, par_data2, 'UniformOutput', false))) = {NaN};

%par_data = [par_data1, par_data2];
%T = cell2table(par_data, 'VariableNames', {'file', 'date', 'par'});

%sqlquery = 'SELECT file, date, par FROM sfl';          % Identifying specific variables somehow loses the variable names when the table is exported.
sqlquery = 'SELECT * FROM sfl';
outfile = [cruise_dir, this_cruise, '/', this_cruise, '_PAR.csv'];
cmd = ['sqlite3 -header -csv ', dbfile, ' "', sqlquery, ';" > ', outfile];
system(cmd);    

T = readtable(outfile);

sqlquery = 'SELECT file, flag FROM outlier';                % QA/QC flags
%flag_data = fetch(conn, sqlquery);
%T_flag = cell2table(flag_data, 'VariableNames', {'file', 'flag'});
outfile = [cruise_dir, this_cruise, '/', this_cruise, '_flag.csv'];
cmd = ['sqlite3 -header -csv ', dbfile, ' "', sqlquery, ';" > ', outfile];
system(cmd);    

T_flag = readtable(outfile);

%close(conn)

% Smooth and plot PAR 

time_dt = datetime(T.date, 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ssXXX','TimeZone','UTC');
T.time = time_dt;
time = datenum(time_dt);
T.par = double(T.par);
region = regionprops(isnan(T.par), 'Area', 'PixelIdxList');         % Find chunks of NaNs
chunk_length = [region.Area];
n_pts_24h = 24*60/3;                % Number of SeaFlow data points in 24 h
ind_chunk = {region.PixelIdxList};
ind_reNaN = {ind_chunk{chunk_length >= n_pts_24h}};

par_fill = fillgaps(T.par, n_pts_24h);                              % Fill in NaNs using autoregressive model for a day's worth of points
%par_fill = fillmissing(T.par, 'movmean', 5);                       % Fill in NaNs using moving mean
par_fill(par_fill < 0) = 0;
par_fill(vertcat(ind_reNaN{:})) = NaN;                                 % Full days of missing data shouldn't be filled

par_smooth = smoothdata(par_fill, 'gaussian', 20); % Smooth the light data
T.par_smooth = par_smooth;

partime = table2timetable(T);
par_hourly = retime(timetable(partime.time, partime.par_smooth, 'VariableNames', {'PAR'}), 'hourly', 'mean');         % Bin by hour

parseries = timeseries(T.par_smooth, time, 'Name', 'PAR_smooth');
ten_min = 10/(60*24);   % Get PAR in 10-min time steps to fit matrix model
t_10min = datenum(par_hourly.Time(1)):ten_min:datenum(par_hourly.Time(end));
par_10min = resample(parseries, t_10min);

colores = viridis(15);

figure(1)
set(1, 'Position', [100 100 1500 750], 'PaperPositionMode', 'auto')

    l11 = line(time, T.par, 'Marker', '.', 'LineStyle', 'none', 'Color', colores(2, :));
    l12 = line(time, par_smooth, 'Color', colores(5, :), 'LineWidth', 2);
    l13 = line(par_10min.Time, par_10min.Data, 'Color', colores(10, :), 'LineWidth', 2);
    l14 = line(datenum(par_hourly.Time), par_hourly.PAR, 'Color', colores(14, :), 'Marker', 'x');
    set(gca, 'FontSize', 15)
    datetick('x')
    xlabel('Date')
    ylabel('PAR (\mumol quanta m^{-2} s^{-1})')
    leg1 = legend('Data', 'Smoothed', '10-min', 'Hourly');
    title(['Photosynthetically Active Radiation, ', strrep(this_cruise, '_', ' ')])
    
print(1, '-dpng', [cruise_dir, this_cruise, '/', this_cruise, '_PAR.png'])

% Remove flagged files from the list

clean_files = char(table2array(T_flag(T_flag.flag == 0, 1)));       % unflagged files only

[~, ia, ~] = intersect(T.file, clean_files);
T_clean = T(ia, :);

writetable(timetable2table(par_hourly), [cruise_dir, this_cruise, '/', this_cruise, '_PAR_hourly_sfl.csv'])

T_10min = table(datetime(par_10min.Time, 'ConvertFrom', 'datenum'), par_10min.Data, 'VariableNames', {'Time', 'PAR_10min'});
writetable(T_10min, [cruise_dir, this_cruise, '/', this_cruise, '_PAR_10min_sfl.csv'])

%% Get size distributions from VCT files

% From SeaFlow databases (22/Feb/2018) looking at Prochlorococcus and
% Synechococcus:
%--- Range of mean diameter:    pro = [0.4987 1.3862] (um)
%                               syn = [0.7190 1.6655]
%--- Range of mean volume:      pro = [0.5194 11.1579] (um^3)
%                               syn = [1.5567 19.3525]
%--- Range of mean Qc:          pro = [0.01430 0.3141] (pg C cell^-1)
%                               syn = [0.04281 1.0310]

m = 57;
size_ind = 1:m;                 % indices of size bins
del_v = 0.125;                  % Logarithmic "space" between size classes; inverse must be an integer

% vmin = 2^(-5);                  % Hunter-Cevera (2014), note typo on page 9853.
vmin = 2^(-3);
log2_bins = log2(vmin) + (size_ind - 1)*del_v;      % Hunter-Cevera (2014) values
vol_bins = 2.^(log2_bins);              % um^3

Qmin = 2^(-8);
log2_bins = log2(Qmin) + (size_ind - 1)*del_v;  
Qc_bins_pro = 2.^(log2_bins);               % pg C cell^-1      

Qmin = 2^(-7);
log2_bins = log2(Qmin) + (size_ind - 1)*del_v;  
Qc_bins_syn = 2.^(log2_bins);               % pg C cell^-1      

vct_dir = [OPP_dir, this_cruise, '/', this_cruise, '_vct/50']; % Specify quantile

t_hour = dateshift(datetime(time(1):(1/24):time(end), 'ConvertFrom', 'datenum', 'TimeZone', 'UTC'), 'start', 'hour');      % Hourly time points
nh = length(t_hour);

vol_dist_pro = zeros(m, nh - 1);
vol_dist_syn = zeros(m, nh - 1);
Qc_dist_pro = zeros(m, nh - 1);
Qc_dist_syn = zeros(m, nh - 1);

for ii = 1:(nh - 1)           % Hourly 
    disp(t_hour(ii))
    file_list = T.file(T.time >= t_hour(ii) & T.time < t_hour(ii + 1));
    
    vol_pro_vec = [];
    vol_syn_vec = [];
    Qc_pro_vec = [];
    Qc_syn_vec = [];
    
    Nvp = NaN*ones(1, m);
    Nvs = NaN*ones(1, m);
    NQp = NaN*ones(1, m);
    NQs = NaN*ones(1, m);
       
    
    for jj = 1:length(file_list)
        this_file = file_list(jj);
        this_file = [vct_dir, '/', char(this_file), '.vct.gz'];     % Seaflow file identifying each particle
        if isfile(char(this_file))    % Confirm that file exists.  Otherwise go to next step
               
            vct_file = gunzip(this_file);
            fid = fopen(vct_file{1});
            vct = textscan(fid, '%f %f %f %f %f %f %s');        % "diam_lwr", "Qc_lwr", "diam_mid", "Qc_mid", "diam_upr", "Qc_upr", "pop"   
            fclose(fid);
            
            diam = vct{3};          % Diameter based on middle value of index of refraction (um)
            vol = (4*pi/3)*diam.^3; % Assume sphere (um^3)
            Qc = vct{4};            % Carbon quota based on middle value of index of refraction (pg C cell^-1)
            pop = vct{7};           % beads, croco, picoeuk, synecho, prochloro
    
            vol_pro = vol(strcmp(pop, 'prochloro'));  % Prochlorococcus volume
            vol_pro_vec = [vol_pro_vec; vol_pro];
            vol_syn = vol(strcmp(pop, 'synecho'));    % Synechococcus
            vol_syn_vec = [vol_syn_vec; vol_syn];
            Qc_pro = Qc(strcmp(pop, 'prochloro'));  % Prochlorococcus carbon quota
            Qc_pro_vec = [Qc_pro_vec; Qc_pro];
            Qc_syn = Qc(strcmp(pop, 'synecho'));    % Synechococcus
            Qc_syn_vec = [Qc_syn_vec; Qc_syn];
        end
    end
    
    [Nvp, ~] = histcounts(vol_pro_vec, [0, vol_bins]);
    vol_dist_pro(:, ii) = Nvp;
        
    [Nvs, ~] = histcounts(vol_syn_vec, [0, vol_bins]);
    vol_dist_syn(:, ii) = Nvs;
        
    [NQp, ~] = histcounts(Qc_pro_vec, [0, Qc_bins_pro]);
    Qc_dist_pro(:, ii) = NQp;
        
    [NQs, ~] = histcounts(Qc_syn_vec, [0, Qc_bins_syn]);
    Qc_dist_syn(:, ii) = NQs;
end

save([cruise_dir, this_cruise, '/', this_cruise, '_size_dist.mat'], 'vol_dist_pro', 'vol_dist_syn', 'Qc_dist_pro', 'Qc_dist_syn', 't_hour', 'vol_bins', 'Qc_bins_pro', 'Qc_bins_syn');

%% Graph the distributions in time
figure(2)
set(2, 'Position', [150 150 1500 1500], 'PaperPositionMode', 'auto', 'Colormap', colormap(viridis()))

ax2 = zeros(4, 1);

ax2(1) = subplot(221);
    pcolor(datenum(t_hour(1:end-1)), vol_bins, vol_dist_pro); shading interp
    datetick('x', 'mmm dd', 'keeplimits')
    ylabel('Volume (\mum^3)')
   
ax2(2) = subplot(222);
    waterfall(vol_bins, datenum(t_hour(1:end-1)), vol_dist_pro')
    datetick('y', 'mmm dd', 'keeplimits')
    xlabel('Volume (\mum^3)')
    
ax2(3) = subplot(223);
    pcolor(datenum(t_hour(1:end-1)), Qc_bins_pro, Qc_dist_pro); shading interp    
    datetick('x', 'mmm dd', 'keeplimits')
    ylabel('Carbon Quota (pg C cell^{-1})')
    
ax2(4) = subplot(224);
    waterfall(Qc_bins_pro, datenum(t_hour(1:end-1)), Qc_dist_pro')
    datetick('y', 'mmm dd', 'keeplimits')
    xlabel('Carbon Quota (pg C cell^{-1})')

set(ax2, 'FontSize', 13)    
[axt2, ht2] = suplabel('\it Prochlorococcus', 't'); 
set(axt2, 'FontSize', 15)

print(2, '-dpng', [cruise_dir, this_cruise, '/', this_cruise, '_Pro_dist.png'])

figure(3)    
set(3, 'Position', [200 200 1500 1500], 'PaperPositionMode', 'auto', 'Colormap', colormap(viridis()))

ax3 = zeros(4, 1);

ax3(1) = subplot(221);
    pcolor(datenum(datenum(t_hour(1:end-1))), vol_bins, vol_dist_syn); shading interp
    datetick('x', 'mmm dd', 'keeplimits')
    ylabel('Volume (\mum^3)')
   
ax3(2) = subplot(222);
    waterfall(vol_bins, datenum(t_hour(1:end-1)), vol_dist_syn')
    datetick('y', 'mmm dd', 'keeplimits')
    xlabel('Volume (\mum^3)')
    
ax3(3) = subplot(223);
    pcolor(datenum(t_hour(1:end-1)), Qc_bins_syn, Qc_dist_syn); shading interp    
    datetick('x', 'mmm dd', 'keeplimits')
    ylabel('Carbon Quota (pg C cell^{-1})')
    
ax3(4) = subplot(224);
    waterfall(Qc_bins_syn, datenum(t_hour(1:end-1)), Qc_dist_syn')
    datetick('y', 'mmm dd', 'keeplimits')
    xlabel('Carbon Quota (pg C cell^{-1})')
    
set(ax3, 'FontSize', 13)        
[axt3, ht3] = suplabel('\it Synechococcus', 't'); 
set(axt3, 'FontSize', 15)

print(3, '-dpng', [cruise_dir, this_cruise, '/', this_cruise, '_Syn_dist.png'])