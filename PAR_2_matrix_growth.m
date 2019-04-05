% PAR_2_MATRIX_GROWTH uses PAR data recorded directly from the ship to get 10-min
% bins in preparation for using the size-structured matrix model (Sosik et
% al, 2003; Hunter-Cevera et al, 2014).  This is used in case SeaFlow did
% not properly record PAR underway.
%
% Other files used:
%   viridis.m
%   MGL-par.txt
%
% Started:  13/Mar/2019 Annette Hynes, UW
% Modified:

close all; clearvars

cruise_dir = '~/Documents/SeaFlow/Cruise_data/';  % Directory containing cruise directories with VCT and other data
this_cruise = 'MGL1704';        

%% Import PAR data

par_file = [cruise_dir, this_cruise, '/MGL-par.txt'];
opts = detectImportOptions(par_file);
opts.Delimiter = {',', '\t'};
opts.VariableTypes = {'char', 'char', 'double', 'double', 'double'};
opts.VariableNames = {'label', 'DateTime', 'PAR', 'thing1', 'thing2'};
T = readtable(par_file, opts);
time = datetime(pad(T.DateTime, 'right', '0'), 'InputFormat', 'u:DDD:HH:mm:ss.SSSS', 'TimeZone', 'UTC');

T.time = time;

%% Get hourly and 10 min means
partime = table2timetable(T);

par_hourly = retime(timetable(partime.time, partime.PAR, 'VariableNames', {'PAR'}), 'hourly', 'mean');         % Bin by hour
par_hourly.PAR(par_hourly.PAR < 0) = 0;       % No negative light

dt = minutes(10);            % 10 minutes
par_10min = retime(timetable(partime.time, partime.PAR, 'VariableNames', {'PAR'}), 'regular', 'mean', 'TimeStep', dt);         % Bin by 10 min
par_10min.PAR(par_10min.PAR < 0) = 0;       % No negative light

%% Plot it

colores = viridis(15);

figure(1)
set(1, 'Position', [100 100 1500 750], 'PaperPositionMode', 'auto')

    l11 = line(time, T.PAR, 'Marker', '.', 'LineStyle', 'none', 'Color', colores(2, :));
    l12 = line(par_10min.Time, par_10min.PAR, 'Color', colores(10, :), 'LineWidth', 2);
    l13 = line(par_hourly.Time, par_hourly.PAR, 'Color', colores(14, :), 'Marker', 'x');
       
    set(gca, 'FontSize', 15)
    datetick('x')
    xlabel('Date')
    ylabel('PAR (\mumol quanta m^{-2} s^{-1})')
    leg1 = legend('1-s Data', '10-min', 'hourly');
    title(['Photosynthetically Active Radiation, ', strrep(this_cruise, '_', ' ')])
   
print(1, '-dpng', [cruise_dir, this_cruise, '/', this_cruise, '_PAR_ship.png']) 

%% Export PAR

T_hourly = table(par_hourly.Time, par_hourly.PAR, 'VariableNames', {'Time', 'PAR'});
writetable(T_hourly, [cruise_dir, this_cruise, '/', this_cruise, '_PAR_hourly.csv'])

T_10min = table(par_10min.Time, par_10min.PAR, 'VariableNames', {'Time', 'PAR_10min'});
writetable(T_10min, [cruise_dir, this_cruise, '/', this_cruise, '_PAR_10min.csv'])