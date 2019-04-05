% MATRIX_GROWTH_RATE_DRIVER applies the size-based matrix model to Pro and
% Syn from SeaFlow data to calculate hourly growth rates.
%
% References:
% 
%   Sosik, et al, 2003.  Limnol. Oceanogr. 48:1756-1765.
%   Hunter-Cevera, et al, 2014.  PNAS. 111:9852-9857.
%   GitHub:  https://github.com/khuntercevera/phyto-division-rate-model
%
% Other files used:
%   *_size_dist.mat
%   *_PAR.csv']);
%   suplabel.m
%   viridis.m
%
% Started:  04/Mar/2019 Annette Hynes, UW
% Modified:

%% Load data

close all; clearvars

data_dir = '~/Documents/SeaFlow/Cruise_data/';  % Directory containing cruise directories with VCT data
this_cruise = 'SCOPE_16';

load([data_dir, this_cruise, '/', this_cruise, '_size_dist.mat']); % Contains vol_dist_pro, vol_dist_syn, Qc_dist_pro, Qc_dist_syn, t_hour, vol_bins, Qc_bins_pro, Qc_bins_syn

%T_PAR_hourly = readtable([data_dir, this_cruise, '/', this_cruise, '_PAR_hourly.csv']);
T_PAR_hourly = readtable([data_dir, this_cruise, '/', this_cruise, '_PAR_hourly_sfl.csv']);
T_PAR_hourly.Time = datetime(T_PAR_hourly.Time, 'TimeZone', 'UTC');      % Need to carry time zone

%T_PAR_10min = readtable([data_dir, this_cruise, '/', this_cruise, '_PAR_10min.csv']);
T_PAR_10min = readtable([data_dir, this_cruise, '/', this_cruise, '_PAR_10min_sfl.csv']);
T_PAR_10min.Time = datetime(T_PAR_10min.Time, 'TimeZone', 'UTC');      % Need to carry time zone

%% Optimize the parameters for size-based matrix model for each day

% dt inverse must be an integer

% From Hunter-Cevera's call_to_opt_mvco.m
%----------------------------------------

%hr1 = 7; hr2 = 25; %time window for the model used in Sosik, 2003 and Hunter-Cevera, 2014
hr1 = 1; hr2 = 25;  % Division allowed all day, every day

restitles = {'day'; 'gmax1'; 'b1'; 'E*1'; 'dmax1'; 'gmax2'; 'b2'; 'E*2'; 'dmax2'; 'proportion'; 'm1'; 'm2'; 'sigma1'; 'sigma2';' s'; '-logL'; 'mu'; 'mu1'; 'mu2'; 'ending proportion 1'; 'ending proportion 2'; 'exitflag'; 'number solver runs'};
ms = MultiStart('Display', 'off', 'TolX', 1e-5, 'UseParallel', 'always', 'StartPointsToRun', 'bounds');
opts = optimset('Display', 'off', 'TolX', 1e-8, 'Algorithm', 'interior-point', 'UseParallel', 'always', 'MaxIter', 3000, 'MaxFunEvals', 10000);
icsTol = 0.2;
tolvec = [0.01 0.01 100 0.005 0.01 0.01 100 0.005 0.01 0.5 0.5 0.5 0.5 10];
lb=-[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 5 5 1 1 1e-4]; %parameter bounds

a1 = -1*eye(14); %set parameter bounds to be interpreted by fmincon
a2 = eye(14);
A = zeros(28,14);
A(1:2:27,:) = a1;
A(2:2:28,:) = a2;

% Step through each day and start at dawn
doy = day(t_hour(1:end-1), 'dayofyear');
day_list = unique(doy);
nd = length(day_list);

doy_par = day(T_PAR_hourly.Time, 'dayofyear');   % May be different from sfl if taken directly from ship PAR files

modelresults = NaN*ones(nd, 23);
allmodelruns = cell(nd, 2);

for ii = 1:nd
    this_day = day_list(ii);
    ind_day_sfl = find(doy == this_day);
    ind_day_par = find(doy_par == this_day);
    
    disp(this_day)
    % set dawn to 0
    
    this_PAR = T_PAR_hourly.PAR(ind_day_par);
        
    if length(this_PAR) > 10                % Want mostly full days, dawn is around 14:00 GMT
        if max(this_PAR) > 300              % Ain't no sunshine when she's gone
            del_par = diff(this_PAR);
            ind_rise = find(del_par > 1);   % Significant jump in PAR
            ind_dawn = ind_rise(1) + 1;     % Index offset due to diff
            ind_start = ind_day_par(ind_dawn);
            t_today = T_PAR_hourly.Time(ind_start:ind_start + 24);
            
            [~, ~, ind_dist] = intersect(t_today, t_hour);  % Match the hourly cell distribution times to the 24 hour of PAR times
            dist_today = Qc_dist_pro(:, ind_dist);
    
            ind_10 = find(T_PAR_10min.Time >= t_today(1) & T_PAR_10min.Time <= t_today(end));
            PAR_today = T_PAR_10min.PAR_10min(ind_10);
       
            ub=[1 15 max(PAR_today) 1 1 15 max(PAR_today) 1 0.5 50 50 15 15 1e4];
    
            B=zeros(27,1);
            B(1:2:27)=lb;
            B(2:2:28)=ub;

            %starting conditions:
            x0 = [0.2*rand 6*rand max(PAR_today)*rand 0.1*rand 0.2*rand 6*rand max(PAR_today)*rand 0.1*rand 0.5*rand 30*rand+20 30*rand+20 10*rand+2 10*rand+2 1e4*rand];

            %random start points:
            tpoints = CustomStartPointSet([0.2*rand(40,1) 6*rand(40,1) max(PAR_today)*rand(40,1) 0.1*rand(40,1) 0.2*rand(40,1) 6*rand(40,1) max(PAR_today)*rand(40,1) 0.1*rand(40,1) 0.5*rand(40,1) 30*rand(40,1)+20 30*rand(40,1)+20 10*rand(40,1)+2 10*rand(40,1)+2 1e4*rand(40,1)]);

            problem = createOptimProblem('fmincon', 'x0', x0, 'objective', @(theta) negloglike_calc(PAR_today, dist_today, theta, vol_bins, hr1, hr2), 'Aineq', A, 'bineq', B, 'options', opts);
            
            tic
            [xmin, fmin, exitflag, ~, soln] = run(ms, problem, tpoints);
            toc
            
            tic
            [mu, mu1, mu2, p1, p2] = growth_rate(PAR_today, vol_bins, dist_today, xmin, hr1, hr2);
            toc
            
            modelresults(ii, :)=[this_day xmin fmin mu mu1 mu2 p1 p2 exitflag NaN];
            %allmodelruns{ii, 1} = modelfits;
            %allmodelruns{ii, 2} = allstarts;

            eval(['save ' data_dir, this_cruise, '/', this_cruise,  ' modelresults allmodelruns'])
        else
            disp('No sun!')
        end
    else
        disp('Partial day.')
    end
end
%---------------------------------------

%% Run model with optimized parameters
savepath = [data_dir, this_cruise, '/'];
theta = [0.0001, 0.0020, 0.4000, 0.0001, 0.0001, 0.0040, 0.6000, 0.0000, 0.0003, 0.0270, 0.0330, 0.0050, 1.5000];
theta_labels = {'gmax1'; 
    'b1';
    'E*1';
    'dmax1';
    'gmax2';
    'b2';
    'E*2';
    'dmax2';
    'proportion';
    'mean 1';
    'mean 2';
    'sigma';
    
    's'};