%                        CMPU250 - Professor Eric Aaron
%                       Drug Simulation - Kyle Patterson
% March 2018

% Constants
dosage_in_mg = 100; % amount of drug in single dose of a pill
dosage_in_ug = dosage_in_mg * 1000;
absorbtion_frac = 0.12; % fraction of drug absorbed by pill
adult_plasma_in_l = 4.7;
adult_plasma_in_ml = adult_plasma_in_l * 1000;
halflife_drug = 22; % in hours
rate_constant = log(2) / halflife_drug; % Exponential decay: 1/2 = e^{-kt}

% Parameters
hours = 75;
interval = 24; % number of hours between when each pill is taken
num_points = 25; % how many points to graph, not including t = 0

% Intermediate calculations
sub_divisions = floor((num_points) / (hours/interval));
% Divisions is approximately equal to `num_points`.
divisions = (hours / interval) * sub_divisions;
delta_time = hours / divisions;

% initial concentration of drug added by each pill
init_conce = dosage_in_ug * absorbtion_frac / adult_plasma_in_ml;

% --- Summation method ---
% Graph generation
times = 0:delta_time:hours; % fill 1D array with numbers of hours passed
conces = zeros([1, length(times)]); % initialize empty array of dim(times)
k = rate_constant;
I = interval;
for i = 1:length(times)
    t = times(i);
    N = floor(t / interval) + 1; % number of pills
    sum_of_decays = exp(-k*t) * (exp(k*I*N) - 1) / (exp(k*I) - 1);
    conces(i) = sum_of_decays * init_conce;
end

% Plot graph
drug_conce = figure();
plot(times, conces, 'b', 'DisplayName', 'Summation');
title('Drug concentration versus time')
xlabel('Time (hours)')
ylabel('Concentration (ug/mL)')
    
% Point calculation (concentration after given number of hours)
N = floor(hours / interval) + 1; % number of pills taken
sum_of_decays = exp(-k*t) * (exp(k*I*N) - 1) / (exp(k*I) - 1);
concentration = sum_of_decays * init_conce;

% There is no equivalent point calculation for the tabular method,
% since one must compute the intermediate values to reach a solution

% --- Tabular method ---
% Graph generation
times2 = 0:delta_time:hours; % fill 1D array with numbers of hours passed
conces2 = zeros([1, length(times2)]); % initialize empty array of dim(times)
conces2(1) = init_conce; % setup initial value
for i = 2:length(times2)
    t = times2(i);
    conces2(i) = (mod(t, interval) == 0) * init_conce + ...
        conces2(i - 1) * exp(-rate_constant * delta_time);
end

% Plot graph
hold;
plot(times2, conces2, '--r', 'DisplayName', 'Tabular');
xlabel('Time (hours)')
ylabel('Concentration (ug/mL)')
legend('show')

