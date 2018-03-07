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
% interval MUST be an integer of function will break
interval = 24; % number of hours between when each pill is taken
num_points = 50; % how many points to graph, including t = 0

% Intermediate calculations
delta_time = hours/(num_points - 1); % size of each time step

% initial concentration of drug added by each pill
init_conce = dosage_in_ug * absorbtion_frac / adult_plasma_in_ml;

% --- Summation method ---
% Graph generation
times = 0:delta_time:hours; % fill 1D array with numbers of hours passed
conces = zeros([1, num_points]); % initialize empty array of dim(times)
syms x
for i = 1:num_points
    t = times(i);
    num_pills = floor(t / interval) + 1;
    sum_of_decays = symsum(exp(-rate_constant * (t - x * interval)) ...
    , x, 0, num_pills - 1);
    conces(i) = sum_of_decays * init_conce;
end

% Plot graph
drug_conce = figure();
plot(times, conces, 'b', 'DisplayName', 'Summation');
title('[Summation] Drug concentration versus time')
xlabel('Time (hours)')
ylabel('Concentration (ug/mL)')
    
% Point calculation
% calculate concentration for given number of hours
% total number of pills taken in time in time interval
num_pills = floor(hours / interval) + 1;
sum_of_decays = symsum(exp(-rate_constant * (hours - x * interval)) ...
    , x, 0, num_pills - 1);
concentration = sum_of_decays * init_conce;

% There is no equivalent point calculation for the tabular method,
% since one must compute the intermediate values to reach a solution

% --- Tabular method ---
% Graph generation
sub_divisions = floor((num_points) / (hours/interval));
% Evaluated region must be subdivided by such that we hit multiples
% of `interval` with `delta_time` in order to add pills.
% Divisions is approximately equal to `num_points`.
divisions = (hours / interval) * sub_divisions;
delta_time = hours / divisions;
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

