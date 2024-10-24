clear; clc; close all;

% Parameters
r_mean = 1.0; % Mean growth rate
S = 2:2:60; % Number of species (only 40 species in this case)
DD = 1e-5; % Migration rate
A_mean_range = [0.0 1]; % Range for mean interaction strength (inhibition alpha)
delta = 0.1; % pH dynamics parameter
beta_range = [0.00 0.1] * delta; % Beta (interaction strength affected by pH)
ystep = 30; % Number of steps for A_mean
T_real = 2000; % Total evolution time
step = 0.1; % Time interval for integration
T = T_real / step; % Number of steps in the simulation
num_sim = 1000; % Number of simulations per parameter set
period = 100; % Time period for calculating species richness
abundance_threshold = 8e-4; % Threshold for survival of species

% Preallocate matrices for storing results
lyapExp_record = zeros(length(S) * ystep, num_sim); % Lyapunov exponents (for chaos detection)
fluc_record = zeros(length(S) * ystep, num_sim); % Record fluctuations before invasion
fluc_record_final = zeros(length(S) * ystep, num_sim); % Record fluctuations after invasion
diversity_record = zeros(length(S) * ystep, num_sim); % Diversity before invasion
invade_record = zeros(length(S) * ystep, num_sim); % Record whether invasion succeeds
invasion_effect = zeros(length(S) * ystep, num_sim); % Long-term invasion effect on community structure

% Generate arrays for mean interaction strength (A_mean) and beta
A_mean = linspace(A_mean_range(1), A_mean_range(2), ystep);
beta = linspace(beta_range(1), beta_range(2), ystep);

% Main simulation loop
for mm = 1:length(S)
    for pp = 1:ystep
        for hh = 1:num_sim
            
            % Compute species composition and invasion success using pH-influenced dynamics
            composition = pH_compute_invasion(S(mm), DD, A_mean(pp), beta(pp), delta, T, step);
            
            % Calculate pre-invasion diversity
            diversity_index = max(composition(T/2 - period + 1:T/2, :)); 
            richness = find(diversity_index > abundance_threshold);
            diversity_record((mm-1) * ystep + pp, hh) = length(richness) / S(mm);
            
            % Record invasion success (if the invader species survived)
            invade_record((mm-1) * ystep + pp, hh) = max(composition(T - period + 1:T, S(mm) + 1)) > abundance_threshold;
            
            % Calculate fluctuation (coefficient of variation, CV) before invasion
            fluc_CV = std(composition(T/2 - 500:T/2 - 1, :)) ./ mean(composition(T/2 - 500:T/2 - 1, :));
            fluc_CV(isnan(fluc_CV)) = 0; % Replace NaN values with 0
            fluc_record((mm-1) * ystep + pp, hh) = mean(fluc_CV(fluc_CV > 0)); % Average CV across species
            
            % Measure the effect of invasion on community structure
            CC = min(mean(composition(T/2 - period + 1:T/2, 1:S(mm))), mean(composition(T - period + 1:T, 1:S(mm))));
            CC = length(find(CC > abundance_threshold)); % Number of species surviving after invasion
            SS = max(mean(composition(T/2 - period + 1:T/2, 1:S(mm))), mean(composition(T - period + 1:T, 1:S(mm))));
            SS = length(find(SS > abundance_threshold)); % Number of species before invasion
            invasion_effect((mm-1) * ystep + pp, hh) = 1 - CC / SS; % Relative change in species composition
        end
    end
end

% Function to simulate community dynamics with pH influence and species invasion
function [composition] = pH_compute_invasion(S, D, A, beta, delta, T, step)

    S = S + 1; % Adding the invader species to the system
    NN = zeros(T, S); % Species abundance over time
    r = ones(S, 1); % Growth rates (set to 1 for simplicity)
    p = zeros(T, 1); % pH dynamics array
    g = (rand(S, 1) - 0.5) * 2; % Growth rate modification by pH
    k = (rand(S, 1) - 0.5) * 2; % pH effect scaling factor for each species

    % Interaction matrix
    AA = rand(S, S) * A * 2; % Random interaction matrix with strength A
    AA = AA - diag(diag(AA)) + eye(S); % Ensure diagonal elements are 1 (self-interaction)

    % Initialize species abundances with small random values
    N0 = rand(1, S) * 0.01;
    N0(S) = 0; % Invader species starts with zero abundance
    NN(1, :) = N0;

    % Simulate dynamics for each time step
    for i = 2:T
        for j = 1:S
            % pH effect on growth rate
            pH_effect = g(j) * p(i-1);

            % Runge-Kutta integration (4th order) for species dynamics
            k1 = r(j) * NN(i-1, j) * (1 - AA(j, :) * (NN(i-1, :)') + pH_effect) * step;
            k2 = r(j) * (NN(i-1, j) + k1 / 2) * (1 - AA(j, :) * (NN(i-1, :)') + pH_effect) * step;
            k3 = r(j) * (NN(i-1, j) + k2 / 2) * (1 - AA(j, :) * (NN(i-1, :)') + pH_effect) * step;
            k4 = r(j) * (NN(i-1, j) + k3) * (1 - AA(j, :) * (NN(i-1, :)') + pH_effect) * step;
            NN(i, j) = NN(i-1, j) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
        end

        % Apply migration
        if i < T/2
            NN(i, 1:S-1) = NN(i, 1:S-1) + D; % Migration before invasion
        else
            NN(i, :) = NN(i, :) + D; % Migration after invasion
        end

        % Update pH dynamics
        p(i) = p(i-1) + step * (delta * (-p(i-1)) + beta * sum(k .* NN(i-1, :)'));
    end

    % Return final composition of species
    composition = NN;
end