clear; clc; close all;

% Parameters
S = 2:2:60; % Number of species considered
DD = 1e-5; % Migration rate
A_mean_range = [0.0 1]; % Range for mean interaction strength (inhibition alpha)
ystep = 30; % Number of steps in the range for A_mean
T_real = 2000; % Total simulation time (real time)
step = 0.1; % Time step for the simulation
T = T_real / step; % Total number of steps
num_sim = 1000; % Number of simulations for each parameter set
period = 100; % Time period for calculating richness after perturbations
abundance_threshold = 8e-4; % Threshold for species survival based on abundance

% Preallocate result matrices
fluc_record = zeros(length(S) * ystep, num_sim); % Record if the resident community fluctuates
diversity_record = zeros(length(S) * ystep, num_sim); % Record pre-invasion diversity
invade_record = zeros(length(S) * ystep, num_sim); % Record invasion success
invasion_effect = zeros(length(S) * ystep, num_sim); % Record long-term effect of invasion on community

% Range of A_mean (interaction strength) values
A_mean = linspace(A_mean_range(1), A_mean_range(2), ystep);

% Main simulation loop
for mm = 1:length(S)
    for pp = 1:ystep
        for hh = 1:num_sim
            
            % Compute species composition and invasion success using Lotka-Volterra dynamics
            [composition, invade_ratio] = LV_compute_invasion(S(mm), DD, A_mean(pp), T, step);
            
            % Calculate diversity before invasion
            diversity_index = max(composition(T/2 - period + 1:T/2, :)); 
            richness = find(diversity_index > abundance_threshold);
            diversity_record((mm-1) * ystep + pp, hh) = length(richness) / S(mm);
            
            % Record invasion success ratio
            invade_record((mm-1) * ystep + pp, hh) = invade_ratio;
            
            % Calculate fluctuation (coefficient of variation, CV)
            fluc_CV = std(composition(T/2 - 500:T/2 - 1, :)) ./ mean(composition(T/2 - 500:T/2 - 1, :));
            fluc_CV(isnan(fluc_CV)) = 0;
            fluc_record((mm-1) * ystep + pp, hh) = mean(fluc_CV(fluc_CV > 0));
            
            % Measure the effect of invasion on community structure
            CC = min(mean(composition(T/2 - period + 1:T/2, 1:S(mm))), mean(composition(T - period + 1:T, 1:S(mm))));
            SS = max(mean(composition(T/2 - period + 1:T/2, 1:S(mm))), mean(composition(T - period + 1:T, 1:S(mm))));
            invasion_effect((mm-1) * ystep + pp, hh) = 1 - CC / SS;
            
        end
    end
end

% Function to simulate community dynamics and invasion success
function [composition, invade_ratio] = LV_compute_invasion(S, D, A, T, step)

    % Parameters
    period = 1000; % Time period for calculating richness
    abundance_threshold = 8e-4; % Threshold for survival species
    S = S + 1; % One extra species for the invader
    N_mean = 1 / (S * A * sqrt(pi / 2)); % Mean initial abundance
    sigma3 = N_mean / 12; % Standard deviation for initial abundance
    NN = zeros(T, S); % Matrix to store species abundance over time
    r = ones(S, 1); % Growth rates (set to 1 for simplicity)

    % Generate interaction matrix
    AA = rand(S, S)*A*2;
    AA = AA - diag(diag(AA)) + eye(S); % Ensure diagonal elements are 1
    
    % Initialize community composition (species abundances)
    N0 = abs(normrnd(N_mean, sigma3, 1, S)); 
    N0(S) = 0; % Initial invader abundance is zero
    NN(1, :) = N0;

    % Simulate resident community dynamics (before invasion)
    for i = 2:T/2
        for j = 1:S
            k1 = r(j) * NN(i-1, j) * (1 - AA(j, :) * (NN(i-1, :)')) * step;
            k2 = r(j) * (NN(i-1, j) + k1 / 2) * (1 - AA(j, :) * (NN(i-1, :)')) * step;
            k3 = r(j) * (NN(i-1, j) + k2 / 2) * (1 - AA(j, :) * (NN(i-1, :)')) * step;
            k4 = r(j) * (NN(i-1, j) + k3) * (1 - AA(j, :) * (NN(i-1, :)')) * step;
            NN(i, j) = NN(i-1, j) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
        end
        NN(i, 1:S-1) = NN(i, 1:S-1) + D; % Migration applied
    end

    % Invasion phase
    num_invade = 1; % Number of invasions
    invade_success = 0;
    for ppp = 1:num_invade
        % Set invader interaction matrix
        AA(S, :) = (rand(1, S) + 0.5) * A;  
        AA(:, S) = (rand(S, 1) + 0.5) * A;
        AA(S, S) = 1; 
        NN(T/2, S) = 1e-3; % Initial invader population

        % Simulate community dynamics with invader present
        for i = T/2+1:T
            for j = 1:S
                k1 = r(j) * NN(i-1, j) * (1 - AA(j, :) * (NN(i-1, :)')) * step;
                k2 = r(j) * (NN(i-1, j) + k1 / 2) * (1 - AA(j, :) * (NN(i-1, :)')) * step;
                k3 = r(j) * (NN(i-1, j) + k2 / 2) * (1 - AA(j, :) * (NN(i-1, :)')) * step;
                k4 = r(j) * (NN(i-1, j) + k3) * (1 - AA(j, :) * (NN(i-1, :)')) * step;
                NN(i, j) = NN(i-1, j) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
            end
            NN(i, :) = NN(i, :) + D; % Migration applied
        end
        
        % Check if the invader successfully establishes
        if max(NN(T-period+1:T, S)) > abundance_threshold
            invade_success = invade_success + 1;
        end
    end
    
    % Return final community composition and invasion success ratio
    composition = NN;
    invade_ratio = invade_success / num_invade;
end
