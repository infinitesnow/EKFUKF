clear all
close all
clc
fclose('all');

global logpath;

addpath('./ukf')
addpath('./ekf')
addpath('./pi')
addpath('./generate')

n_sigmas = 150;
sigma_start = 3;
sigma_end = -4;
sigma_noises = [ 0 0.005 0.08 0.17 0.3];
qs = [ 1e-3 1e-4 10^(-4.5) 1e-5 10^(-5.5) 1e-6 10^(-6.5) 1e-7 10^(-7.5) 1e-8 1e-9];
sigma_omega_noise = 0;
initialization_noise_sigma = 0.001;
n_iterations = 1;
convergence_threshold = 10;
LOGSPACE = true;

%% PI
% Parameters
initial_omega = pi/4;
step_profile = [pi/20 -pi/10 pi/8 -pi/6.5 pi/5.5 -pi/5];
step_length = 1000;
t_transient = [100 200 300 400 500 500];

% Perform analysis
if LOGSPACE==true
    sigmas = logspace(sigma_start,sigma_end,n_sigmas);
else
    sigmas = linspace(sigma_start,sigma_end,n_sigmas);
end

for sigma_noise = sigma_noises
    for q = qs
        % Create folder
        namestring = sprintf ('sn_%1.3f_q%1.2e',sigma_noise,q);
        logpath = strcat(pwd,'\data\pi\',namestring,'\');
        if ~exist(logpath,'dir'), mkdir(logpath), end
        
        % Generate signal
        [signal, omega] = ...
            generate_signal_step( ...
            step_length,initial_omega,step_profile,sigma_noise,sigma_omega_noise);

        % Compute analysis
        [pi_curve_ekf, pi_curve_ukf] = ...
            pi_analysis(signal, omega, ...
            step_length, t_transient, q, initialization_noise_sigma, sigmas, n_iterations);
        
        % Save results
        save(strcat(logpath,'pi_curve_ekf.mat'),'pi_curve_ekf')
        save(strcat(logpath,'pi_curve_ukf.mat'),'pi_curve_ukf')
        
        pi_plot
    end
end
