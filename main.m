clf
clc
clear all

%% Initial values
initial_omega = pi/8;
tolerance = 10;
n_simulations = 100;
step_profile = [0.1 0.3 0.5 0.8];

[mse_ekf, mse_ukf] = generate_predictions(n_simulations, initial_omega, step_profile);
ekf_accuracy = sum(mse_ekf<tolerance)/n_simulations
ukf_accuracy = sum(mse_ukf<tolerance)/n_simulations

function [mse_ekf, mse_ukf] = generate_predictions(n_simulations, initial_omega, step_profile)
mse_ekf = zeros(n_simulations, 1);
mse_ukf = zeros(n_simulations, 1);

for ii = 1:n_simulations
    %% Generate signal
    % Generate signal with a given frequency profile
    [signal, omega]=generate_signal_step(3,initial_omega,step_profile,0.01); 

    % Generate signal with harmonics
    % [signal omega]=generate_signal_harmonics(5,[initial_omega 2; initial_omega*4 0.01],0); 

    %% Track
    % We initialize filter with states to 0 and around the right initial frequency with a given variance
    % This value is also assumed as a known uncertainty in the initialization of P for the EKF and UKF
    sigma_init=0.5*initial_omega; 
    x_pred_0=[0 0 normrnd(initial_omega,sigma_init)];

    pred_vec_ekf=ekf( ...
        1, ...%lambda
        x_pred_0, ...%x_pred_0
        sigma_init, ...%sigma_init
        signal ...%signal
    );

    pred_vec_ukf=ukf( ...
        1, ...%alpha
        2, ...%beta
        2, ...%k
        1e-5,...%q
        1e-10,...%r
        1e-10, ...%sigma_init
        x_pred_0, ...%x_pred_0
        signal ...%signal
    );
    
    mse_ekf(ii) = sum((omega-pred_vec_ekf(3,:)).^2);
    mse_ukf(ii) = sum((omega-pred_vec_ukf(3,:)).^2);
end
end