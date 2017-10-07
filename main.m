clf
clc
clear all

%% Initial values
initial_omega=pi/8;

%% Generate signal

% Generate signal with a given frequency, no noise
[signal omega]=generate_signal_step(10,initial_omega,[],0.1); 

% Generate signal with a given frequency profile
%[signal omega]=generate_signal_step(3,initial_omega,[0.1 0.3 0.5 0.8],0.1); 

% Generate signal from simulation 
%[signal omega]=generate_signal_simulation(1,0,initial_omega,0.01,0.1,1000); 

%% Track
lambda=100;
% We initialize filter with states to 0 and around the right initial frequency with a given variance
sigma_init=1*initial_omega; % This value is also used for the initialization of P
x_pred_0=[0 0 initial_omega+normrnd(0,sigma_init)];


[pred_vec K_vec e_vec P_vec]=ekf(lambda,x_pred_0,sigma_init,signal,omega);

%% Plot results
% dynamic_plot(signal,omega,pred_vec,K_vec,e_vec,P_vec) % Plot dynamically
static_plot(signal,omega,pred_vec,K_vec,e_vec) % Plot dynamically