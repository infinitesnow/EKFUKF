clf
clc
clear all

%% Initial values
initial_omega=pi/8;

%% Generate signal

% Generate signal with a given frequency, no noise
% [signal omega]=generate_signal_step(10,initial_omega,[],0.1); 

% Generate signal with a given frequency profile
[signal omega]=generate_signal_step(3,initial_omega,[0.1 0.3 0.5 0.8],0.1); 

% Generate signal with harmonics
% [signal omega]=generate_signal_harmonics(5,[initial_omega 2; initial_omega*4 0.01],0);

% Generate signal from simulation 
% [signal omega]=generate_signal_simulation(1,0,initial_omega,0.01,0.1,1000); 

%% Track
lambda=1;
% We initialize filter with states to 0 and around the right initial frequency with a given variance
sigma_init=1*initial_omega; % This value is also used for the initialization of P
x_pred_0=[0 0 initial_omega+normrnd(0,sigma_init)];

tic
% [pred_vec K_vec e_vec P_vec]=ekf(lambda,x_pred_0,sigma_init,signal);
pred_vec=ukf(1,2,2,0.1,0.01,x_pred_0,sigma_init,signal);
toc

%% Plot results
dynamic_plot(omega,pred_vec,0.001) % Plot dynamically only the frequency estimationnon 
% dynamic_plot_complete(signal,omega,pred_vec,K_vec,e_vec,P_vec) % Plot dynamically all values
% static_plot(signal,omega,pred_vec,K_vec,e_vec) % Plot statically