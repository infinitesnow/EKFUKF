clc
clear all
close all

addpath('./ukf/');
addpath('./ekf/');
addpath('./pi/');

%% Read signal from audio
freq = 1760;          % in hertz, known a priori
[ signal, sr ] = audioread(strcat('audio/',num2str(freq),'hz.wav'));
% Crop signal
n_samples = 500;      % samples
signal = signal(1:n_samples);
% Compute correct pulsation for verification
freq = 1/sr*freq;    % in samples/sec
initial_omega = 2*pi*freq;       % in rad/sec
omega = ones(1,length(signal))*initial_omega;
q = 10^(-11);
sigma = 10^(9);
t_transient = 250;

%% Track
% We initialize filter with states to 0 and around the right initial frequency with a given variance
sigma_init=0.003; % This value is also used for the initialization of P for the EKF
x_pred_0=[0 0 normrnd(initial_omega,sigma_init*initial_omega)];

disp('Computing EKF estimation')
tic
pred_vec_ekf = ekf( ... 
    signal, ...%signal
    x_pred_0, ...%x_pred_0
    sigma_init, ...%sigma_init
    q*sigma, ... %r,
    q ... %q
);
toc

disp('Computing UKF estimation')
tic
pred_vec_ukf = ukf( ... 
    signal, ...%signal
    x_pred_0, ...%x_pred_0
    sigma_init, ...%sigma_init
    q*sigma, ... %r,
    q ... %q
);
toc

%% Plot statically
figure()
hold on
t=1:length(signal);
subplot(3,1,1)
plot(t,omega,'k-',t,pred_vec_ekf(3,:),'r-',t,pred_vec_ukf(3,:),'b');
legend(['True (',num2str(omega(1)),')'],'Estimated EKF','Estimated UKF');
title('Omega estimation');
subplot(3,1,2)
plot(t,signal,'ko',t,pred_vec_ekf(1,:),'b-',t,pred_vec_ekf(2,:),'r-');
legend('signal','1','2')
title('States EKF');
subplot(3,1,3)
plot(t,signal,'ko',t,pred_vec_ukf(1,:),'b-',t,pred_vec_ukf(2,:),'r-');
legend('signal','1','2')
title('States UKF');

[pi_ekf(1), pi_ekf(2)] = compute_pi(pred_vec_ekf(3,:),omega,n_samples,t_transient,false);
[pi_ukf(1), pi_ukf(2)] = compute_pi(pred_vec_ukf(3,:),omega,n_samples,t_transient,false);
fprintf('PI for EKF:\n e_t: %e, e_ss: %e \n',pi_ekf(1),pi_ekf(2));
fprintf('PI for UKF:\n e_t: %e, e_ss: %e \n',pi_ukf(1),pi_ukf(2));