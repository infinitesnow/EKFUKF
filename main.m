clc
clear all

%% Initial values for PI
initial_omega = pi/4;
step_profile = [pi/20 -pi/10 pi/8 -pi/6.5 pi/5.5 -pi/5];
step_length = 1000;
t_transient = [100 200 300 400 500 500];
t_steadystate = step_length-t_transient;
symmetric = false; % Should the profile be symmetric?

%% Generate signal
[signal, omega]=generate_signal_step(step_length,initial_omega,step_profile,false); 

%% Track    
sigma_init = 0.01;
x_pred_0 = [0, 0, initial_omega]; % Perfect initialization

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
    1e-2,...%q
    1e-5,...%r
    sigma_init, ...%sigma_init
    x_pred_0, ...%x_pred_0
    signal ...%signal
);

% Take only the state we are interested in
pred_omega_ekf = pred_vec_ekf(3,:); 
pred_omega_ukf = pred_vec_ukf(3,:); 

%% Compute PI
n_steps = length(step_profile+1);    % We add the first omega which is not a step

% EKF
[mse_ekf_transient, mse_ekf_steadystate] = compute_pi(pred_omega_ekf, omega, step_length, t_transient)

% UKF
[mse_ukf_transient, mse_ukf_steadystate] = compute_pi(pred_omega_ukf, omega, step_length, t_transient)


%% Plot ground truth and predictions
figure(1)
title('Omega profile');
t = 1:length(omega);
stairs(t, omega, 'k');
xlabel('Samples')
ylabel('Omega')
set(gca,'YTick',pi/8:pi/16:3/8*pi);
set(gca,'YTickLabel',{'pi/8','3/8*pi','pi/4','5/8*pi'});
hold on
plot(t,pred_omega_ekf,'ro');
plot(t,pred_omega_ukf,'bx');



function [mse_transient, mse_steadystate] = compute_pi(prediction, omega, step_length, t_transient)
omega = omega(step_length+1:end); % Skip first step
prediction = prediction(step_length+1:end); % Skip first step
n_steps = floor(length(omega)/step_length)-1;
mse_transient = zeros(1,n_steps);
mse_steadystate = mse_transient;
for ii = 1:n_steps
    this_pred_omega = prediction(1+(ii-1)*step_length:ii*step_length);
    this_pred_omega_transient = this_pred_omega(1:t_transient(ii));
    this_pred_omega_steadystate = this_pred_omega(t_transient(ii)+1:end);
    this_omega = omega(ii*step_length+1);
    mse_transient(ii) = sum((this_pred_omega_transient-this_omega).^2);
    mse_steadystate(ii) = sum((this_pred_omega_steadystate-this_omega).^2);
end
mse_transient = mean(mse_transient);
mse_steadystate = mean(mse_steadystate);
end