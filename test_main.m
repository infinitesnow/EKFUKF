clc
clear all

%% Initial values
initial_omega=pi/8;

%% Generate signal

% Generate signal with a given frequency, no noise
% [signal omega]=generate_signal_step(10,initial_omega,[],0.1); 

% Generate signal with a given frequency profile
[signal omega]=generate_signal_step(3,initial_omega,[0.1 0.3 0.5 0.8],0.01); 

% Generate signal with harmonics
% [signal omega]=generate_signal_harmonics(5,[initial_omega 2; initial_omega*4 0.01],0);

% Generate signal from simulation 
% [signal omega]=generate_signal_simulation(1,0,initial_omega,0.01,0.1,1000); 

%% Track
% We initialize filter with states to 0 and around the right initial frequency with a given variance
sigma_init=0.3*initial_omega; % This value is also used for the initialization of P for the EKF
x_pred_0=[0 0 normrnd(initial_omega,sigma_init)];

disp('Computing EKF estimation')
tic
pred_vec_ekf=ekf( ...
    1, ...%lambda
    x_pred_0, ...%x_pred_0
    sigma_init, ...%sigma_init
    signal ...%signal
);
toc

disp('Computing UKF estimation')
tic
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
toc

%% Plot statically
figure()
hold on
static_plot(signal,omega,pred_vec_ekf,pred_vec_ukf) 

function static_plot(signal,omega,pred_vec_ekf,pred_vec_ukf)
    t=1:length(signal);
    subplot(3,1,1)
    plot(t,omega,'k-',t,pred_vec_ekf(3,:),'r-',t,pred_vec_ukf(3,:),'b');
    legend('True','Estimated EKF','Estimated UKF');
    title('Omega estimation');
    subplot(3,1,2)
    plot(t,signal,'ko',t,pred_vec_ekf(1,:),'b-',t,pred_vec_ekf(2,:),'r-');
    title('States EKF');
    subplot(3,1,3)
    plot(t,signal,'ko',t,pred_vec_ukf(1,:),'b-',t,pred_vec_ukf(2,:),'r-');
    title('States UKF');
end