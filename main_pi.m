clear all
clc

q = 1e-7;
n_sigmas = 500;
sigma_start = 3;
sigma_end = 0;
initialization_noise_sigma = 0.001;
n_iterations = 2;
convergence_threshold = 10;
COMPUTE = false;
LOGSPACE = true;
PLOT_CLOUD = true;

%% PI
% Parameters
initial_omega = pi/4;
step_profile = [pi/20 -pi/10 pi/8 -pi/6.5 pi/5.5 -pi/5];
step_length = 1000;
t_transient = [100 200 300 400 500 500];
sigma_noise = 0.001;

% Generate signal
[signal, omega]=generate_signal_step(step_length,initial_omega,step_profile,sigma_noise);

% Perform analysis
if COMPUTE
    if LOGSPACE==true
        sigmas = logspace(sigma_start,sigma_end,n_sigmas);
    else
        sigmas = linspace(sigma_start,sigma_end,n_sigmas);
    end
    
    [pi_curve_ekf, pi_curve_ukf] = ...
        pi_analysis(signal, omega, step_length, t_transient, q, initialization_noise_sigma, sigmas, n_iterations);
else    
    if exist('pi_curve_ekf.mat','file') && exist('pi_curve_ukf.mat','file')
        load('pi_curve_ekf.mat','pi_curve_ekf')
        load('pi_curve_ukf.mat','pi_curve_ukf')
    else
        error('No saved curves. Recompute.')
    end
end

%% Remove iterations which didn't converge
not_converged_ekf_indices = union( ...
    find(pi_curve_ekf(1,:)>convergence_threshold), find(pi_curve_ekf(2,:)>convergence_threshold) );
converged_ekf_indices = setdiff(1:size(pi_curve_ekf,2),not_converged_ekf_indices);
fprintf('EKF didn''t converge %d times\n',length(not_converged_ekf_indices))
not_converged_ukf_indices = union( ...
    find(pi_curve_ukf(1,:)>convergence_threshold), find(pi_curve_ukf(2,:)>convergence_threshold) );
converged_ukf_indices = setdiff(1:size(pi_curve_ekf,2),not_converged_ukf_indices);
fprintf('UKF didn''t converge %d times\n',length(not_converged_ukf_indices))

converged_ekf_bool = zeros(1,size(pi_curve_ekf,2));
converged_ukf_bool = zeros(1,size(pi_curve_ukf,2));
converged_ekf_bool(converged_ekf_indices) = 1;
converged_ukf_bool(converged_ukf_indices) = 1;

pi_curve_ekf_converged = pi_curve_ekf(:,converged_ekf_indices);
pi_curve_ukf_converged = pi_curve_ukf(:,converged_ukf_indices);


%% LS regression

fit_model = @(b,x) b(1) + b(2)./(x + b(3)); % Generalised Hyperbola

ls_fit = @(b,pi_curve) sum((pi_curve(1,:)-fit_model(b,pi_curve(2,:))).^2);
ls_fit_ekf = @(b) ls_fit(b,pi_curve_ekf_converged);
ls_fit_ukf = @(b) ls_fit(b,pi_curve_ukf_converged);

B0 = [0 0 0];
opt = optimset('TolFun',1e-8,'TolX',1e-8);
B_ekf = fminsearch(ls_fit_ekf, B0, opt);
B_ukf = fminsearch(ls_fit_ukf, B0, opt);    

figure(2)
title('PI curves')
pi_curve_ekf_et = pi_curve_ekf_converged(1,:);
pi_curve_ekf_ess = pi_curve_ekf_converged(2,:);
pi_curve_ekf_sigma = pi_curve_ekf_converged(3,:);
pi_curve_ukf_et = pi_curve_ukf_converged(1,:);
pi_curve_ukf_ess = pi_curve_ukf_converged(2,:);
pi_curve_ukf_sigma = pi_curve_ukf_converged(3,:);
hold on
plot(pi_curve_ekf_ess,fit_model(B_ekf,pi_curve_ekf_ess),'r-');
plot(pi_curve_ukf_ess,fit_model(B_ukf,pi_curve_ukf_ess),'b-');
if (PLOT_CLOUD)
    scatter(pi_curve_ekf_ess,pi_curve_ekf_et,'ro')
    scatter(pi_curve_ukf_ess,pi_curve_ukf_et,'bx')
    legend('EKF','EKF fit','UKF','UKF fit')
else
    legend('EKF fit','UKF fit')
end
xlabel('MSE steady state')
ylabel('MSE transient')

figure(3)
subplot(1,2,1)
semilogx(pi_curve_ekf_sigma,pi_curve_ekf_et);
title('EKF Transient error')
xlabel('Sigma')
ylabel('Error')
subplot(1,2,2)
semilogx(pi_curve_ekf_sigma,pi_curve_ekf_ess);
xlabel('Sigma')
ylabel('Error')
title('EKF Steady state error')

figure(4)
subplot(1,2,1)
semilogx(pi_curve_ukf_sigma,pi_curve_ukf_et);
xlabel('Sigma')
ylabel('Error')
title('UKF Transient error')
subplot(1,2,2)
semilogx(pi_curve_ukf_sigma,pi_curve_ukf_ess);
xlabel('Sigma')
ylabel('Error')
title('UKF Steady state error')

%% Plot convergence
k = 20;
sigma = 10;
gaussian_filter=fspecial('gauss',[1, k], sigma);
derivative_gaussian_filter = gradient(gaussian_filter); 

figure(5)
subplot(1,2,1)
convergence_ekf = conv(cumsum(converged_ekf_bool),derivative_gaussian_filter);
convergence_ekf = convergence_ekf(1:end-k+1);
semilogx(pi_curve_ekf(3,:),convergence_ekf,'b')
xlim([0 +inf]) 
title('Convergence ratio EKF')
subplot(1,2,2)
convergence_ukf = conv(cumsum(converged_ukf_bool),derivative_gaussian_filter);
convergence_ukf = convergence_ukf(1:end-k+1);
semilogx(pi_curve_ukf(3,:),convergence_ukf,'b')
xlim([0 +inf])
title('Convergence ratio UKF')