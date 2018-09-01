%% RI
clear all
clc

PLOT_PREDICTION = false;

% Parameters
initial_omega = pi/2;
step_profile = [pi/13 -pi/7 pi/5 -pi/3 pi/2 -pi/1.7 pi/1.5]; % Approximation of Fig.2 on the paper
step_length = 400;    
sigma = 1e-4;
q = 1e-7;
r = sigma*q;
n_simulations = 20;
threshold = 50;
sigma_noise = 0.03;

ri_ekf = zeros(1,n_simulations);
ri_ukf = ri_ekf;
for ii = 1:n_simulations
    fprintf('Iteration %d\n',ii)
    % Generate signal
    [signal, omega]=generate_signal_step(step_length,initial_omega,step_profile,sigma_noise);

    %% Track    
    sigma_init = 0.01;
    x_pred_0 = [0, 0, omega(1)]; % Perfect initialization

    pred_vec_ekf = ekf( ... 
        signal, ...%signal
        x_pred_0, ...%x_pred_0
        sigma_init, ...%sigma_init
        r, ... %r,
        q ... %q
    );
    pred_vec_ukf = ukf( ... 
        signal, ...%signal
        x_pred_0, ...%x_pred_0
        sigma_init, ...%sigma_init
        r, ... %r,
        q ... %q
    );

    % Take only the state we are interested in
    pred_omega_ekf = pred_vec_ekf(3,:); 
    pred_omega_ukf = pred_vec_ukf(3,:); 
    
    disp('************ EKF *************');
    ri_ekf(ii) = compute_ri(pred_omega_ekf,omega,step_length,threshold);
    disp('************ UKF *************');
    ri_ukf(ii) = compute_ri(pred_omega_ukf,omega,step_length,threshold);

    % Plot ground truth and predictions
    if PLOT_PREDICTION
        figure(1)
        clf
        title('Predictions for RI');
        t = 1:length(omega);
        stairs(t, omega, 'k');
        xlabel('Samples')
        ylabel('Omega')
        set(gca,'YLim',[0,pi]);
        set(gca,'YTick',0:pi/4:pi);
        set(gca,'YTickLabel',{'0' 'pi/4','pi/2','3/4pi','pi'});
        hold on
        plot(t,pred_omega_ekf,'ro');
        plot(t,pred_omega_ukf,'bx');
        legend('True')%,'EKF','UKF')
    end
end

n_steps = length(step_profile);
steps_axis = 1:n_steps;
psi_ekf = hist(ri_ekf, steps_axis);
psi_ukf = hist(ri_ukf, steps_axis);

figure()
title('RI')
bar(steps_axis,[psi_ekf', psi_ukf'])
legend('EKF','UKF')