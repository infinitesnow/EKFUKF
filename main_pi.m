clear all

q = 1e-7;
n_sigmas = 100;
sigma_start = +3;
sigma_end = -8;
COMPUTE = true;
PLOT_PREDICTIONS = false;

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
sigmas = logspace(sigma_start,sigma_end,n_sigmas);
pi_curve_ekf = zeros(2,length(sigmas));
pi_curve_ukf = pi_curve_ekf;
ii = 1;
if COMPUTE
    for sigma = sigmas
        fprintf('Iteration %d\n',ii);
        [pi_ekf, pred_omega_ekf] = pi_analysis_ekf(signal, omega, step_length, t_transient, q*sigma, q);
        pi_curve_ekf(:,ii) = pi_ekf';
        [pi_ukf, pred_omega_ukf] = pi_analysis_ukf(signal, omega, step_length, t_transient, q*sigma, q);
        pi_curve_ukf(:,ii) = pi_ukf';
        ii=ii+1;
        fprintf('\n');

        % Plot ground truth and predictions
        if PLOT_PREDICTIONS
            figure(1)
            title('Predictions for PI');
            t = 1:length(omega);
            stairs(t, omega, 'k');
            xlabel('Samples')
            ylabel('Omega')
            set(gca,'YLim',[pi/8,3/8*pi]);
            set(gca,'YTick',pi/8:pi/16:3/8*pi);
            set(gca,'YTickLabel',{'pi/8','3/16*pi','pi/4','5/16*pi','3/8pi'});
            hold on
            plot(t,pred_omega_ekf,'ro');
            plot(t,pred_omega_ukf,'bx');
            legend('True','EKF','UKF')
        end
    end

    save('pi_curve_ekf.mat','pi_curve_ekf')
    save('pi_curve_ukf.mat','pi_curve_ukf')
else
    if exist('pi_curve_ekf.mat','file') && exist('pi_curve_ukf.mat','file')
        load('pi_curve_ekf.mat','pi_curve_ekf')
        load('pi_curve_ukf.mat','pi_curve_ukf')
    else
        error('No saved curves. Recompute.')
    end
end

figure(2)
title('PI curves')
scatter(pi_curve_ekf(2,:),pi_curve_ekf(1,:),'ro')
hold on
scatter(pi_curve_ukf(2,:),pi_curve_ukf(1,:),'bx')
legend('EKF','UKF')
xlabel('MSE steady state')
ylabel('MSE transient')