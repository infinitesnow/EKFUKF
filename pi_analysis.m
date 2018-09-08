function [pi_curve_ekf, pi_curve_ukf] = pi_analysis(signal, omega, step_length, t_transient, q, initialization_noise_sigma, sigmas, n_iterations)
n_sigmas = length(sigmas);
pi_curve_ekf = zeros(3,n_sigmas*n_iterations);
pi_curve_ukf = pi_curve_ekf;

for ii = 0:n_iterations-1
    fprintf('\n\n\n********** Iteration %d **********\n',ii+1);
    [   pi_curve_ekf(:,1+ii*n_sigmas : (ii+1)*n_sigmas), ...
        pi_curve_ukf(:,1+ii*n_sigmas : (ii+1)*n_sigmas) ...
        ] = pi_analysis_aux(signal,omega,step_length,t_transient,q,initialization_noise_sigma,sigmas);
end

pi_curve_ekf = sortrows(pi_curve_ekf',3)';
pi_curve_ukf = sortrows(pi_curve_ukf',3)';

save('pi_curve_ekf.mat','pi_curve_ekf')
save('pi_curve_ukf.mat','pi_curve_ukf')
end

function [pi_curve_ekf, pi_curve_ukf] = pi_analysis_aux(signal, omega, step_length, t_transient, q, initialization_noise_sigma, sigmas)
PLOT_PREDICTIONS = false;
pi_curve_ekf = zeros(3,length(sigmas));
pi_curve_ukf = pi_curve_ekf;
ii = 1;
initial_omega = normrnd(omega(1),initialization_noise_sigma*omega(1));
initial_sigma = initialization_noise_sigma*omega(1); %assume sigma is known a priori
fprintf('Initializing with initial omega %f, real omega %f (init. sigma %e)\n\n'...
    ,initial_omega,omega(1),initialization_noise_sigma);

for sigma = sigmas
    fprintf('Evaluating sigma no. %d...\n',ii);
    
    [pi_ekf, pred_omega_ekf] = compute_pi_ekf_aux( ...
        signal, omega, initial_omega, initial_sigma, ...
        step_length, t_transient, q*sigma, q);
    pi_curve_ekf(:,ii) = [ pi_ekf'; sigma ];
    
    [pi_ukf, pred_omega_ukf] = compute_pi_ukf_aux( ...
        signal, omega, initial_omega, initial_sigma, ...
        step_length, t_transient, q*sigma, q);
    pi_curve_ukf(:,ii) = [ pi_ukf'; sigma ];
    
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
end