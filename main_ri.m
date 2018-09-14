%% RI
clear all
close all
clc
fclose('all');

global logpath;

addpath('./ukf')
addpath('./ekf')
addpath('./generate')
addpath('./ri')

PLOT_PREDICTION = false;
PLOT_PROFILE = false;
COMPUTE = true;

% Parameters
initial_omega = pi/2;
step_profile = [pi/11 -pi/5.5 pi/4 -pi/3.4 pi/2.8 -pi/2.5 pi/2 -pi/1.6 pi/1.4]; % Approximation of Fig.2 on the paper
step_length = 400;    
sigma = 1e+1;
q = 1e-5;
r = sigma*q;
n_simulations = 200;
initialization_noise_sigma = 0.001;
threshold = 50;
sigma_noise = 0.02;
sigma_noise_omega = 0.001;

if (PLOT_PREDICTION)
    ri_figure = figure('visible','off');
    set(ri_figure, 'Position', [0, 0, 720, 720],'Resize','off')
end

if (COMPUTE)
    ri_ekf = zeros(1,n_simulations);
    ri_ukf = ri_ekf;
    for ii = 1:n_simulations
        fprintf('***** Iteration %d *****\n',ii)
        % Generate signal
        [signal, omega]=generate_signal_step(step_length,initial_omega,step_profile,sigma_noise,sigma_noise_omega);

        %% Track    
        x_pred_0 = [0, 0, normrnd(omega(1),omega(1)*initialization_noise_sigma)]; % Initialization

        pred_vec_ekf = ekf( ... 
            signal, ...%signal
            x_pred_0, ...%x_pred_0
            initialization_noise_sigma, ...%sigma_init
            r, ... %r,
            q ... %q
        );
        pred_vec_ukf = ukf( ... 
            signal, ...%signal
            x_pred_0, ...%x_pred_0
            initialization_noise_sigma, ...%sigma_init
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
            set(0, 'currentfigure', ri_figure);
            clf;
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
            legend('True','EKF','UKF')
            RI(ii) = getframe;
        end
    end
    
    if (PLOT_PREDICTION)
        generate_video('ri',RI);
    end

    save('ri_ekf.mat','ri_ekf')
    save('ri_ukf.mat','ri_ukf')
else
    if exist('ri_ekf.mat','file') && exist('ri_ukf.mat','file')
        load('ri_ekf.mat','ri_ekf')
        load('ri_ukf.mat','ri_ukf')
    else
        error('No saved curves. Recompute.')
    end
end

n_steps = length(step_profile);
steps_axis = 1:n_steps;
psi_ekf = hist(ri_ekf, steps_axis);
psi_ukf = hist(ri_ukf, steps_axis);

if PLOT_PROFILE 
    figure(1)
    hold off
    stairs(omega)
    xlabel('Samples')
    ylabel('Omega')
    set(gca,'YLim',[initial_omega-pi/2,initial_omega+pi/2]);
    set(gca,'YTick',0:pi/4:pi);
    set(gca,'YTickLabel',{'0' 'pi/4','pi/2','3/4pi','pi'});
end

ri_figure = figure(2)
title('RI')
bar(steps_axis,[psi_ekf', psi_ukf'])
legend('EKF','UKF')
string = sprintf('\\ri\\ri_figure_q%1.2e_s%1.2e_sn%1.2e_sno%1.2e.png',q,sigma,sigma_noise,sigma_noise_omega);
saveas(ri_figure, strcat(pwd,string), 'png');