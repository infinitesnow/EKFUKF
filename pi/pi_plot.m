PLOT_FIT = false;
path = strcat(pwd,'\data\pi\',namestring,'\');
if ~exist(path,'dir'), mkdir(path), end

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

pi_curve_ekf_et = pi_curve_ekf_converged(1,:);
pi_curve_ekf_ess = pi_curve_ekf_converged(2,:);
pi_curve_ekf_sigma = pi_curve_ekf_converged(3,:);
pi_curve_ukf_et = pi_curve_ukf_converged(1,:);
pi_curve_ukf_ess = pi_curve_ukf_converged(2,:);
pi_curve_ukf_sigma = pi_curve_ukf_converged(3,:);

pi_curves_figure = figure('visible','off');
title('PI curves')
hold on
if (PLOT_FIT)
    plot(pi_curve_ekf_ess,fit_model(B_ekf,pi_curve_ekf_ess),'r-');
    plot(pi_curve_ukf_ess,fit_model(B_ukf,pi_curve_ukf_ess),'b-');
else
    scatter(pi_curve_ekf_ess,pi_curve_ekf_et,'ro')
    scatter(pi_curve_ukf_ess,pi_curve_ukf_et,'bx')
end
xlabel('MSE steady state')
ylabel('MSE transient')
if(~PLOT_FIT)
    legend('EKF','UKF')
    saveas(pi_curves_figure, strcat(path,'pi_cloud_',namestring,'.png'), 'png');
else
    legend('EKF fit','UKF fit')
    saveas(pi_curves_figure, strcat(path,'pi_fit_',namestring,'.png'), 'png');
end

sigma_error_figure = figure('visible','off');
subplot(2,2,1)
semilogx(pi_curve_ekf_sigma,pi_curve_ekf_et);
title('EKF Transient error')
xlabel('Sigma')
ylabel('Error')
subplot(2,2,2)
semilogx(pi_curve_ekf_sigma,pi_curve_ekf_ess);
xlabel('Sigma')
ylabel('Error')
title('EKF Steady state error')

subplot(2,2,3)
semilogx(pi_curve_ukf_sigma,pi_curve_ukf_et);
xlabel('Sigma')
ylabel('Error')
title('UKF Transient error')
subplot(2,2,4)
semilogx(pi_curve_ukf_sigma,pi_curve_ukf_ess);
xlabel('Sigma')
ylabel('Error')
title('UKF Steady state error')

saveas(sigma_error_figure, strcat(path,'sigma_error_figure_',namestring,'.png'), 'png');

%% Plot convergence
k = floor(n_sigmas/10);
sigma = 10;
gaussian_filter=fspecial('gauss',[1, k], sigma);
derivative_gaussian_filter = gradient(gaussian_filter); 

convergence_figure = figure('visible','off');
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

saveas(convergence_figure, strcat(path,'convergence_',namestring,'.png'), 'png');