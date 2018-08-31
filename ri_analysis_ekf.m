function [ri_ekf, pred_omega] = ri_analysis_ekf(signal, omega, lambda)
VERBOSE = true;

%% RI analysis
%% Track    
sigma_init = 0.01;
x_pred_0 = [0, 0, omega(1)]; % Perfect initialization

pred_vec=ekf( ...
    signal, ...%signal
    x_pred_0, ...%x_pred_0
    sigma_init, ...%sigma_init
    lambda ... %lambda
);

% Take only the state we are interested in
pred_omega = pred_vec(3,:); 

%% Compute RI

%% Output values
if VERBOSE
    disp('************ EKF *************');
    fprintf('Lambda: %e\n',lambda);
    fprintf('');
end

end