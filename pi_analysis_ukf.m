function [pi_ukf, pred_omega] = pi_analysis_ukf(signal, omega, step_length, t_transient, r, q)
VERBOSE = true;

%% PI analysis
t_steadystate = step_length-t_transient;

%% Track    
sigma_init = 0.01;
x_pred_0 = [0, 0, omega(1)]; % Perfect initialization

pred_vec = ukf( ... 
    signal, ...%signal
    x_pred_0, ...%x_pred_0
    sigma_init, ...%sigma_init
    r, ... %r,
    q ... %q
);

% Take only the state we are interested in
pred_omega = pred_vec(3,:); 

%% Compute PI
[mse_transient, mse_steadystate] = compute_pi(pred_omega, omega, step_length, t_transient);
pi_ukf = [mse_transient, mse_steadystate];

%% Output values
if VERBOSE
    disp('************ UKF *************');
    fprintf('Sigma: %e \t( r=%e, q=%e)\n',r/q,r,q); % Not only dependent on sigma?
    fprintf('Transient: %e, SS: %e\n',mse_transient,mse_steadystate);
end

end