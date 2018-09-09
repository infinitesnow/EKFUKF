function [pi_ekf, pred_omega] = compute_pi_ekf_aux(signal, omega, initial_omega, initial_sigma, step_length, t_transient, r, q)
VERBOSE = true;

%% PI analysis
t_steadystate = step_length-t_transient;

%% Track    
x_pred_0 = [0, 0, initial_omega];

pred_vec=ekf( ...
    signal, ...%signal
    x_pred_0, ...%x_pred_0
    initial_sigma, ... %initialization of P is done with the a-priori known sigma of initialization noise
    r, ... %r
    q ... %q
);

% Take only the state we are interested in
pred_omega = pred_vec(3,:); 

%% Compute PI
[mse_transient, mse_steadystate] = compute_pi(pred_omega, omega, step_length, t_transient, true);
pi_ekf = [mse_transient, mse_steadystate];

%% Output values
if VERBOSE
    disp('************ EKF *************');
    fprintf('Lambda: %e \t( r=%e, q=%e)\n',r/q,r,q);
    fprintf('Transient: %e, SS: %e\n',mse_transient,mse_steadystate);
end

end