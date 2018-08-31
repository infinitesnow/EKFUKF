function [mse_transient, mse_steadystate] = compute_pi(prediction, omega, step_length, t_transient)
omega = omega(step_length+1:end); % Skip first step
prediction = prediction(step_length+1:end); % Skip first step
n_steps = floor(length(omega)/step_length);
mse_transient = zeros(1,n_steps);
mse_steadystate = mse_transient;

for ii = 1:n_steps
    this_pred_omega = prediction(1+(ii-1)*step_length:ii*step_length);
    this_pred_omega_transient = this_pred_omega(1:t_transient(ii));
    this_pred_omega_steadystate = this_pred_omega(t_transient(ii)+1:end);
    this_omega = omega(ii*step_length+1);
    mse = sum((this_pred_omega-this_omega).^2);
    mse_transient(ii) = sum((this_pred_omega_transient-this_omega).^2);
    mse_steadystate(ii) = sum((this_pred_omega_steadystate-this_omega).^2);
    
end
mse_transient = mean(mse_transient);
mse_steadystate = mean(mse_steadystate);
end