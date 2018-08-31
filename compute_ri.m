function ri = compute_ri(prediction, omega, step_length,threshold)
omega = omega(step_length+1:end); % Skip first step
prediction = prediction(step_length+1:end); % Skip first step
n_steps = floor(length(omega)/step_length);
ri = n_steps;

%% Compute performance indexes

for ii = 1:n_steps
    this_pred_omega = prediction(1+(ii-1)*step_length:ii*step_length);
    this_omega = omega(ii*step_length+1);
    mse = sum((this_pred_omega-this_omega).^2);
    if mse < threshold
        ri = ii;
    end    
end

end