function ri = compute_ri(prediction, omega, step_length,threshold)
VERBOSE = false;

omega = omega(step_length+1:end); % Skip first step
prediction = prediction(step_length+1:end); % Skip first step
n_steps = floor(length(omega)/step_length);
ri = n_steps;

%% Compute performance indexes

for ii = 1:n_steps
    this_pred_omega = prediction(1+(ii-1)*step_length:ii*step_length);
    this_omega = omega(ii*step_length+1);
    mse = sum((this_pred_omega-this_omega).^2);
    if (VERBOSE)
        fprintf('MSE %e ',mse);
    end
    if mse < threshold
        ri = ii;
        if(VERBOSE)
            fprintf('OK, updating RI to %d\n',ri);
        end
    else
        if(VERBOSE)
            fprintf('Not updating RI at step %d\n',ii);
        end
        break
    end
end

end