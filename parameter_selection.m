clc
clear all

%% Initial values
initial_omega = pi/8;
tolerance = 20;
n_simulations = 20;
step_profile = [0.1 0.3 0.5 0.8];

mse_ukf = zeros(n_simulations, 1);
ukf_accuracy = zeros([4,4]);

n_parameter=6;
offset=-1;
for ii=1:n_parameter
    for jj=1:n_parameter
            if (ii>jj)
                continue;
            end
            for ss = 1:n_simulations
                %% Generate signal
                % Generate signal with a given frequency profile
                [signal, omega]=generate_signal_step(3,initial_omega,step_profile,0.01); 

                % Generate signal with harmonics
                % [signal omega]=generate_signal_harmonics(5,[initial_omega 2; initial_omega*4 0.01],0); 

                %% Track
                % We initialize filter with states to 0 and around the right initial frequency with a given variance
                % This value is also assumed as a known uncertainty in the initialization of P for the EKF and UKF
                sigma_init=0.5*initial_omega; 
                x_pred_0=[0 0 normrnd(initial_omega,sigma_init)];

                pred_vec_ukf=ukf( ...
                    1, ...%alpha
                    2, ...%beta
                    2, ...%k
                    10^(-ii+offset),...%q
                    10^(-jj+offset),...%r
                    sigma_init, ...%sigma_init
                    x_pred_0, ...%x_pred_0
                    signal ...%signal
                );

                mse_ukf(ss) = sum((omega-pred_vec_ukf(3,:)).^2);
            end
        ukf_accuracy(ii,jj) = sum(mse_ukf<tolerance)/n_simulations;
        mse_ukf=zeros(n_simulations,1);
    end
end