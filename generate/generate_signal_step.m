function [signal, instantaneous_omega] = ...
    generate_signal_step(n_samples_step,initial_omega,steps,sigma_noise,sigma_omega_noise)
PLOT_SIGNAL = false;  

%% Generate omega profile
profile = zeros(1,length(steps)+1);
profile(1) = initial_omega;
ii = 2;
for step = steps
    profile(ii) = profile(ii-1) + step;
    ii = ii+1;
end

n_steps = length(profile);
signal=zeros(1,n_samples_step*n_steps);
instantaneous_omega=signal;
t = 1:n_samples_step*n_steps;

for ii=1:n_steps
    omega = profile(ii);
    instantaneous_omega(1+(ii-1)*n_samples_step:ii*n_samples_step) = ...
        omega*ones(1,n_samples_step);
end

Phi = cumsum(instantaneous_omega);

for ii=1:length(t)
    signal(ii) = sin(Phi(ii) + normrnd(0,sigma_omega_noise*instantaneous_omega(ii)));
end

if(sigma_noise~=0)
    t = 1:length(signal);
    noise = randn(1,length(signal));
    noise = noise + (sin(randi(100)/100*pi*t)).^2;
    noise = noise + (sin(randi(100)/100*pi*t)).^3;
    noise = noise + (sin(randi(100)/100*pi*t)).^4;
    noise = noise / max(noise);
    signal = signal + sigma_noise*noise;
    if (PLOT_SIGNAL)
        plot(noise);
        title('Generated signal')
        hold on;
        plot(signal);
        legend('Noise','Signal')
    end
end
end