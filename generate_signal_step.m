function [noisy_signal, instantaneous_omega] = generate_signal_step(n_periods_step,initial_omega,steps,sigma_error)
    profile = generate_profile(initial_omega,steps);  

    % We generate the input signal
    signal=[];
    instantaneous_omega=[];
    for omega = profile
        T=2*pi/omega;
        t = 0:1:(T*n_periods_step);
        signal = [signal sin(omega*t)]; % for each step, we generate a sinusoid with the corresponding frequency
        instantaneous_omega = [instantaneous_omega omega*ones(1,length(t))]; % Return the real frequency
    end
    noisy_signal = signal+normrnd(0,sigma_error,1,length(signal)); % We add noise
end 

function profile = generate_profile(initial_omega,steps)
    profile=initial_omega;
    for step=steps
        profile=[profile initial_omega+initial_omega*step initial_omega-initial_omega*step];
    end
end
