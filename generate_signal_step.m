function [signal, instantaneous_omega] = generate_signal_step(n_samples_step,initial_omega,steps,symmetric)
    
    %% Generate omega profile
    if (symmetric)
        profile=initial_omega;
        for step = steps
            profile=[profile initial_omega+initial_omega*step initial_omega-initial_omega*step];
        end
    else 
        profile = zeros(1,length(steps));
        profile(1) = initial_omega;
        jj = 2;
        for step = steps
            profile(jj) = profile(jj-1) + step;
            jj = jj+1;
        end
    end
    
    signal=[];
    instantaneous_omega=[];
    for omega = profile
        t = 0:n_samples_step;
        signal = [signal sin(omega*t)]; % For each step, we generate a sinusoid
        instantaneous_omega = ...
            [instantaneous_omega omega*ones(1,length(t))]; % In the end we 
                            % return the real omegas for testing purposes
    end
end