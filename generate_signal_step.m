function [signal, instantaneous_omega] = generate_signal_step(n_samples_step,initial_omega,steps,varargin)
PLOT_SIGNAL = false;

%% Parse arguments
% only want 2 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 2
    error('Too many inputs');
end

% Set defaults for optional inputs
optargs = {...
    0, ...%sigma_noise
    false ...%symmetric
};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in variable names
[sigma_noise,symmetric] = optargs{:};    


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

sigma_omega = initial_omega * sigma_noise;

signal=[];
instantaneous_omega=[];
for omega = profile
    t = 0:n_samples_step;
    signal = [signal sin(normrnd(omega,sigma_omega)*t)]; % For each step, we generate a sinusoid
    instantaneous_omega = ...
        [instantaneous_omega omega*ones(1,length(t))]; % In the end we 
                        % return the real omegas for testing purposes
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