function [noisy_signal, instantaneous_omega] = generate_signal_harmonics(n_periods,harmonics,sigma_error)
    main_harmonic=harmonics(1,1);
    [n_harmonics,~]=size(harmonics);
    T=2*pi/main_harmonic;
    t = 0:1:(T*n_periods);
    simulation_length=length(t);
    signal = zeros(1,simulation_length);
    for ii = 1:n_harmonics
        signal = signal+harmonics(ii,2)*sin(harmonics(ii,1)*t);
    end
    noisy_signal = signal+normrnd(0,sigma_error,1,simulation_length); % We add noise
    instantaneous_omega=main_harmonic*ones(simulation_length);
end 