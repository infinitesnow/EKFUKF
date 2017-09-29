function [signal, noisy_signal, instantaneous_omega] = generate_signal_pi(n_samples_step,n_periods_step,sigma)
    do_plot=0
    
    t = linspace(1,(2*n_periods_step+1)*pi,n_samples_step);

    % We create the test profile for the Performance Index
    steps=[pi/20 -pi/10 pi/8 -pi/6.5 pi/5.5 -pi/5];

    omega = pi/4;
    for i = 1:length(steps)
        omega=[omega omega(i)+steps(i)];
    end

    % We generate the input signal
    signal=[];
    instantaneous_omega=[];
    for i = 1:length(omega)
        signal = [signal sin(omega(i)*t)]; % for each step, we generate a sinusoid with the corresponding frequency
        instantaneous_omega = [instantaneous_omega omega(i)*ones(1,length(t))]; % Return the real frequency
    end
    noisy_signal = signal+normrnd(0,sigma,1,length(signal)); % We add noise
    
    
    if do_plot==1
        % We plot the signal
        clf
        figure(1)
        x=1:length(signal);
        plot(x,signal,'b--',x,noisy_signal,'r-')
        title('Signal')
        xlabel('t')

        % We plot the test profile
        figure(2)
        stairs(1:length(omega)+1,[omega omega(length(omega))]) % We duplicate the last sample to reproduce the plot in the paper
        ylim([pi/8 3/8*pi])
        title('Omega (PI)')
        ylabel('Omega')
        xlabel('ith step')
    end
end 

    