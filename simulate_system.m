function out = simulate_system(x1_0,x2_0,x3_0,q,r,simulation_length)
    do_plot=0;

    %% Create animated vector plot
    if do_plot==1 
        figure(1)
        hold on
        title('System simulation')
        xlim([-1 1]) % Fit plot to the unit circle
        ylim([-1 1])
        pbaspect([1 1 1]) % 1:1 aspect ratio
    end

    %% Simulate
    out=[];
    freq=[];
    x1=x1_0;
    x2=x2_0;
    x3=x3_0;
    for t = 1:simulation_length
        % x3=omega(t-1) % Use frequencies from profile vector
        x1_new=cos(x3)*x1 - sin(x3)*x2;
        x2_new=sin(x3)*x1 + cos(x3)*x2;
        x3_new=x3+normrnd(0,q);
        y=x1_new+normrnd(0,r);

        % Update variables
        x1=x1_new; 
        x2=x2_new;
        x3=x3_new;

        % Store system states
        out=[out y];
        freq=[freq x3];

        % Plot
        if do_plot==1
            quiver(0,0,x1,x2)
            drawnow
            pause(0.01)
            cla
        end
    end

    %% Plot signal
    if do_plot==1 
        figure(2)
        subplot(2,1,1)
        plot(2:simulation_length,out)
        title('Simulated output measurement')
        subplot(2,1,2)
        plot(2:simulation_length,freq)
        hline=refline([0 x3_0])
        hline.Color='k'
        title('Actual frequency values')
        ylabel('rad/s')
    end
end