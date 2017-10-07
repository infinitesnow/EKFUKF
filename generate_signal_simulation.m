function [y omega] = generate_signal_simulation(x1_0,x2_0,x3_0,omega_variance,measurement_variance,simulation_length)
    y=[];
    omega=[];
    x1=x1_0;
    x2=x2_0;
    x3=x3_0;
    for t = 1:simulation_length
        x1_new=cos(x3)*x1 - sin(x3)*x2;
        x2_new=sin(x3)*x1 + cos(x3)*x2;
        x3_new=x3+normrnd(0,omega_variance);
        measurement=x1_new+normrnd(0,measurement_variance);

        % Update variables
        x1=x1_new; 
        x2=x2_new;
        x3=x3_new;

        % Store system states
        y=[y measurement];
        omega=[omega x3];
    end
end