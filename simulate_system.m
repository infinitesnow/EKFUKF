% Initial conditions
x1_0=1;
x2_0=0;
x3_0=1/10;
epsilon=0;
do_plot=0;

% Noise variances
q=0.001;
r=0.01;

simulation_length = 1000;

% We create a time varying frequency profile
omega = [1/10*ones(1,simulation_length/4)...
    1/5*ones(1,simulation_length/4)...
    1/2*ones(1,simulation_length/4)...
    1/3*ones(1,simulation_length/4) ];

% Create animated plot
figure(1)
hold on
title('System simulation')
xlim([-1 1]) % Fit plot to the unit circle
ylim([-1 1])
pbaspect([1 1 1]) % 1:1 aspect ratio

% Simulate
out=[];
freq=[];
x1=x1_0;
x2=x2_0;
x3=x3_0;
for t = 2:simulation_length
    % x3=omega(t-1) % Use frequencies from our vector
    x1_new=cos(x3)*x1 - sin(x3)*x2;
    x2_new=sin(x3)*x1 + cos(x3)*x2;
    x3_new=(1-epsilon)*x3+normrnd(0,q);
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

% Plot signal
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