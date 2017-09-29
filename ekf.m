clf
do_plot=1

%% Input signal
sigma=0.1

%% Generate signal with PI profile
% n_samples_step = 2000;
% n_periods_step=10;
% [original_signal,y,omega]=generate_signal_pi(n_samples_step,n_periods_step,sigma);


%% Generate signal with single frequency
% n_samples=10000
% n_periods=10
% true_pulsation=10

%% Generate system from simulation
x1_0=1
x2_0=0
x3_0=0.05
n_samples=1000
y=simulate_system(x1_0,x2_0,x3_0,0.001,0.01,n_samples)

%% System data
H=[1 0 0];
I3=[0 0 0
    0 0 0
    0 0 1];

%%% Initial values
simulation_length=length(y);
lambda=10
x_pred_0=[x1_0+0.3 x2_0-0.1 x3_0-0.06]';
K_0=[0 0 0]';
P_0=0.1*ones(3);

%% Predict
%%% Set initial values
r=1;
q=lambda*r;
x_pred=x_pred_0
K=K_0
P=P_0 
pred_vec=[];

%%% Initialize plot
figure(1)
subplot(2,1,1)
h=animatedline
j=animatedline
k=animatedline
l=animatedline
m=animatedline
subplot(2,1,1)
hold on
subplot(2,1,2)
hold on
for i = 1:simulation_length
    x_pred_model=f(x_pred); % Compute f in the state predicted by the Kalman filter in the previous instant
    e = y(i)-H*x_pred_model; % Compute the innovation
    x_pred=x_pred_model+K*e; % Compute the new Kalman prediction
    
    pred_vec = [pred_vec x_pred]; % Save new prediction

    F=compute_F(x_pred); % Compute F in the predicted state
    P_new=F*(P-K*H*P)*F'+q*I3; % Compute the new P
    P=P_new % Update P
    K=P*H'*inv(H*P*H'+r) % Compute K
    
    %%% Plot dynamically
    if do_plot==1
        subplot(2,1,2)
        addpoints(h,i,x3_0)
        addpoints(j,i,x_pred(3),'r-')
        subplot(2,1,2)
        addpoints(k,i,y(i),'k-o')
        addpoints(l,i,x_pred(1),'b--')
        addpoints(m,i,x_pred(2),'r--')
        pause(0.1)
        drawnow
    end
end

function x_pred_model = f(x) % Use the model to predict next state
    x_pred_model = [    cos(x(3))*x(1)-sin(x(3))*x(2)
                        sin(x(3))*x(1)+cos(x(3))*x(2)
                        x(3)                                  ];
end
function F = compute_F(x) % Compute F matrix in a given state 
    F=[     cos(x(3))   -sin(x(3))  -sin(x(3))*x(1)-cos(x(3))*x(2)
            sin(x(3))   cos(x(3))    cos(x(3))*x(1)-sin(x(3))*x(2)
            0           0            1                               ];
end

