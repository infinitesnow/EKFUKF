clf

% Input signal
n_samples_step = 2000;
n_periods_step=10;
[original_signal,y,omega]=generate_signal_pi(n_samples_step,n_periods_step,0.001);
simulation_length=length(y);

% System data
H=[1 0 0];
I3=[0 0 0
    0 0 0
    0 0 1];

% Initial values
lambda=1
x_pred=[0 0 0]';
K=[0 0 0]';
P=0.1*eye(1);

r=1
q=lambda*r

pred_vec=[];
for i = 1:simulation_length
    x_pred_model=f(x_pred); % Compute f in the state predicted by the Kalman filter in the previous instant
    e = y(i)-H*x_pred_model; % Compute the innovation
    x_pred=x_pred_model+K*e; % Compute the new Kalman prediction
    
    pred_vec = [pred_vec x_pred]; % Save new prediction

    F=compute_F(x_pred); % Compute F in the predicted state
    P_new=F*(P-K*H*P)*F'+q*I3; % Compute the new P
    P=P_new; % Update P
    K=P*H'*inv(H*P*H'+r); % Compute K
end


figure(1)
subplot(3,1,1)
hold on
stairs(omega)
plot(pred_vec(3,:))
subplot(3,1,2)
hold on
plot(pred_vec(1,:))
plot(pred_vec(2,:))
subplot(3,1,3)
hold on
plot(y)
plot(original_signal,'r--')

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

