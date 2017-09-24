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
q=0.005;
r=0.01;
epsilon=0;
x_pred=[1/2 sqrt(2)/2 pi/4]';
K=[0.01 0.01 0.01]';
P=1*eye(1);

pred_vec=[];
for i = 1:simulation_length
    x_pred_model = [    cos(x_pred(3))*x_pred(1)-sin(x_pred(3))*x_pred(2)
                        sin(x_pred(3))*x_pred(1)+cos(x_pred(3))*x_pred(2)
                        (1-epsilon)*x_pred(3)                                  ]; % Use the model to predict next state
    e = y(i)-H*x_pred_model; % Compute the innovation
    x_pred=x_pred_model + K*e; % Compute the new Kalman prediction
    pred_vec = [pred_vec x_pred]; % Add new prediction of omega

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

function F = compute_F(x)
    F=[     cos(x(3))   -sin(x(3))  -sin(x(3))*x(1)-cos(x(3))*x(2)
            sin(x(3))   cos(x(3))    cos(x(3))*x(1)-sin(x(3))*x(2)
            0           0            1                               ];
end

