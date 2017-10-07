function [pred_vec K_vec e_vec P_vec]=ekf(lambda,x_pred_0,sigma_init,signal,omega)

    %% Set initial values
    simulation_length=length(signal);
    
    r=1;
    q=lambda*r;
    
    x_pred=x_pred_0';
    K=[0 0 0]';
    P=sigma_init*eye(3);
    
    H=[1 0 0];
    I3=[1 0 0
        0 1 0
        0 0 1];
    
    %% Start tracking
    pred_vec=[];
    K_vec=[];
    e_vec=[];
    P_vec=[];
    for y=signal
        e = y-H*compute_f(x_pred); % Compute the innovation
        x_pred=compute_f(x_pred)+K*e; % Compute the new Kalman prediction

        F=compute_F(x_pred); % Compute F in the predicted state
        P=F*P*F'+q*I3-F*P*H'*inv(H*P*H'+r)*H*P*F'; % Compute the new P
        K=P*H'*inv(H*P*H'+r); % Compute K

        pred_vec = [pred_vec x_pred]; % Save new prediction
        K_vec = [K_vec K];
        e_vec = [e_vec e];
        P_vec = [P_vec, P];
    end
end

function x_pred_model = compute_f(x) % Use the model to predict next state
    x_pred_model = [    cos(x(3))*x(1)-sin(x(3))*x(2)
                        sin(x(3))*x(1)+cos(x(3))*x(2)
                        x(3)                            ]; 
end
function F = compute_F(x) % Compute F matrix in a given state 
    F=[     cos(x(3))   -sin(x(3))  -sin(x(3))*x(1)-cos(x(3))*x(2)
            sin(x(3))   cos(x(3))    cos(x(3))*x(1)-sin(x(3))*x(2)
            0           0            1                               ];
end