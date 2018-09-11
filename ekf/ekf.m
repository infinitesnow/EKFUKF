function pred_vec=ekf(signal,x_pred_0,initial_sigma,varargin)
    window_size = 50;
    SAVE_EKF_PLOT = false;
    VERBOSE = false;
        
    %% Parse arguments
     % only want 2 optional inputs at most
    numvarargs = length(varargin);
    if numvarargs > 2
        error('Too many inputs');
    end

    % Set defaults for optional inputs
    optargs = {...
        1e-3, ...%r
        1, ...%q
    };

    % Now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;

    % Place optional args in variable names
    [r,q] = optargs{:};
    
    %% Set initial values
    simulation_length=length(signal);
    
    x_pred=x_pred_0';
    K=[0 0 0]';
    P=initial_sigma*eye(3);
    
    H=[1 0 0];
    I3=eye(3);
    
    if (SAVE_EKF_PLOT)
        initialize_plot_ekf
    end
    
    %% Start tracking
    pred_vec=[];
    
    for ii=1:simulation_length
        y = signal(ii);
        e = y-H*compute_f(x_pred); % Compute the innovation
        x_pred=compute_f(x_pred)+K*e; % Compute the new Kalman prediction

        F=compute_F(x_pred); % Compute F in the predicted state
        P=F*P*F'+q*I3-F*P*H'*inv(H*P*H'+r)*H*P*F'; % Compute the new P
        K=P*H'*inv(H*P*H'+r); % Compute K

        pred_vec = [pred_vec x_pred]; % Save new prediction
        
        if (VERBOSE)
            fprintf('Iteration %d...\n',ii);
        end
        if (SAVE_EKF_PLOT && ii~=simulation_length)
            plot_ekf
        end
    end
    
    if (SAVE_EKF_PLOT)
        addpath('./generate/')
        generate_video('ekf',EKF);
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