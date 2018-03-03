function [ pred_vec K_vec e_vec P_vec ]= ukf(alpha,beta,k,q,r,x_pred_0,sigma_init,signal)
    %% Set initial values
    L=3;
    simulation_length=length(signal);
    
    lambda=alpha^2*(L+k)-L;
    gamma=sqrt(L+lambda);
    
    I3=eye(3);
    Q=q*I3;
    
    x_pred=x_pred_0';
    P=sigma_init*eye(L);
        
    %% Compute weights
    % Weights to compute the mean of the transformed sigma points
    Wm0=lambda/(L+lambda)
    Wmi=1/(2*(L+lambda))*ones(1,2*L);
    Wm=[Wm0, Wmi];
    
    % Weights to compute the covariance of the transformed sigma points
    Wc0=lambda/(L+lambda)+1-alpha^2+beta;
    Wci=Wmi;
    Wc=[Wc0, Wci]; 
    
    % Check that dimensions are valid (both row, columns)
    assert(all(size(Wm)==[1 2*L+1])); 
    assert(all(size(Wc)==[1 2*L+1]));
    
    pred_vec=[]
    K_vec=[]
    e_vec=[]
    P_vec=[]
    
    %% Run filter
    for k=1:simulation_length
        
        sqrtP=chol(P,'lower'); % The matrix square root of P
                
        %% Compute the sigma points
        Sigma_points_p = zeros(L,L);
        Sigma_points_m = zeros(L,L);
        for ii=1:L
            Sigma_points_p(:,ii) = x_pred + gamma*sqrtP(:,ii);
            Sigma_points_m(:,ii) = x_pred - gamma*sqrtP(:,ii);
        end
        Sigma_points = [x_pred Sigma_points_p Sigma_points_m];
        n_sigma_points = size(Sigma_points,2);
        assert(n_sigma_points==2*L+1);
        
        %% Time update
        % Propagate the sigma points through the model
        transformed_sigma_points_model=zeros(L,n_sigma_points);
        for ii = 1:n_sigma_points
            transformed_sigma_points_model(:,ii)=compute_f(Sigma_points(:,ii));
        end
        
        % Compute the transformed sigma points mean
        mean_sigma_points_model=transformed_sigma_points_model*Wm';
        
        % Compute the transformed sigma points covariance
        Covariance_sigma_points_model=zeros(L);
        for ii=1:length(Wc)
            deviation_ii = transformed_sigma_points_model(:,ii)-mean_sigma_points_model;
            Covariance_sigma_points_model=Covariance_sigma_points_model+Wc(ii)*(deviation_ii)*(deviation_ii)';
        end
        Covariance_sigma_points_model = Covariance_sigma_points_model + Q;
        
        %% Measurements update
        
        % Propagate the sigma points through the measurement function
        % In our case, we have no measurement function, so h is the
        % identity function. 
        transformed_sigma_points_measurement=zeros(1,n_sigma_points);
        for ii = 1:n_sigma_points
            transformed_sigma_points_measurement(:,ii)=compute_h(transformed_sigma_points_model(:,ii));
        end
        
        % Compute the transformed sigma points mean
        mean_sigma_points_measurement=transformed_sigma_points_measurement*Wm';
        
        % Compute the transformed sigma points covariance
        Covariance_sigma_points_measurement=zeros(1);
        for ii=1:length(Wc)
            deviation_ii = transformed_sigma_points_measurement(:,ii)-mean_sigma_points_measurement;
            Covariance_sigma_points_measurement=Covariance_sigma_points_measurement+Wc(ii)*(deviation_ii)*(deviation_ii)';
        end
        Covariance_sigma_points_measurement = Covariance_sigma_points_measurement + r;
        
        %% Covariance between measurement and state 
        Covariance_sigma_points_both=zeros(L,1);
        for ii=1:length(Wc)
            deviation_model_ii = transformed_sigma_points_model(:,ii)-mean_sigma_points_model;
            deviation_measurement_ii = transformed_sigma_points_measurement(:,ii)-mean_sigma_points_measurement;
            Covariance_sigma_points_both=Covariance_sigma_points_both+Wc(ii)*(deviation_model_ii)*(deviation_measurement_ii)';
        end

        K = Covariance_sigma_points_both*inv(Covariance_sigma_points_measurement);
        
        e = signal(k)-mean_sigma_points_measurement;
        x_pred = mean_sigma_points_model+K*(e);

        P = Covariance_sigma_points_model-K*Covariance_sigma_points_measurement*K';
        P(P<2*eps)=0;
        
        pred_vec=[pred_vec, x_pred];
        K_vec = [K_vec, K];
        P_vec = [P_vec, P];
        e_vec = [e_vec, e];
    end
    
end

function x_pred_model = compute_f(x) % Use the model to predict next state
    x_pred_model = [    cos(x(3))*x(1)-sin(x(3))*x(2)
                        sin(x(3))*x(1)+cos(x(3))*x(2)
                        x(3)                            ]; 
end

function y_bar = compute_h(x) % Compute measurement
    y_bar = x(1);
end