function pred_vec = ukf(signal,x_pred_0,sigma_init,varargin)
    %% Settings
    global window_size 
    window_size = 50;
    PLOT=false;
    PLOT_K=false;
    PLOT_e=false;
    PLOT_P=false;    
    SHOW_SIGMA_POINTS = false;

    if (PLOT) 
        initialize_plot();
    end
    if (PLOT_K)
        initialize_plot_K();
    end
    if (PLOT_e)
        initialize_plot_e();
    end
    if (PLOT_P)
        initialize_plot_P();
    end
    if (SHOW_SIGMA_POINTS)
        initialize_sp_plot();
    end
    
    %% Parse arguments
    % only want 5 optional inputs at most
    numvarargs = length(varargin);
    if numvarargs > 5
        error('Too many inputs');
    end

    % Set defaults for optional inputs
    optargs = {...
        1e-13, ...%r
        1e-10, ...%q
        1, ...%alpha
        2, ...%beta
        2 ...%k
    };

    % Now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;

    % Place optional args in variable names
    [r,q,alpha,beta,k] = optargs{:};
            
    %% Set initial values
    L = size(compute_f(x_pred_0'),1); % The state space size;
    L_m = size(compute_h(x_pred_0'),1); % The measurement state space size
    simulation_length=length(signal);
    
    lambda=alpha^2*(L+k)-L;
    gamma=sqrt(L+lambda);
    
    I3=eye(L);
    Q=q*I3;
    
    x_pred=x_pred_0';
    P=sigma_init*eye(L);
        
    %% Compute weights
    % Weights to compute the mean of the transformed sigma points
    Wm0=lambda/(L+lambda);
    Wmi=1/(2*(L+lambda))*ones(1,2*L);
    Wm=[Wm0, Wmi];
%     Wm=Wm/sum(Wm); % Weight normalization?
    
    % Weights to compute the covariance of the transformed sigma points
    Wc0=Wm0+1-alpha^2+beta;
    Wci=Wmi;
    Wc=[Wc0, Wci]; 
%     Wc=Wc/sum(Wc); % Weight normalization?
    
    % Check that dimensions are valid (both row, columns)
    assert(all(size(Wm)==[1 2*L+1])); 
    assert(all(size(Wc)==[1 2*L+1]));
    
    pred_vec=[];
    
    %% Run prediction
    for k=1:simulation_length
        
        sqrtP=chol(P,'lower'); % The matrix square root of P
                
        %% Compute the sigma points
        Sigma_points_plus = zeros(L,L);
        Sigma_points_minus = zeros(L,L);
        for ii=1:L
            Sigma_points_plus(:,ii)  = x_pred + gamma*sqrtP(:,ii);
            Sigma_points_minus(:,ii) = x_pred - gamma*sqrtP(:,ii);
        end
        Sigma_points = [x_pred Sigma_points_plus Sigma_points_minus];
        n_sigma_points = size(Sigma_points,2);
        assert(n_sigma_points==2*L+1);
        
        %% Model update
        % Propagate the sigma points through the model
        transformed_sigma_points_model=zeros(L,n_sigma_points);
        for ii = 1:n_sigma_points
            transformed_sigma_points_model(:,ii)=compute_f(Sigma_points(:,ii));
        end
        
        % Compute the model-propagated sigma points mean
        mean_sigma_points_model=transformed_sigma_points_model*Wm';
        
        % Compute the model-propagated sigma points covariance
        Covariance_sigma_points_model=zeros(L);
        for ii=1:length(Wc)
            deviation_ii = transformed_sigma_points_model(:,ii)-mean_sigma_points_model;
            Covariance_sigma_points_model=Covariance_sigma_points_model+Wc(ii)*(deviation_ii)*(deviation_ii)';
        end
        Covariance_sigma_points_model = Covariance_sigma_points_model + Q;
             
        %% Measurement update
        
        % Propagate the sigma points through the measurement function
        % In our case, we have no measurement function, so h is the
        % identity function. 
        transformed_sigma_points_measurement=zeros(L_m,n_sigma_points);
        for ii = 1:n_sigma_points
            transformed_sigma_points_measurement(:,ii)=compute_h(transformed_sigma_points_model(:,ii));
        end
        
        % Compute the measurement-propagated sigma points mean
        mean_sigma_points_measurement=transformed_sigma_points_measurement*Wm';
                
        % Compute the measurement-propagated sigma points covariance
        Covariance_sigma_points_measurement=zeros(L_m);
        for ii=1:length(Wc)
            deviation_ii = transformed_sigma_points_measurement(:,ii)-mean_sigma_points_measurement;
            Covariance_sigma_points_measurement=Covariance_sigma_points_measurement+Wc(ii)*(deviation_ii)*(deviation_ii)';
        end
        Covariance_sigma_points_measurement = Covariance_sigma_points_measurement + r;
                
        % Covariance between measurement and state 
        Covariance_sigma_points_statemeasurement=zeros(L,L_m);
        for ii=1:length(Wc)
            deviation_model_ii = transformed_sigma_points_model(:,ii)-mean_sigma_points_model;
            deviation_measurement_ii = transformed_sigma_points_measurement(:,ii)-mean_sigma_points_measurement;
            Covariance_sigma_points_statemeasurement=Covariance_sigma_points_statemeasurement+Wc(ii)*(deviation_model_ii)*(deviation_measurement_ii)';
        end

        K = Covariance_sigma_points_statemeasurement*inv(Covariance_sigma_points_measurement);
        
        e = signal(k)-mean_sigma_points_measurement;
        x_pred = mean_sigma_points_model+K*(e);

        P = Covariance_sigma_points_model-K*Covariance_sigma_points_measurement*K';
        
        pred_vec=[pred_vec, x_pred];
        
        if (PLOT && ii~=simulation_length)
            do_plot(x_pred,signal(k+1),k)
        end
        if (SHOW_SIGMA_POINTS)
            show_sigma_points(Sigma_points)
        end
        if (PLOT_K) 
            plot_K(K,k);
        end
        if (PLOT_e)
            plot_e(e,k);
        end
        if (PLOT_P)
            plot_P(P);
        end
        if (PLOT || PLOT_K || PLOT_e || PLOT_P)
            drawnow
        end
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


function initialize_plot()
    figure(1)
    clf
    set(figure(1), 'position', [ 0 0 500 300]) 
    movegui('northeast')
    subplot(6,3,[1 2 3])
    hold on
    xlim([0 20])
    ylim auto
    title('State estimation');
    global x1_line
    global x1true_line
    global x2_line
    x1_line=animatedline('Color','r','LineStyle','-');
    x2_line=animatedline('Color','b','LineStyle','-');
    x1true_line=animatedline();
    legend('x(1)','x(2)','Signal')
    subplot(6,3,[4 5 6])
    global x3_line
    x3_line=animatedline('Color','g','LineStyle','-');
    legend('estimated omega')
end
function do_plot(x,x1true,ii)
    global window_size
    subplot(6,3,[1 2 3])
    xlim([ii-window_size ii])
    global x1_line
    global x1true_line
    global x2_line
    addpoints(x1_line,ii,x(1))
    addpoints(x1true_line,ii,x1true)
    addpoints(x2_line,ii,x(2))
    subplot(6,3,[4 5 6])
    global x3_line
    xlim([ii-window_size ii])
    addpoints(x3_line,ii,x(3))
end
function initialize_plot_K()
    subplot(6,3,[7 8 9])
    hold on
    xlim([0 20])
    ylim auto
    title('Kalman gain');
    global k1_line
    global k2_line
    global k3_line
    k1_line=animatedline('Color','r','LineStyle','-');
    k2_line=animatedline('Color','b','LineStyle','-');
    k3_line=animatedline('Color','g','LineStyle','-');
    legend('K(1)','K(2)','K(3)')
end
function plot_K(K,ii)
    global window_size
    subplot(6,3,[7 8 9])
    global k1_line
    global k2_line
    global k3_line
    xlim([ii-window_size ii])
    addpoints(k1_line,ii,K(1))
    addpoints(k2_line,ii,K(2))
    addpoints(k3_line,ii,K(3))
end
function initialize_plot_e()
    subplot(6,3,[10 11 12])
    hold on
    xlim([0 20])
    ylim auto
    title('Innovation');
    global e_line
    e_line=animatedline('Color','k','LineStyle','-');
end
function plot_e(e,ii)
    global window_size
    subplot(6,3,[10 11 12])
    xlim([ii-window_size ii])
    global e_line
    addpoints(e_line,ii,e)
end
function initialize_plot_P()
    subplot(6,3,[13 14 15])
end
function initialize_sp_plot()
    figure(2)
    clf
    set(figure(2), 'position', [ 0 0 500 250]) 
    movegui('southeast')
    title('Sigma points')
    view([0.9,1.1,0.9]);
    figure(1)
end
function plot_P(P)
    subplot(6,3,[13 14])
    imagesc(P)
    colorbar
end
function show_sigma_points(Sigma_points)
    figure(2)
    hold on
    for ii = 1:size(Sigma_points,2)
        scatter3(Sigma_points(1,ii),Sigma_points(2,ii),Sigma_points(3,ii))
    end
    figure(1)
end