function pred_vec = ukf(signal,x_pred_0,initial_sigma,varargin)
%% Settings
window_size = 50;
VERBOSE = false;
SAVE_UKF_PLOT = false;   
SAVE_SP_PLOT = false;

if (SAVE_UKF_PLOT) 
    initialize_plot_ukf
end
if (SAVE_SP_PLOT)
    initialize_sp_plot
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
P=initial_sigma*eye(L);

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

    if (VERBOSE)
        fprintf('Iteration %d...\n',k);
    end
    if (SAVE_UKF_PLOT && ii~=simulation_length)
        plot_ukf
    end
    if (SAVE_SP_PLOT && ii~=simulation_length)
        plot_sp
    end
end

if (SAVE_UKF_PLOT | SAVE_SP_PLOT)
    addpath('./generate/')
    if (SAVE_UKF_PLOT)
        generate_video('ukf',UKF);
    end
    if (SAVE_SP_PLOT)
        SP(1) = [];
        generate_video('sp',SP);
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