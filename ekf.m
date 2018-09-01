function pred_vec=ekf(signal,x_pred_0,sigma_init,varargin)
    global window_size 
    window_size = 50;
    PLOT=false;
    PLOT_K=false;
    PLOT_e=false;
    PLOT_P=false;
    
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
    P=sigma_init*eye(3);
    
    H=[1 0 0];
    I3=eye(3);
    
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
        
        if (PLOT && ii~=simulation_length)
            do_plot(x_pred,signal(ii+1),ii)
        end
        if (PLOT_K) 
            plot_K(K,ii);
        end
        if (PLOT_e)
            plot_e(e,ii);
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
function F = compute_F(x) % Compute F matrix in a given state 
    F=[     cos(x(3))   -sin(x(3))  -sin(x(3))*x(1)-cos(x(3))*x(2)
            sin(x(3))   cos(x(3))    cos(x(3))*x(1)-sin(x(3))*x(2)
            0           0            1                               ];
end
function initialize_plot()
    subplot(6,3,[1 2 3])
    hold on
    xlim([0 20])
    ylim auto
    global x1_line
    global x1true_line
    x1_line=animatedline('Color','r','LineStyle','-');
    x1true_line=animatedline();
    legend('x(1)','Signal')
    subplot(6,3,[4 5 6])
    global x2_line
    x2_line=animatedline('Color','b','LineStyle','-');
    legend('x(2)')
    subplot(6,3,[7 8 9])
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
    addpoints(x1_line,ii,x(1))
    addpoints(x1true_line,ii,x1true)
    subplot(6,3,[4 5 6])
    xlim([ii-window_size ii])
    global x2_line
    addpoints(x2_line,ii,x(2))
    subplot(6,3,[7 8 9])
    xlim([ii-window_size ii])
    global x3_line
    addpoints(x3_line,ii,x(3))
end
function initialize_plot_K()
    subplot(6,3,[10 11 12])
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
    subplot(6,3,[10 11 12])
    global k1_line
    global k2_line
    global k3_line
    xlim([ii-window_size ii])
    addpoints(k1_line,ii,K(1))
    addpoints(k2_line,ii,K(2))
    addpoints(k3_line,ii,K(3))
end
function initialize_plot_e()
    subplot(6,3,[13 14 15])
    hold on
    xlim([0 20])
    ylim auto
    title('Innovation');
    global e_line
    e_line=animatedline('Color','k','LineStyle','-');
end
function plot_e(e,ii)
    global window_size
    subplot(6,3,[13 14 15])
    xlim([ii-window_size ii])
    global e_line
    addpoints(e_line,ii,e)
end
function initialize_plot_P()
    subplot(6,3,[16 17 18])
end
function plot_P(P)
    subplot(6,3,[16 17 18])
    imagesc(P)
    colorbar
end