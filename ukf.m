function pred_vec = ukf(alpha,beta,k,q,r,x_pred_0,sigma_init,signal)
    %% Set initial values
    L=3;
    simulation_length=length(signal);
    
    lambda=alpha^2*(L+k);
    gamma=L+lambda;
    
    I3=eye(3);
    Q=q*I3;
    
    x_pred=x_pred_0';
    P=sigma_init*eye(L);
    
    pred_vec=[]
    
    %% Compute weights
    Wm0=lambda/(L+lambda)+1-alpha^2+beta;
    Wmi=1/(2*(L+lambda))*ones(1,2*L);
    Wc0=lambda/(L+lambda);
    Wci=Wmi;
    Wm=[Wm0, Wmi];
    Wc=[Wc0, Wci];
    
    for k=1:simulation_length

        %% Compute propagation of sigma points
        sP=chol(gamma*P,'lower'); % The matrix square root of P
        X=zeros(L,2*L+1);
        X(:,1)=compute_f(x_pred);
        for i=1:L
            X(:,i+1)=compute_f(x_pred+sP(:,i));
        end
        for i=1:L
            X(:,i+L+1)=compute_f(x_pred-sP(:,i));
        end

        x_pred_bar=X*Wm';

        %% Compute P
        Px=zeros(L);
        for i=1:length(Wc)
            Px=Px+Wc(i)*(X(:,i)-x_pred_bar)*(X(:,i)-x_pred_bar)';
        end
        Px=Px+Q;

        Y=X(1,:);

        y_pred_bar=Y*Wm';

        Py=0;
        for i=1:length(Wc)
            Py=Py+Wc(i)*(Y(i)-y_pred_bar)*(Y(i)-y_pred_bar)';
        end
        Py=Py+r;

        Pxy=zeros(L,1);
        for i=1:length(Wc)
            Pxy=Pxy+Wc(i)*(X(:,i)-x_pred_bar)*(Y(i)-y_pred_bar)';
        end

        K=Pxy*inv(Py);
        
        x_pred=x_pred_bar+K*(signal(k)-y_pred_bar);

        P=Px-K*Py*K';
        
        pred_vec=[pred_vec, x_pred_bar];
    end
    
end

function x_pred_model = compute_f(x) % Use the model to predict next state
    x_pred_model = [    cos(x(3))*x(1)-sin(x(3))*x(2)
                        sin(x(3))*x(1)+cos(x(3))*x(2)
                        x(3)                            ]; 
end