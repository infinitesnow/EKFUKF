set(0, 'currentfigure', ukf_figure);

%% State estimation
subplot(5,3,[1 2 3])
xlim([k-window_size k])
addpoints(x1_line,k,x_pred(1))
addpoints(x1true_line,k,signal(k))
addpoints(x2_line,k,x_pred(2))

%% Omega estimation
subplot(5,3,[4 5 6])
xlim([k-window_size k])
addpoints(x3_line,k,x_pred(3))

%% Kalman gain
subplot(5,3,[7 8 9])
xlim([k-window_size k])
addpoints(k1_line,k,K(1))
addpoints(k2_line,k,K(2))
addpoints(k3_line,k,K(3))

%% Innovation
subplot(5,3,[10 11 12])
xlim([k-window_size k])
addpoints(e_line,k,e)

%% P
subplot(5,3,[ 13 14 15 ])
imagesc(P)
colorbar

UKF(k) = getframe(ukf_figure);