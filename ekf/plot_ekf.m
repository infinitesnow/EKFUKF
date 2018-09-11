%% Internal states estimation
subplot(6,3,[1 2 3])
xlim([ii-window_size ii])
addpoints(x1_line,ii,x_pred(1))
addpoints(x1true_line,ii,signal(ii))

%% Omega estimation
subplot(6,3,[4 5 6])
xlim([ii-window_size ii])
addpoints(x2_line,ii,x_pred(2))

%% Kalman gain
subplot(6,3,[7 8 9])
xlim([ii-window_size ii])
addpoints(x3_line,ii,x_pred(3))

%% Innovation
subplot(6,3,[10 11 12])
xlim([ii-window_size ii])
addpoints(k1_line,ii,K(1))
addpoints(k2_line,ii,K(2))
addpoints(k3_line,ii,K(3))
subplot(6,3,[13 14 15])
xlim([ii-window_size ii])
addpoints(e_line,ii,e)

%% P
subplot(6,3,[16 17 18])
imagesc(P)
colorbar

EKF(ii) = getframe(ekf_figure);