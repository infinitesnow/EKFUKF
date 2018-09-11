%% Initialize plot
ekf_figure = figure('visible','off');
set(ekf_figure, 'Position', [0, 0, 720, 720],'Resize','off')

%% Internal states estimation
subplot(6,3,[1 2 3])
hold on
xlim([0 20])
ylim auto
x1_line=animatedline('Color','r','LineStyle','-');
x1true_line=animatedline();
legend('x(1)','Signal')
subplot(6,3,[4 5 6])
x2_line=animatedline('Color','b','LineStyle','-');
legend('x(2)')

%% Omega estimation
subplot(6,3,[7 8 9])
x3_line=animatedline('Color','g','LineStyle','-');
legend('estimated omega')

%% Kalman gain
subplot(6,3,[10 11 12])
hold on
xlim([0 20])
ylim auto
title('Kalman gain');
k1_line=animatedline('Color','r','LineStyle','-');
k2_line=animatedline('Color','b','LineStyle','-');
k3_line=animatedline('Color','g','LineStyle','-');
legend('K(1)','K(2)','K(3)')

%% Innovation
subplot(6,3,[13 14 15])
hold on
xlim([0 20])
ylim auto
title('Innovation');
e_line=animatedline('Color','k','LineStyle','-');

%% P
subplot(6,3,[16 17 18])