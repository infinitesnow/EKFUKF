%% Initialize plot
ukf_figure = figure('visible','off');
set(ukf_figure, 'Position', [0, 0, 720, 720],'Resize','off')

%% Internal states estimation
subplot(5,3,[1 2 3])
hold on
xlim([0 20])
ylim auto
title('State estimation');
x1_line=animatedline('Color','r','LineStyle','-');
x2_line=animatedline('Color','b','LineStyle','-');
x1true_line=animatedline();
legend('x(1)','x(2)','Signal')

%% Omega estimation
subplot(5,3,[4 5 6])
x3_line=animatedline('Color','g','LineStyle','-');
legend('estimated omega')

%% Kalman gain
subplot(5,3,[7 8 9])
hold on
xlim([0 20])
ylim auto
title('Kalman gain');
k1_line=animatedline('Color','r','LineStyle','-');
k2_line=animatedline('Color','b','LineStyle','-');
k3_line=animatedline('Color','g','LineStyle','-');
legend('K(1)','K(2)','K(3)')

%% Innovation
subplot(5,3,[10 11 12])
hold on
xlim([0 20])
ylim auto
title('Innovation');
e_line=animatedline('Color','k','LineStyle','-');

%% P
subplot(5,3,[13 14 15])