%% Initialize plot
sp_figure = figure('visible','off');
set(sp_figure, 'Position', [0, 0, 720, 720],'Resize','off')

title('Sigma points')
view([0.9,1.1,0.9]);
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
zlim([-0.5 0.5]);
hold on