function dynamic_plot_complete(signal,omega,pred_vec)
    %% Initialize dynamic plot
    figure(1);
    window_size=50;
    
    subplot(4,2,[1 2])
    hold on
    xlim([0 20])
    ylim auto
    trueomega_line=animatedline('Color','r');
    omega_line=animatedline;
    title('Omega estimation')
    legend('True','Estimated')
    
    subplot(4,2,[3 4])
    hold on
    xlim([0 20])
    ylim auto
    title('States')
    truesignal_line=animatedline('Color','k','LineStyle','-','Marker','o');
    x1_line=animatedline('Color','b','LineStyle','--');
    x2_line=animatedline('Color','r','LineStyle','--');
    legend('Signal','x1','x2');
    
    %% Plot
    for ii = 1:length(signal)
        x_pred=pred_vec(:,ii);
        
        subplot(2,2,[1 2])
        xlim([ii-window_size ii])
        addpoints(trueomega_line,ii,omega(ii))
        addpoints(omega_line,ii,x_pred(3))
        
        subplot(2,2,[3 4])
        xlim([ii-window_size ii])
        addpoints(truesignal_line,ii,signal(ii))
        addpoints(x1_line,ii,x_pred(1))
        addpoints(x2_line,ii,x_pred(2))
        
        drawnow
        
    end
end