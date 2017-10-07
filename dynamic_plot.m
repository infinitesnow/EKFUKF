function dynamic_plot(omega,pred_vec,pausetime)
    %% Initialize dynamic plot
    figure(1);
    window_size=50;
    
    hold on
    xlim([0 20])
    ylim auto
    trueomega_line=animatedline('Color','r');
    omega_line=animatedline;
    title('Omega estimation')
    legend('True','Estimated')
    
    %% Plot
    for ii = 1:length(omega)
        pred=pred_vec(3,ii);
        
        xlim([ii-window_size ii])
        addpoints(trueomega_line,ii,omega(ii))
        addpoints(omega_line,ii,pred)
        
        pause(pausetime)
        drawnow
    end
end