function dynamic_plot_complete(signal,omega,pred_vec,K_vec,e_vec,P_vec)
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
    
    subplot(4,2,[5 6])
    hold on
    xlim([0 20])
    ylim auto
    title('Kalman gain');
    k1_line=animatedline('Color','r','LineStyle','-');
    k2_line=animatedline('Color','b','LineStyle','-');
    k3_line=animatedline('Color','g','LineStyle','-');
    legend('K(1)','K(2)','K(3)')
    
    subplot(4,2,7)
    hold on
    xlim([0 20])
    ylim auto
    title('Innovation');
    e_line=animatedline('Color','k','LineStyle','-');
    
    %% Plot
    for ii = 1:length(signal)
        x_pred=pred_vec(:,ii);
        K=K_vec(:,ii);
        e=e_vec(ii);
        P=P_vec(:,(3*ii-2):(3*ii));
        
        subplot(4,2,[1 2])
        xlim([ii-window_size ii])
        addpoints(trueomega_line,ii,omega(ii))
        addpoints(omega_line,ii,x_pred(3))
        
        subplot(4,2,[3 4])
        xlim([ii-window_size ii])
        addpoints(truesignal_line,ii,signal(ii))
        addpoints(x1_line,ii,x_pred(1))
        addpoints(x2_line,ii,x_pred(2))
        
        subplot(4,2,[5 6])
        xlim([ii-window_size ii])
        addpoints(k1_line,ii,K(1))
        addpoints(k2_line,ii,K(2))
        addpoints(k3_line,ii,K(3))
        
        subplot(4,2,7)
        xlim([ii-window_size ii])
        addpoints(e_line,ii,e)
        
        subplot(4,2,8)
        imagesc(P)
        colorbar
        title('P');
        
        drawnow
        
    end
end