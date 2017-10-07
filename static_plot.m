function static_plot(signal,omega,pred_vec,K_vec,e_vec)
    t=1:length(signal);
    subplot(4,1,1)
    plot(t,omega,'r-',t,pred_vec(3,:),'k-');
    title('Omega estimation');
    legend('True','Estimated');
    subplot(4,1,2)
    plot(t,signal,'ko',t,pred_vec(1,:),'b-',t,pred_vec(2,:),'r-');
    title('States');
    subplot(4,1,3)
    plot(t,K_vec(1,:),'r-',t,K_vec(2,:),'g-',t,K_vec(3,:),'b-');
    title('Kalman gain');
    legend('K(1)','K(2)','K(3)');
    subplot(4,1,4)
    plot(t,e_vec,'k-');
    title('Innovation');
end