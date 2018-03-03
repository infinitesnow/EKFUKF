function static_plot(signal,omega,pred_vec)
    t=1:length(signal);
    subplot(2,1,1)
    plot(t,omega,'r-',t,pred_vec(3,:),'k-');
    title('Omega estimation');
    legend('True','Estimated');
    subplot(2,1,2)
    plot(t,signal,'ko',t,pred_vec(1,:),'b-',t,pred_vec(2,:),'r-');
    title('States');
end