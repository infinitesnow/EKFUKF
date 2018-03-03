function static_plot(signal,omega,pred_vec_ekf,pred_vec_ukf)
    t=1:length(signal);
    subplot(3,1,1)
    plot(t,omega,'k-',t,pred_vec_ekf(3,:),'r-',t,pred_vec_ukf(3,:),'b');
    title('Omega estimation');
    legend('True','Estimated EKF','Estimated UKF');
    subplot(3,1,2)
    plot(t,signal,'ko',t,pred_vec_ekf(1,:),'b-',t,pred_vec_ekf(2,:),'r-');
    title('States EKF');
    subplot(3,1,3)
    plot(t,signal,'ko',t,pred_vec_ukf(1,:),'b-',t,pred_vec_ukf(2,:),'r-');
    title('States UKF');
end