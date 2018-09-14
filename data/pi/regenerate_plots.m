files = dir;
directoryNames = {files([files.isdir]).name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

for directory = directoryNames
directory = cell2mat(directory);
cd(directory);
    
load('pi_curve_ekf.mat','pi_curve_ekf');
load('pi_curve_ukf.mat','pi_curve_ukf');

convergence_threshold = 10;

%% Remove iterations which didn't converge
not_converged_ekf_indices = union( ...
    find(pi_curve_ekf(1,:)>convergence_threshold), find(pi_curve_ekf(2,:)>convergence_threshold) );
converged_ekf_indices = setdiff(1:size(pi_curve_ekf,2),not_converged_ekf_indices);
fprintf('EKF didn''t converge %d times\n',length(not_converged_ekf_indices))
not_converged_ukf_indices = union( ...
    find(pi_curve_ukf(1,:)>convergence_threshold), find(pi_curve_ukf(2,:)>convergence_threshold) );
converged_ukf_indices = setdiff(1:size(pi_curve_ekf,2),not_converged_ukf_indices);
fprintf('UKF didn''t converge %d times\n',length(not_converged_ukf_indices))

converged_ekf_bool = zeros(1,size(pi_curve_ekf,2));
converged_ukf_bool = zeros(1,size(pi_curve_ukf,2));
converged_ekf_bool(converged_ekf_indices) = 1;
converged_ukf_bool(converged_ukf_indices) = 1;

pi_curve_ekf_converged = pi_curve_ekf(:,converged_ekf_indices);
pi_curve_ukf_converged = pi_curve_ukf(:,converged_ukf_indices);

pi_curves_figure = figure('visible','off');
title('PI curves')
pi_curve_ekf_et = pi_curve_ekf_converged(1,:);
pi_curve_ekf_ess = pi_curve_ekf_converged(2,:);
pi_curve_ekf_sigma = pi_curve_ekf_converged(3,:);
pi_curve_ukf_et = pi_curve_ukf_converged(1,:);
pi_curve_ukf_ess = pi_curve_ukf_converged(2,:);
pi_curve_ukf_sigma = pi_curve_ukf_converged(3,:);
hold on
scatter(pi_curve_ekf_ess,pi_curve_ekf_et,'ro')
scatter(pi_curve_ukf_ess,pi_curve_ukf_et,'bx')
legend('EKF','UKF')
xlabel('MSE steady state')
ylabel('MSE transient')
saveas(pi_curves_figure, strcat('pi_cloud_',directory,'.png'), 'png');

cd ..
end