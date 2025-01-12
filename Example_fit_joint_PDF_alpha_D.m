
clear;
 close all;
%% load simulated trajectories with three components
load('simulated_trajectory_example.mat','Trajectory','list_alpha_true','list_D_true','list_p_true','dt')



%%
lagtime_max=5; % maximum lag time (tau_max in the article) for foitting TAMSD
n_peak=3; % number of component of the true parameters p(alpha,D)
n_bin=20;% number of bins of the histograms of alpha and y for fitting
n_bin_scan=50;% number of bins for precomputing the moments of estimated parameters as a function of the true ones before interpolation.


%% function for fitting
[fit_alpha_D_p,list_D,list_alpha,resnorm] = fit_joint_PDF_alpha_D_from_trajectory(Trajectory,dt,lagtime_max,n_peak,n_bin,n_bin_scan);

disp(['squared residual of the fitting ',num2str(resnorm)])


%% comparison to groundtruth
ground_truth=[list_alpha_true;list_D_true;list_p_true];
IDX=zeros(1,n_peak);
for n=1:n_peak
I=sum((ground_truth-repmat(fit_alpha_D_p(:,n),1,n_peak)).^2,1);
[~,IDX(n)]=min(I);
end
fit_alpha_D_p=fit_alpha_D_p(:,IDX);

results_interlaced= reshape ([ ground_truth' ;fit_alpha_D_p'], size(ground_truth,1), [] );

table_result=array2table(results_interlaced,...
    'VariableNames',{'alpha true','alpha fitted','D true','D fitted','proportion true','propotion fitted'})
