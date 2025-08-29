clear;
close all;
%% load simulated trajectories with three components
load('simulated_trajectory_example_continuous_alpha_D.mat','Trajectory','table_ground_truth');
% 'list_true_anomalous_diffusion_parameters' is a structure containing the input parameters of
% the simulated data
% Example dataset contains M=400 trajectories of N=100 points
dt=0.03; % time interval between two position recordings (time resolution)

%%

parameter_analysis.lagtime_max=5; % maximum lag time (tau_max in the article) for fitting TAMSD
parameter_analysis.n_peak=2; % number of component of the true parameters continuous distribution p(alpha,D)
parameter_analysis.n_bin=20;% number of bins of the histograms of alpha and y for fitting
parameter_analysis.n_bin_scan=50;% number of bins for precomputing the moments of estimated parameters as a function of the true ones before interpolation.
parameter_analysis.n_bin_true=20;% number of bins to discretize the theoretical true distribution 
parameter_analysis.pdf_marginal_D='beta_rice';%'beta_lognormal';%
% joint PDF of the true anoomalous diffusion parameters. 
% Two options: 'beta_lognormal' and 'beta_rice'. 
% For the anomalous exponent alpha : 
% beta distribution with support on the interval [0,2] in any case  
% For D:
% Lognormal distribution is approximately Gaussian when its mean is larger than its
% standard deviation (std), it not very stable when std is larger than mean.
% Rice distribution is approximately Gaussian when the mean-3std>0, it is
% stable over the whole range of parameters.
% Both distribution have strictly positive support.
%  The behavior near zeros is different. Rice distribution is advised first.


%% function for fitting and plotting results
[fit_alpha_mean_std_D_mean_std_rho_p,list_D,list_alpha,resnorm] = fit_joint_PDF_alpha_D_continuous_distribution_v2(Trajectory,dt,parameter_analysis);
%% fitting results
alpha_mean_fit=fit_alpha_mean_std_D_mean_std_rho_p(:,1);
alpha_std_fit=fit_alpha_mean_std_D_mean_std_rho_p(:,2);

if strcmp(parameter_analysis.pdf_marginal_D,'beta_rice')
    L_half=@(x)exp(x/2+log((1-x).*besseli(0,-x/2)-x.*besseli(1,-x/2) ));
    nu=fit_alpha_mean_std_D_mean_std_rho_p(:,3);
    sigma=fit_alpha_mean_std_D_mean_std_rho_p(:,4);
D_mean_fit=sigma.*sqrt(pi/2).*L_half(-nu.^2./(2.*sigma.^2));
D_std_fit=(2*sigma.^2+nu.^2-pi*sigma.^2/2.*L_half(-nu.^2./(2*sigma.^2)).^2).^0.5;

else
D_mean_fit=fit_alpha_mean_std_D_mean_std_rho_p(:,3);
D_std_fit=fit_alpha_mean_std_D_mean_std_rho_p(:,4);

end
rho_fit=fit_alpha_mean_std_D_mean_std_rho_p(:,5);
p_fit=fit_alpha_mean_std_D_mean_std_rho_p(:,6);

table_result=array2table(fit_alpha_mean_std_D_mean_std_rho_p,...
    'VariableNames',{'alpha mean','alpha std','D mean','D std','correlation','proportion'})

