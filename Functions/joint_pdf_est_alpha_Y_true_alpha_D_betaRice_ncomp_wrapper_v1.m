function [mat_Y_alpha_est] = joint_pdf_est_alpha_Y_true_alpha_D_betaRice_ncomp_wrapper_v1(x,xData,n_bin_true,X_alpha_stat_scan)


%wrapper for fitting joint pdf of alpha and Y
%%
x_Y_est=xData(:,:,1);
x_alpha_est=xData(:,:,2);
%% joint PDF parameters of true alpha and D
alpha_mean=x(1:6:end);
alpha_std=x(2:6:end);
alpha_min=0;
alpha_max=2;
D_mean=x(3:6:end);
D_std=x(4:6:end);
rho=x(5:6:end);
p=x(6:6:end);
%%
%%
% % % [mat_Y_alpha_est] =joint_pdf_alpha_Y_est_fct_true_alpha_D_betalognorm(x_Y_est,x_alpha_est,x_D_true,x_alpha_true,a,b,alpha_min,alpha_max,D_mulognorm,D_stdlognorm,rho,X_alpha_stat_scan);
% % [mat_Y_alpha_est] =joint_pdf_alpha_Y_est_fct_true_alpha_D_betalognorm_v2(x_Y_est,x_alpha_est,n_bin_true,a,b,alpha_min,alpha_max,D_mulognorm,D_stdlognorm,rho,X_alpha_stat_scan);
% [mat_Y_alpha_est] =joint_pdf_alpha_Y_est_fct_true_alpha_D_betalognorm_v3(x_Y_est,x_alpha_est,n_bin_true,alpha_mean,alpha_std,alpha_min,alpha_max,D_mulognorm,D_stdlognorm,rho,X_alpha_stat_scan);

[mat_Y_alpha_est]=joint_pdf_alpha_Y_est_fct_true_alpha_D_betaRice_ncomp_v1(x_Y_est,x_alpha_est,n_bin_true,alpha_mean,alpha_std,alpha_min,alpha_max,D_mean,D_std,rho,p,X_alpha_stat_scan);

end