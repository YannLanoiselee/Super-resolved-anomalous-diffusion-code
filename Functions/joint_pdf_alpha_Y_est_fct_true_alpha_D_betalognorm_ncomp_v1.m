function [mat_Y_alpha_est,x_alpha_true_edges,x_D_true_edges] =joint_pdf_alpha_Y_est_fct_true_alpha_D_betalognorm_ncomp_v1(x_Y_est,x_alpha_est,n_bin_true,alpha_mean,alpha_std,alpha_min,alpha_max,D_mulognorm,D_stdlognorm,rho,p,X_alpha_stat_scan)

% x_alpha_true and x D true are precalculated
[x_alpha_true_edges,x_D_true_edges] = find_bounds_true_alpha_D_beta_logn_v2(alpha_mean,alpha_std,alpha_min,alpha_max,D_mulognorm,D_stdlognorm,n_bin_true,n_bin_true);

%[x_alpha_true_edges,x_D_true_edges] = find_bounds_true_alpha_D_beta_logn_Lebesgue_v1(alpha_mean,alpha_std,alpha_min,alpha_max,D_mulognorm,D_stdlognorm,p,n_bin_true,n_bin_true);

%% joint pdf of estimated Y and alpha given that true alpha and D have a beta-lognormal distribution

x_D_true_centers=(x_D_true_edges(2:end)-x_D_true_edges(1:end-1))./(log(x_D_true_edges(2:end))-log(x_D_true_edges(1:end-1)));%D
%x_D_true_centers=(x_D_true_edges(1:end-1)+x_D_true_edges(2:end))/2;%D

x_alpha_centers=(x_alpha_true_edges(1:end-1)+x_alpha_true_edges(2:end))/2;%alpha




%%
[x_true,y_true]=meshgrid(x_D_true_centers,x_alpha_centers);

[dx_true,dy_true]=meshgrid(x_D_true_edges,x_alpha_true_edges);
dx_true=diff(dx_true,1,2);
dy_true=diff(dy_true,1,1);
dx_true=dx_true(1:end-1,1:end);
dy_true=dy_true(1:end,1:end-1);


%%
% [pdf_true_alpha_D] = copula_beta_Gauss_lognD_alpha_true(x_alpha_centers,x_D_true_centers,a,b,alpha_min,alpha_max,D_mulognorm,D_stdlognorm,rho);
[pdf_true_alpha_D] = copula_beta_Gauss_lognD_alpha_true_mean_var_ncomp(x_alpha_centers,x_D_true_centers,alpha_mean,alpha_std,alpha_min,alpha_max,D_mulognorm,D_stdlognorm,rho,p);


%%

list_H_D=[y_true(:)/2,x_true(:)];
p=pdf_true_alpha_D'.*dx_true.*dy_true;

p=p(:);
%sum(p)
I=p<10^-4;
list_H_D(I,:)=[];
p(I,:)=[];
%sum(p)

[mat_Y_alpha_est] = joint_pdf_alpha_Y_n_point_interp_v4(x_Y_est,x_alpha_est,list_H_D,p,X_alpha_stat_scan);
end