function [Y_alpha_stat_scan] = Y_alpha_statistics_scan_nocorr_v1(list_alpha,list_D,n_bin,N,lagtime_max,dt)
list_alpha_scan=linspace(max([0.001,min(list_alpha)]),min([1.95,max(list_alpha)]),n_bin);
%  list_D_scan=linspace(max([0,min(list_D)]),max(list_D),n_bin);
list_D_scan=logspace(log10(max([0,min(list_D)])),log10(max(list_D)),n_bin);

[x,y]=meshgrid(list_alpha_scan,list_D_scan);

% [Y_alpha_stat] = covar_Y_alpha_fBm_noise_mean_var_v4(x(:)',y(:)',N,lagtime_max,sig_error,dt);
[Y_alpha_stat] = covar_Y_alpha_fBm_noise_mean_var_nocorr_v1(x(:)',y(:)',N,lagtime_max,dt);
Y_alpha_stat_scan=struct;
Y_alpha_stat_scan.alpha_mean_th=reshape(Y_alpha_stat.alpha_mean_th',n_bin,n_bin);
Y_alpha_stat_scan.alpha_var_th=reshape(Y_alpha_stat.alpha_var_th',n_bin,n_bin);
Y_alpha_stat_scan.covar_Y_alpha_th=reshape(Y_alpha_stat.covar_Y_alpha_th',n_bin,n_bin);
Y_alpha_stat_scan.Y_mean_th=reshape(Y_alpha_stat.Y_mean_th',n_bin,n_bin);
Y_alpha_stat_scan.Y_var_th=reshape(Y_alpha_stat.Y_var_th',n_bin,n_bin);
Y_alpha_stat_scan.list_alpha_scan=x;
Y_alpha_stat_scan.list_D_scan=y;

end