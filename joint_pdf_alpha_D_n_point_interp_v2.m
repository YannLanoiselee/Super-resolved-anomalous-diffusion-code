function [mat] = joint_pdf_alpha_D_n_point_interp_v2(x_D_est,x_alpha_est,list_H_D,p,Y_alpha_stat,fit_factor)
%% new function that computes all at once (faster)
% [Y_alpha_stat] = covar_Y_alpha_fBm_noise_mean_var_v3(2*list_H_D(:,1),list_H_D(:,2),N,lagtime_max,sig_error);

% statistics are precalculated to improve speed (much faster)

if numel(Y_alpha_stat.list_alpha_scan)>2
    if size(Y_alpha_stat.list_alpha_scan,1)>2
Y_mean_th= interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_D_scan,Y_alpha_stat.Y_mean_th,2*list_H_D(:,1),list_H_D(:,2),'cubic');
alpha_mean_th= interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_D_scan,Y_alpha_stat.alpha_mean_th,2*list_H_D(:,1),list_H_D(:,2),'cubic');
Y_var_th = interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_D_scan,Y_alpha_stat.Y_var_th,2*list_H_D(:,1),list_H_D(:,2),'cubic');
alpha_var_th = interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_D_scan,Y_alpha_stat.alpha_var_th,2*list_H_D(:,1),list_H_D(:,2),'cubic');
covar_Y_alpha_th = interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_D_scan,Y_alpha_stat.covar_Y_alpha_th,2*list_H_D(:,1),list_H_D(:,2),'cubic');
    else
Y_mean_th= interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_D_scan,Y_alpha_stat.Y_mean_th,2*list_H_D(:,1),list_H_D(:,2),'linear');
alpha_mean_th= interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_D_scan,Y_alpha_stat.alpha_mean_th,2*list_H_D(:,1),list_H_D(:,2),'linear');
Y_var_th = interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_D_scan,Y_alpha_stat.Y_var_th,2*list_H_D(:,1),list_H_D(:,2),'linear');
alpha_var_th = interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_D_scan,Y_alpha_stat.alpha_var_th,2*list_H_D(:,1),list_H_D(:,2),'linear');
covar_Y_alpha_th = interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_D_scan,Y_alpha_stat.covar_Y_alpha_th,2*list_H_D(:,1),list_H_D(:,2),'linear');

    end
else
Y_mean_th= Y_alpha_stat.Y_mean_th;
alpha_mean_th= Y_alpha_stat.alpha_mean_th;
Y_var_th = Y_alpha_stat.Y_var_th;
alpha_var_th = Y_alpha_stat.alpha_var_th;
covar_Y_alpha_th = Y_alpha_stat.covar_Y_alpha_th;

end

    x1=[log(2*x_D_est(:)) x_alpha_est(:)]';
mat=zeros(size(x_D_est,1),size(x_alpha_est,2));
for n_pair=1:size(list_H_D,1)
    mn=[Y_mean_th(n_pair),alpha_mean_th(n_pair)]';
    %  C_th=[Y_alpha_stat.Y_var_th(n_pair),Y_alpha_stat.covar_Y_alpha_th(n_pair);...
    %        Y_alpha_stat.covar_Y_alpha_th(n_pair),Y_alpha_stat.alpha_var_th(n_pair)];
    C_th=[Y_var_th(n_pair),covar_Y_alpha_th(n_pair);...
        covar_Y_alpha_th(n_pair),alpha_var_th(n_pair)];

    %%
    % [x,y]=meshgrid(X',Y');

    %%%multivar Gassiaan
    mulGau= 1/(2*pi*det(C_th)^(1/2))*exp(-0.5.*(x1-mn)'*inv(C_th)*(x1-mn));
    G=1./x_D_est.*reshape(diag(mulGau),size(x_D_est,1),size(x_alpha_est,1));
    %% copula
% [pdf] = copula_beta_Gauss_D_alpha(Y,X,a,b,beta_min,beta_max,mu,sigma,rho);


    mat=mat+p(n_pair)*G;
end
mat=mat*fit_factor;
end