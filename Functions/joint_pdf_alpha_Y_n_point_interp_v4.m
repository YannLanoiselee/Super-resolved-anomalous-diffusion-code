function [mat] = joint_pdf_alpha_Y_n_point_interp_v4(x_Y_est,x_alpha_est,list_H_D,p,Y_alpha_stat)
%% Interpolation of the moments of Y and alpha
Y_alpha_stat.list_Y_scan=log(2*Y_alpha_stat.list_D_scan);

if numel(Y_alpha_stat.list_alpha_scan)>2
    if size(Y_alpha_stat.list_alpha_scan,1)>2
Y_mean_th= interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_Y_scan,Y_alpha_stat.Y_mean_th,2*list_H_D(:,1),log(2*list_H_D(:,2)),'cubic');
alpha_mean_th= interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_Y_scan,Y_alpha_stat.alpha_mean_th,2*list_H_D(:,1),log(2*list_H_D(:,2)),'cubic');
Y_var_th = interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_Y_scan,Y_alpha_stat.Y_var_th,2*list_H_D(:,1),log(2*list_H_D(:,2)),'cubic');
alpha_var_th = interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_Y_scan,Y_alpha_stat.alpha_var_th,2*list_H_D(:,1),log(2*list_H_D(:,2)),'cubic');
covar_Y_alpha_th = interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_Y_scan,Y_alpha_stat.covar_Y_alpha_th,2*list_H_D(:,1),log(2*list_H_D(:,2)),'cubic');
    else
Y_mean_th= interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_Y_scan,Y_alpha_stat.Y_mean_th,2*list_H_D(:,1),log(2*list_H_D(:,2)),'linear');
alpha_mean_th= interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_Y_scan,Y_alpha_stat.alpha_mean_th,2*list_H_D(:,1),log(2*list_H_D(:,2)),'linear');
Y_var_th = interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_Y_scan,Y_alpha_stat.Y_var_th,2*list_H_D(:,1),log(2*list_H_D(:,2)),'linear');
alpha_var_th = interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_Y_scan,Y_alpha_stat.alpha_var_th,2*list_H_D(:,1),log(2*list_H_D(:,2)),'linear');
covar_Y_alpha_th = interp2(Y_alpha_stat.list_alpha_scan,Y_alpha_stat.list_Y_scan,Y_alpha_stat.covar_Y_alpha_th,2*list_H_D(:,1),log(2*list_H_D(:,2)),'linear');
    end
else
Y_mean_th= Y_alpha_stat.Y_mean_th;
alpha_mean_th= Y_alpha_stat.alpha_mean_th;
Y_var_th = Y_alpha_stat.Y_var_th;
alpha_var_th = Y_alpha_stat.alpha_var_th;
covar_Y_alpha_th = Y_alpha_stat.covar_Y_alpha_th;
end

    x1=[x_Y_est(:) x_alpha_est(:)]';
mat=zeros(size(x_Y_est,1),size(x_alpha_est,2));
for n_pair=1:size(list_H_D,1)
%% theoretical mean and covariance of Y and alpha
    mn=[Y_mean_th(n_pair),alpha_mean_th(n_pair)]';
    C_th=[Y_var_th(n_pair),covar_Y_alpha_th(n_pair);...
        covar_Y_alpha_th(n_pair),alpha_var_th(n_pair)];

    %% multivariate Gaussian
%inv_C_th=1/(C_th(1)*C_th(4)-C_th(2)*C_th(3))*[C_th(4),-C_th(2);-C_th(3),C_th(1)];
 %   mulGau= 1/(2*pi*det(C_th)^(1/2))*exp(-0.5.*(x1-mn)'*inv_C_th*(x1-mn));
  %  G=reshape(diag(mulGau),size(x_Y_est,1),size(x_alpha_est,1));
      mulGau= mvnpdf(x1',mn',C_th);
    G=reshape((mulGau),size(x_Y_est,2),size(x_alpha_est,2));
   mat=mat+p(n_pair)*G;
end

end
