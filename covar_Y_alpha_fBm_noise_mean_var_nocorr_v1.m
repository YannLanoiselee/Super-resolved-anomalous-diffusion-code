function [Y_alpha_stat] = covar_Y_alpha_fBm_noise_mean_var_nocorr_v1(list_alpha,list_D,N,lagtime_max,dt)

covar_Y_alpha_th=zeros(numel(list_alpha),1);
alpha_var_th=zeros(numel(list_alpha),1);
alpha_mean_th=zeros(numel(list_alpha),1);
Y_mean_th=zeros(numel(list_alpha),1);
Y_var_th=zeros(numel(list_alpha),1);
 covar_D_alpha_th=zeros(numel(list_alpha),1);
 D_mean_th=zeros(numel(list_alpha),1);
 D_var_th=zeros(numel(list_alpha),1);

%% init stuff for all values of D andalpha
R_tau=cell(1,lagtime_max);
for tau=1:lagtime_max
    [R_tau{1,tau}] = R_matrix(N,tau);
end

%% loop on D and alpha
for cp_pair=1:numel(list_alpha)
    H=list_alpha(cp_pair)/2;
    D=list_D(cp_pair);

      [cov_mat] =cov_mat_fBm_dt(N,H,D,dt);


%% single lag

    mat=cell(1,lagtime_max);
    for tau=1:lagtime_max
        mat{tau}=R_tau{1,tau}*cov_mat;
    end

    %% Mean TAMSD
    M_n_mean=zeros(1,lagtime_max);
    for tau=1:lagtime_max
        M_n_mean(1,tau)=trace(mat{tau});
    end

    %% derivatives mean TAMSd over lag-times
    % mean lag time
    mean_log_tau=mean(log((1:lagtime_max)*dt));
    h_mu=sum(log(M_n_mean))/lagtime_max;
    dh=1./(lagtime_max*M_n_mean)';% first derivative
    d2h=-1./(lagtime_max*M_n_mean(1:end).^2)';% second derivative

    % anomalous exponent
    Coeff=log((1:lagtime_max))/sum(log(1:lagtime_max).^2);
    f_mu=sum(Coeff.*log(M_n_mean./M_n_mean(1)));

    df=zeros(lagtime_max,1);% first derivative
    df(1)=-sum(Coeff)/M_n_mean(1);
    df(2:end)=Coeff(2:end)./M_n_mean(2:end);

    d2f=zeros(lagtime_max,1);% second derivative
    d2f(1)=sum(Coeff)/M_n_mean(1)^2;
    d2f(2:end)=-Coeff(2:end)./M_n_mean(2:end).^2;

    % Gaussian X
    g_mu=h_mu-f_mu*mean_log_tau;

    dg=zeros(lagtime_max,1);% first derivative
    dg(1)=dh(1)-df(1)*mean_log_tau;
    dg(2:end)=dh(2:end)-df(2:end)*mean_log_tau;

    d2g=zeros(lagtime_max,1);% second derivative
    d2g(1)=d2h(1)-d2f(1)*mean_log_tau;
    d2g(2:end)=d2h(2:end)-d2f(2:end)*mean_log_tau;

    % D
gD_mu=1/2*exp(g_mu);
% dgD=zeros(lagtime_max,1);% first derivative
dgD=gD_mu.*dg;
d2gD=gD_mu*dg*dg';





    %% Covariance - variances

    [xixk_star]=covar_TAMSD_centered(lagtime_max,R_tau,cov_mat);
    %  [xixk_star]=covar_TAMSD_centered(lagtime_max,mat,mat_error);
    
    alpha_mean_th(cp_pair)=f_mu+1/2*trace(diag(d2f)'*xixk_star);
    alpha_var_th(cp_pair)=df'*xixk_star*df;% variance of alpha

    covar_Y_alpha_th(cp_pair)=dh'*xixk_star*df...
        -alpha_var_th(cp_pair)*mean_log_tau; % covariance of alpha and X
    
    Y_mean_th(cp_pair)=g_mu+1/2*trace(diag(d2g)'*xixk_star);
    Y_var_th(cp_pair)=dg'*xixk_star*dg;% variance of X
    
    covar_D_alpha_th(cp_pair)=dgD'*xixk_star*dg;% variance of X
    
    D_mean_th(cp_pair)=gD_mu+1/2*trace(d2gD'*xixk_star);% variance of X
    D_var_th(cp_pair)=dgD'*xixk_star*dgD;% variance of X

end

Y_alpha_stat=struct;
Y_alpha_stat.covar_Y_alpha_th=covar_Y_alpha_th;
Y_alpha_stat.alpha_mean_th=alpha_mean_th;
Y_alpha_stat.alpha_var_th=alpha_var_th;
Y_alpha_stat.Y_mean_th=Y_mean_th;
Y_alpha_stat.Y_var_th=Y_var_th;

Y_alpha_stat.covar_D_alpha_th=covar_D_alpha_th;
Y_alpha_stat.D_var_th=D_var_th;
Y_alpha_stat.D_mean_th=D_mean_th;


%% Function
    function [xixk_star]=covar_TAMSD_centered(lagtime_max,R_tau,cov_mat)
        xixk_star=zeros(lagtime_max ,lagtime_max );
        cov_mat2=cov_mat;
        for j=1:lagtime_max
            for k=1:lagtime_max
                xixk_star(j,k)=2*trace((R_tau{1,j}*cov_mat2)*(R_tau{1,k}*cov_mat2));%
            end
        end
    end

end
