function [fit_alpha_mean_std_D_mean_std_rho_p,Y_alpha_stat_scan,mat_res,fitted_mat,x_Y,x_alpha,resnorm] = fun_fit_joint_PDF_alpha_Y_est_nocorr_continuous_beta_Rice_v1(list_alpha,list_D,parameter)
n_bin=parameter.n_bin;
n_bin_true=parameter.n_bin_true;
n_bin_scan=parameter.n_bin_scan;
n_peak=parameter.n_peak;
N=parameter.N;
lagtime_max=parameter.lagtime_max;
dt=parameter.dt;



%%
list_Y=log(2*list_D);
%%
x_Y=linspace(min(list_Y),max(list_Y),n_bin);
x_alpha=linspace(min(list_alpha),max(list_alpha),n_bin);
[x,y]=meshgrid(x_Y,x_alpha);% Y,alpha
%%
[mat_res] = histcounts2(list_Y',list_alpha',x_Y,x_alpha,'Normalization','pdf');



%% Centers of histrogam 
X=(x_Y(1:end-1)+x_Y(2:end))/2;%Y
Y=(x_alpha(1:end-1)+x_alpha(2:end))/2;%alpha

[x,y]=meshgrid(X,Y);% Y,alpha
%% parameter scan
border_alpha=[max([0,mean(list_alpha)-6*std(list_alpha)]),min([2,mean(list_alpha)+6*std(list_alpha)])];
border_D=[min(list_D)/100,max(list_D)*10];
[Y_alpha_stat_scan] = Y_alpha_statistics_scan_nocorr_v1(border_alpha,border_D,n_bin_scan,N,lagtime_max,dt);

%% GMM fitting
GMModels = cell(n_peak,1); % Preallocation
Mu=cell(n_peak,1);
Sigma=cell(n_peak,1);

options = statset('MaxIter',4*10^3,'MaxFunEvals',4*10^3,'TolFun',10^-10);
for j =n_peak% 1:n_peak_max
    GMModels{j} = fitgmdist([list_alpha',list_Y'],j,'SharedCovariance',false,'Options',options);
    Mu{j} = GMModels{j}.mu;
   Sigma{j}= [squeeze(GMModels{j}.Sigma(1,1,:)), squeeze(GMModels{j}.Sigma(2,2,:))];
  p_GMM{j}=GMModels{j}.ComponentProportion';  %Mu{j}(:,2)=Mu{j}(:,2);
end



%%
%% Initial guess for alpha mean and D mean taking in account noise correction

alpha_mean_I=Mu{n_peak}(:,1);%}mean(list_alpha);
Y_mean_I=Mu{n_peak}(:,2);
x_alpha_test=Y_alpha_stat_scan.list_alpha_scan(1,:);
x_D_test=Y_alpha_stat_scan.list_D_scan(:,1);
[x_alpha_test,x_D_test]=meshgrid(x_alpha_test,x_D_test);
n_interp=200;
x_alpha_interp=linspace(min(x_alpha_test(1,:)),max(x_alpha_test(1,:)),n_interp);
x_D_interp=logspace(log10(min(x_D_test(:,1))),log10(max(x_D_test(:,1))),n_interp);
[x_alpha_interp,x_D_interp]=meshgrid(x_alpha_interp,x_D_interp);
alpha_mean_th_interp= interp2(x_alpha_test,x_D_test,Y_alpha_stat_scan.alpha_mean_th,x_alpha_interp,x_D_interp,'spline');
X_mean_th_interp= interp2(x_alpha_test,x_D_test,Y_alpha_stat_scan.Y_mean_th,x_alpha_interp,x_D_interp,'spline');

alpha_var_th_interp= interp2(x_alpha_test,x_D_test,Y_alpha_stat_scan.alpha_var_th,x_alpha_interp,x_D_interp,'spline');
X_var_th_interp= interp2(x_alpha_test,x_D_test,Y_alpha_stat_scan.Y_var_th,x_alpha_interp,x_D_interp,'spline');


for n_comp=1:n_peak
test_alpha=abs(alpha_mean_th_interp-alpha_mean_I(n_comp));
test_D=abs(X_mean_th_interp-Y_mean_I(n_comp));

test_alpha=abs(test_alpha)./max(abs(test_alpha(:)));
test_D=abs(test_D)./max(abs(test_D(:)));

test=abs(test_alpha).^2+abs(test_D).^2;
 
[row,col]=find(test==min(test(:)));

D_mean0(n_comp,1)=x_D_interp(row,1)*exp(X_var_th_interp(row,col)/2);
alpha_mean0(n_comp,1)=x_alpha_interp(1,col);


Y_std0=X_var_th_interp(row,col)^0.5;
alpha_est_std0=alpha_var_th_interp(row,col)^0.5;

end




%% initial conditions for fitting
alpha_std0=Sigma{n_peak}(:,1).^0.5-alpha_est_std0/2;%std(list_alpha);%/10;
mu0=(alpha_mean0-0)/(2-0);
phi0=(2-0)^2*mu0.*(1-mu0)./alpha_std0.^2-1;
a0=mu0.*phi0;
b0=(1-mu0).*phi0;
%D_mean0=exp(Mu{n_peak}(:,2)+Sigma{n_peak}(:,2).^2/2)/2;%mean(list_D);
%D_std0=sqrt(4*(exp((Sigma{n_peak}(:,2).^0.5-Y_std0).^2)-1).*D_mean0.^2);%std(list_D);
D_std0=sqrt(1/2*(exp((Sigma{n_peak}(:,2).^0.5-0*Y_std0).^2)-1).*D_mean0.^2);%std(list_D);

p0=p_GMM{n_peak};


%% fitting
% a=x(1);
% b=x(2);
% alpha_min=0;
% alpha_max=2;
% D_mulognorm=x(3);
% D_stdlognorm=x(4);
% rho=x(5);

x0=[alpha_mean0,alpha_std0,D_mean0,D_std0,zeros(n_peak,1),p0];
x0=x0';
x0=x0(:)';
lb=[0,0,10^-6,10^-6,-0.99,0];
ub=[2,inf,inf,inf,0.99,1];

lb=repmat(lb,1,n_peak);
ub=repmat(ub,1,n_peak);

xdata(:,:,1)=x;
xdata(:,:,2)=y;
 fun=@(x,xdata)joint_pdf_est_alpha_Y_true_alpha_D_betaRice_ncomp_wrapper_v1(x,xdata,n_bin_true,Y_alpha_stat_scan);
options = optimoptions('lsqcurvefit','TolFun',10^-6,'TolX',10^-6,'MaxFunctionEvaluations',2*10^3);
Aeq=zeros(1,numel(x0));
Aeq(6:6:end)=1;
[fit_param,resnorm] = lsqcurvefit(fun,x0,xdata,mat_res',lb,ub,options);%,[], [],Aeq,1);
fit_param(6:6:end)=fit_param(6:6:end)/sum(fit_param(6:6:end));
mat_fit=fun(fit_param,xdata);





alpha_mean_fit=fit_param(1:6:end);
alpha_std_fit=fit_param(2:6:end);
alpha_min=0;
alpha_max=2;
D_mean_fit=fit_param(3:6:end);
D_std_fit=fit_param(4:6:end);
rho_fit=fit_param(5:6:end);
p_fit=fit_param(6:6:end);


%%
fit_alpha_mean_std_D_mean_std_rho_p=[alpha_mean_fit',alpha_std_fit',D_mean_fit',D_std_fit',rho_fit',p_fit'];
fitted_mat=fun(fit_param,xdata);
end