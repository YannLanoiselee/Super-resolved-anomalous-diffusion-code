function [fit_H_D_p,Y_alpha_stat_scan,mat_res,fitted_mat,x_Y,x_alpha,resnorm] = fun_fit_joint_PDF_alpha_Y_est_nocorr_v3(list_alpha,list_D,parameter)
n_bin=parameter.n_bin;
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
%%
[mat_res] = histcounts2(list_Y',list_alpha',x_Y,x_alpha,'Normalization','pdf');
%% Centers of histrogam 
X=(x_Y(1:end-1)+x_Y(2:end))/2;%Y
Y=(x_alpha(1:end-1)+x_alpha(2:end))/2;%alpha

%% parameter scan
border_alpha=[max([0,mean(list_alpha)-6*std(list_alpha)]),min([2,mean(list_alpha)+6*std(list_alpha)])];
border_D=[min(list_D)/100,max(list_D)*10];
[Y_alpha_stat_scan] = Y_alpha_statistics_scan_nocorr_v1(border_alpha,border_D,n_bin_scan,N,lagtime_max,dt);

%% GMM fitting
GMModels = cell(n_peak,1); % Preallocation
Mu=cell(n_peak,1);
options = statset('MaxIter',2*10^3,'TolFun',10^-10);
for j =n_peak% 1:n_peak_max
    GMModels{j} = fitgmdist([list_alpha'/2,list_Y'],j,'SharedCovariance',true,'Options',options);
    Mu{j} = GMModels{j}.mu;
    Mu{j}(:,2)=Mu{j}(:,2);
end

C=[Mu{j},GMModels{j}.ComponentProportion'];
%% initial estimates
xa=Y_alpha_stat_scan.list_alpha_scan(1,:);
xD=Y_alpha_stat_scan.list_D_scan(:,1);
list_init=zeros(size(C,1),size(C,2));
for j=1:n_peak
    ua=abs(Y_alpha_stat_scan.alpha_mean_th-C(j,1)*2);
    uD=abs(Y_alpha_stat_scan.Y_mean_th-C(j,2));
    u=(ua/mean(ua(:)))+(uD/mean(uD(:)));

    [idx_D,idx_a]=find(u==min(u(:)));

    list_init(j,1)=xa(idx_a)/2;
    list_init(j,2)=xD(idx_D);
    list_init(j,3)=C(j,3);
end
list_init=list_init';
x0=list_init(:)';
%% least-square fitting
[x,y]=meshgrid(X,Y);
xdata(:,:,1)=x;
xdata(:,:,2)=y;

lb=[border_alpha(1)/2,border_D(1),0.1]; %  lower bounds for fitting
ub=[1.8/2,border_D(2),1]; %  upper bounds for fitting
lb=repmat(lb,1,n_peak);
ub=repmat(ub,1,n_peak);
p=1/n_peak*ones(1,n_peak);
%% fit with FMinCon
%             Aeq=[zeros(n_peak,2),ones(n_peak,1)]';
%             Aeq=Aeq(:)';
%             beq=1;
%             % options = optimoptions('fmincon','FunctionTolerance',10^-14,'StepTolerance',10^-14);%,'MaxIterations',1000,'MaxFunctionEvaluations',1000);
%             options = optimoptions('fmincon','FunctionTolerance',10^-3,'StepTolerance',10^-3);%,'MaxIterations',1000,'MaxFunctionEvaluations',1000);
% 
%             fun_opt=@(x)sum((joint_pdf_alpha_Y_n_point_interp_wrapper_v2(x,xdata,Y_alpha_stat_scan,fit_factor)-mat_res'*fit_factor).^2,'all');
%             fit_param = fmincon(fun_opt,x0,[],[],Aeq,beq,lb,ub,[],options);
%% fit with lsqcurvefit
opts = optimset('Display','off');
[fit_param,resnorm] = lsqcurvefit(@(x,xData)...
    joint_pdf_alpha_Y_n_point_interp_wrapper_v3(x,xData,Y_alpha_stat_scan),x0,xdata,mat_res',lb,ub,opts);%
%%
fit_H_D_p=[fit_param(1:3:end);fit_param(2:3:end);fit_param(3:3:end)];
fitted_mat=joint_pdf_alpha_Y_n_point_interp_wrapper_v3(fit_param,xdata,Y_alpha_stat_scan);
end