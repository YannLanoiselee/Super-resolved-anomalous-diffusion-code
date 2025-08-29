function [fit_alpha_mean_std_D_mean_std_rho_p,list_D,list_alpha,resnorm] = fit_joint_PDF_alpha_D_continuous_distribution_v2(Trajectory,dt,parameter_analysis)
%%
% Input:
% Trajectory : N x M array of M trajectories of N steps
% lagtime_max: lag_time max for  fitting (tau_max in the article)
% n_peak: number of components for fitting the joint pdf of estimated
% \hat{alpha}, \hat{D}




% Output:

%%
lagtime_max=parameter_analysis.lagtime_max; % maximum lag time (tau_max in the article) for fitting TAMSD
n_peak=parameter_analysis.n_peak; % number of component of the true parameters continuous distribution p(alpha,D)
n_bin=parameter_analysis.n_bin;% number of bins of the histograms of alpha and y for fitting
n_bin_scan=parameter_analysis.n_bin_scan;% number of bins for precomputing the moments of estimated parameters as a function of the true ones before interpolation.
n_bin_true=parameter_analysis.n_bin_true;% number of bins to discretize the theoretical true distribution 

pdf_marginal_D=parameter_analysis.pdf_marginal_D;% PDF of the true diffucion coefficient. Two options: 'lognormal' and 'rice'. Both distribution have strictly positive support. 


%%
ffont=15; % fontsize of figures
LW=2; % line width of the curves and meshes

color=lines;
color(1:3,:)=[216,27,96;30,136,229;255,193,7]/255; % colors for displaying fitted components

%%

%%
warning('off','MATLAB:griddedInterpolant:CubicUniformOnlyWarnId');
%%

N=size(Trajectory,1);
M=size(Trajectory,2);

[list_TAMSD ] = TAMSD_Trajectory( Trajectory,lagtime_max );
list_alpha=sum(repmat(log(1:lagtime_max)',1,M).*log(list_TAMSD./repmat(list_TAMSD(1,:),lagtime_max,1)),1)...
    ./sum(log(1:lagtime_max).^2);

list_meanlagTASMD=mean(log(list_TAMSD),1);
B=mean(log((1:lagtime_max)'*dt));
list_Y=list_meanlagTASMD-list_alpha*B;
list_D=1/2*exp(list_Y);

%% D vs alpha

parameter_fit=struct;
parameter_fit.n_peak=n_peak;
parameter_fit.n_bin=n_bin;
parameter_fit.n_bin_true=n_bin_true;

parameter_fit.n_bin_scan=n_bin_scan;
parameter_fit.N=N;
parameter_fit.lagtime_max=lagtime_max;
parameter_fit.dt=dt;
alpha_min=0;
alpha_max=2;
%% fitting
if strcmp(pdf_marginal_D,'beta_lognormal')
[fit_alpha_mean_std_D_mean_std_rho_p,Y_alpha_stat_scan,mat_res,fitted_mat,x_Y,x_alpha,resnorm] = fun_fit_joint_PDF_alpha_Y_est_nocorr_continuous_beta_lognorm_v1(list_alpha,list_D,parameter_fit);
elseif strcmp(pdf_marginal_D,'beta_rice')
[fit_alpha_mean_std_D_mean_std_rho_p,Y_alpha_stat_scan,mat_res,fitted_mat,x_Y,x_alpha,resnorm] = fun_fit_joint_PDF_alpha_Y_est_nocorr_continuous_beta_Rice_v1(list_alpha,list_D,parameter_fit);

elseif strcmp(pdf_marginal_D,'discrete')

end    
%%*
alpha_mean_fit=fit_alpha_mean_std_D_mean_std_rho_p(:,1);
alpha_std_fit=fit_alpha_mean_std_D_mean_std_rho_p(:,2);
D_mean_fit=fit_alpha_mean_std_D_mean_std_rho_p(:,3);
D_std_fit=fit_alpha_mean_std_D_mean_std_rho_p(:,4);
rho_fit=fit_alpha_mean_std_D_mean_std_rho_p(:,5);
p_fit=fit_alpha_mean_std_D_mean_std_rho_p(:,6);

delta_alpha=2;

mu_fit=(alpha_mean_fit-alpha_min)./delta_alpha;
phi_fit=delta_alpha^2.*mu_fit.*(1-mu_fit)./alpha_std_fit.^2-1;
a_fit=mu_fit.*phi_fit;
b_fit=(1-mu_fit).*phi_fit;



if strcmp(pdf_marginal_D,'beta_lognormal')
D_mulognorm_fit=log(D_mean_fit.^2./sqrt(D_mean_fit.^2+D_std_fit.^2));
D_stdlognorm_fit=log(1+D_std_fit.^2./D_mean_fit.^2);
[x_alpha_true_fit,x_D_true_fit] = find_bounds_true_alpha_D_beta_logn_v2(alpha_mean_fit,alpha_std_fit,alpha_min,alpha_max,D_mulognorm_fit,D_stdlognorm_fit,100,100);

elseif strcmp(pdf_marginal_D,'beta_rice')
    [x_alpha_true_fit,x_D_true_fit] = find_bounds_true_alpha_D_beta_Rice_v1(alpha_mean_fit,alpha_std_fit,alpha_min,alpha_max,D_mean_fit,D_std_fit,100,100);

end


X=(x_Y(1:end-1)+x_Y(2:end))/2;%Y

Y=(x_alpha(1:end-1)+x_alpha(2:end))/2;%alpha

%% plot joint pdf estimated vs fitting

histogram2('XBinEdges',x_Y,'YBinEdges',x_alpha,'BinCounts',mat_res,'FaceColor','flat');
hold on

xlabel('$$\ln(2\hat{D})$$','Interpreter','latex')
ylabel('$$\hat{\alpha}$$','Interpreter','latex')
zlabel('joint pdf')


hold on
n_bin_plot=2*numel(x_Y);
x_Y_plot=linspace((min(list_Y(~isinf(list_Y)))),(max(list_Y)),n_bin_plot);

 x_alpha_plot=linspace(min(list_alpha),max(list_alpha),n_bin_plot);

%x_alpha_plot=linspace(-0.5,1.8,n_bin_plot);

[x_plot1,y_plot1]=meshgrid(x_Y_plot,x_alpha_plot);

xData(:,:,1)=x_plot1;
xData(:,:,2)=y_plot1;


% fitted_mat= joint_pdf_alpha_Y_n_point_interp_wrapper_v3(fit_H_D_p(:)',xData,Y_alpha_stat_scan);
% 
fit_param_fin=fit_alpha_mean_std_D_mean_std_rho_p';
if strcmp(pdf_marginal_D,'beta_lognormal')
 fitted_mat = joint_pdf_est_alpha_Y_true_alpha_D_betalogn_ncomp_wrapper_v1(fit_param_fin(:)',xData,n_bin_true,Y_alpha_stat_scan);

elseif strcmp(pdf_marginal_D,'beta_rice')
   fitted_mat = joint_pdf_est_alpha_Y_true_alpha_D_betaRice_ncomp_wrapper_v1(fit_param_fin(:)',xData,n_bin_true,Y_alpha_stat_scan);
end

fitted_mat(fitted_mat<min(mat_res(mat_res>0)))=nan;
mesh(x_Y_plot,x_alpha_plot,fitted_mat,'facealpha',0,'EdgeColor',[255,29,141]/255,'LineWidth',2)


view(75,45)

set(findall(gcf,'-property','FontSize'),'FontSize',ffont)
drawnow

%% figure end results
%%
color=[1,0,1;... % true
    0,1,1;... % est
    0.7,0.7,0.7;...% histogram true
    139/255,0,0]; % data points true
put_label=0;
color(1:2,:)=color(1:2,:)*2/3;

color(1,:)=2/3*[0,1,1];
color(2,:)=[255,29,141]/255;
color(4,:)=[132,82,129]/255;
color_hist_est=[1,1,0]*3/4;


mix=repmat((0:63)'/63,1,3);
mix=mix(end:-1:1,:).^2;
base_color=repmat([1,0,1],64,1);
color_hist_true_grad=( base_color+(mix).*(ones(64,3)- repmat([1,0,1],64,1)));

% base_color=repmat(color(2,:),64,1);
% color_hist_true_grad=( base_color+(mix).*(ones(64,3)- repmat(color(2,:),64,1)));


% color_hist_true_grad=exp(-color_hist_true_grad);
color_hist_true_grad=color_hist_true_grad;

n_bin_est_circ=200; % number of pixels of the contour plot of the fitted joint pdf of estimated parameters (the higher the smoother)
perc_peak=0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_summary=figure('color','white');
%% 2d plot alpha vs D
subplot(3,3,[4,5,7,8]); 
colormap(color_hist_true_grad)

%%
semilogy(list_alpha,list_D,'marker','.','Color','black','linestyle','none','DisplayName','estimated')
hold on
x_D_est_circ=1/2*exp(linspace(min(list_Y),max(list_Y),n_bin_est_circ));
x_alpha_est_circ=linspace(min(list_alpha),max(list_alpha),n_bin_est_circ);

x_D_est_plot_int=(x_D_est_circ(2:end)-x_D_est_circ(1:end-1))./(log(x_D_est_circ(2:end))-log(x_D_est_circ(1:end-1)));%D
x_alpha_est_plot_int=(x_alpha_est_circ(1:end-1)+x_alpha_est_circ(2:end))/2;%alpha

x_D_est=1/2*exp(linspace(min(list_Y),max(list_Y),n_bin));
x_alpha_est=linspace(min(list_alpha),max(list_alpha),n_bin);
x_D_est_plot=x_D_est;%D
x_alpha_est_plot=x_alpha_est;%alpha

%%
 [x,y]=meshgrid(x_D_est_plot,x_alpha_est_plot);% D,alpha

if strcmp(pdf_marginal_D,'beta_lognormal')
[mat_D_alpha_est_fit] =joint_pdf_alpha_Y_est_fct_true_alpha_D_betalognorm_ncomp_v1(log(2*x),y,n_bin_true,alpha_mean_fit,alpha_std_fit,alpha_min,alpha_max,D_mulognorm_fit,D_stdlognorm_fit,rho_fit,p_fit,Y_alpha_stat_scan);

elseif strcmp(pdf_marginal_D,'beta_rice')
[mat_D_alpha_est_fit] =joint_pdf_alpha_Y_est_fct_true_alpha_D_betaRice_ncomp_v1(log(2*x),y,n_bin_true,alpha_mean_fit,alpha_std_fit,alpha_min,alpha_max,D_mean_fit,D_std_fit,rho_fit,p_fit,Y_alpha_stat_scan);


end
mat_D_alpha_est_fit=mat_D_alpha_est_fit./x;
mat_D_alpha_est_fit_plot=mat_D_alpha_est_fit;

[X_circ,Y_circ]=meshgrid(x_D_est_plot,x_alpha_est_plot);
[Xq_circ,Yq_circ]=meshgrid(x_D_est_plot_int,x_alpha_est_plot_int);


dx_D_est_circ=diff(x_D_est_circ);
dx_alpha_est_circ=diff(x_alpha_est_circ);
[dx_int,dy_int]=meshgrid(dx_D_est_circ,dx_alpha_est_circ);




%
 n_bin_true_circ=200;
x_D_true_circ=logspace(log10(min(list_D)),log10(max(list_D)),n_bin_true_circ);
x_alpha_true_circ=linspace(max([0,min(list_alpha)]),max(list_alpha),n_bin_true_circ);
x_D_true_circ_mean=(x_D_true_circ(2:end)-x_D_true_circ(1:end-1))./(log(x_D_true_circ(2:end))-log(x_D_true_circ(1:end-1)));
x_alpha_true_circ_mean=1/2*(x_alpha_true_circ(2:end)+x_alpha_true_circ(1:end-1));
[x_D_true_circ_mat,x_alpha_true_circ_mat]=meshgrid(x_D_true_circ_mean,x_alpha_true_circ_mean);

dx_D_true_circ=diff(x_D_true_circ);
dx_alpha_true_circ=diff(x_alpha_true_circ);
[dx_int_true,dy_int_true]=meshgrid(dx_D_true_circ,dx_alpha_true_circ);
x_D_true_circ=((x_D_true_circ(2:end))-(x_D_true_circ(1:end-1)))./(log(x_D_true_circ(2:end))-log(x_D_true_circ(1:end-1)));
x_alpha_true_circ=1/2*(x_alpha_true_circ(2:end)+x_alpha_true_circ(1:end-1));

[pdf_true_alpha_D_circ] = copula_beta_Rice_D_alpha_true_mean_var_ncomp(x_alpha_true_circ,x_D_true_circ,alpha_mean_fit,alpha_std_fit,alpha_min,alpha_max,D_mean_fit,D_std_fit,rho_fit,p_fit);
pdf_true_alpha_D_circ=pdf_true_alpha_D_circ.*x_D_true_circ_mat';
% imagesc(x_alpha_true_circ,x_D_true_circ,(pdf_true_alpha_D_circ/max(pdf_true_alpha_D_circ(:))).^(4/4) )
% imagesc(x_alpha_true_circ,x_D_true_circ,(pdf_true_alpha_D_circ.*dx_int_true.*dy_int_true)')
imagesc(x_alpha_true_circ,x_D_true_circ,(pdf_true_alpha_D_circ))

% % I=pdf_true_alpha_D_circ.*dx_int_true.*dy_int_true;
% %  int1_true=min(I(I>0));
% % 
% % int2_true=max(I(I>0));
% % 
% % val_true_plot_contour=perc_peak*(int2_true-int1_true)+int1_true;
% % contour(x_alpha_true_circ,x_D_true_circ,I>val_true_plot_contour,1/2*[1,1],'color',color(2,:),'linewidth',2,'DisplayName',['inferred']);%, axis equal 
% % % contour(x_alpha_true_circ,x_D_true_circ,I>0.1^,1/2*[1,1],'color',color(2,:),'linewidth',2,'DisplayName',['inferred']);%, axis equal 


%%


mat_D_alpha_est_fit_circ= interp2(X_circ,Y_circ,mat_D_alpha_est_fit,Xq_circ,Yq_circ,'cubic');

I=mat_D_alpha_est_fit_circ.*dx_int.*dy_int;
int1_I=min(I(I>0));
int2_I=max(I(I>0));
val_alpha_D_est2=int1_I+perc_peak*(int2_I-int1_I);%(int2_I-int1_I)/(log(int2_I)-log(int1_I));%sqrt(M)/10;%int1+1/n_bin_est*(int2-int1);
contour(x_alpha_est_plot_int,x_D_est_plot_int,I'>val_alpha_D_est2,1/2*[1,1],'color',color(1,:),'linewidth',2,'DisplayName',['fitted']);%, axis equal 
%%
%
hold on
plot(list_alpha,list_D,'marker','.','Color','black','linestyle','none','DisplayName','estimated')
%set(gca,'YScale','log')
% set(gca,'ColorScale','log')
ylim([min(list_D),max(list_D)])
xlim([min(list_alpha),max(list_alpha)])

xlabel('Anomalous exponent','FontSize',ffont)
ylabel('Diffusion coefficient','FontSize',ffont)

legend('$$\lbrace\hat{\alpha},\hat{D}\rbrace$$','fitted','inferred','interpreter','latex','box','off','location','southeast')

legend('show','box','off')

%% histogram of alpha
subplot(3,3,1:2)
h_alpha=histogram(list_alpha,x_alpha_est,'Normalization','pdf','FaceColor',color_hist_est,'FaceAlpha',0.2,'DisplayName','estimated');
h_alpha.DisplayName='estimated';
hold on
%pdf_alpha_est_fit=sum(mat_D_alpha_est_fit.*dx_est,2);
pdf_alpha_est_fit=trapz(x_D_est_plot,mat_D_alpha_est_fit,2);


% pdf_alpha_est_fit=pdf_alpha_est_fit./trapz(x_alpha_est_plot,pdf_alpha_est_fit);
dx=(x_alpha_est_plot(2)-x_alpha_est_plot(1))/2;
plot(x_alpha_est_plot,pdf_alpha_est_fit,'linewidth',2,'Color',color(1,:),'displayname','fitted')

pdf_true_alpha_fit=zeros(1,numel(x_alpha_true_fit));
for n_comp=1:n_peak

    pdf_true_alpha_fit=pdf_true_alpha_fit+p_fit(n_comp)*1/2*betapdf(x_alpha_true_fit/2,a_fit(n_comp),b_fit(n_comp));

end
plot(x_alpha_true_fit,pdf_true_alpha_fit,'linewidth',2,'Color',color(2,:),'displayname','inferred')

xlim([min(list_alpha),max(list_alpha)])
if put_label==1
% ylabel('$$p_{\hat{\alpha}}(\hat{\alpha})$$','Interpreter','latex','FontSize',ffont)
ylabel('Marginal PDF','FontSize',ffont)
end
% set(gca,'yScale','log')
%legend('show','box','off')
%% histogram of D
subplot(3,3,[6,9])
hD=histogram(list_D,x_D_est,'Normalization','pdf','FaceColor',color_hist_est,'FaceAlpha',0.2);
hD.DisplayName='sim';
hold on


%pdf_D_est_fit=sum(mat_D_alpha_est_fit.*dy_est,1);
pdf_D_est_fit=trapz(x_alpha_est_plot,mat_D_alpha_est_fit,1);
% pdf_D_est_fit=pdf_D_est_fit./trapz(x_D_est_plot,pdf_D_est_fit);
hold on
plot(x_D_est_plot,pdf_D_est_fit,'linewidth',2,'Color',color(1,:))

pdf_D_true_fit=zeros(1,numel(x_D_true_fit));
for n_comp=1:n_peak
    if strcmp(pdf_marginal_D,'beta_lognormal')

pdf_D_true_fit=pdf_D_true_fit+p_fit(n_comp)*lognpdf(x_D_true_fit,D_mulognorm_fit(n_comp),D_stdlognorm_fit(n_comp));
elseif strcmp(pdf_marginal_D,'beta_rice')
rad = makedist('Rician','s',D_mean_fit(n_comp),'sigma',D_std_fit(n_comp));

pdf_D_true_fit=pdf_D_true_fit+p_fit(n_comp)*pdf(rad,x_D_true_fit);

end
    end
plot(x_D_true_fit,pdf_D_true_fit,'linewidth',2,'Color',color(2,:))
xlim([min(list_D),max(list_D)])
set(gca,'xScale','log')
set(gca,'view',[90,90])
set(gca,'xDir','reverse')
if put_label==1
% ylabel('$$p_{\hat{D}}(\hat{D})$$','Interpreter','latex','FontSize',ffont)
ylabel('Marginal PDF','Interpreter','latex','FontSize',ffont)
end
% set(gca,'yScale','log')

%% legend

subplot(3,3,3)
histogram([nan,nan],[0,1],'Normalization','pdf','FaceColor',color_hist_est,'FaceAlpha',0.2,'DisplayName','estimated');
hold on
%histogram([nan,nan],[0,1],'Normalization','pdf','FaceColor',color_hist_true,'FaceAlpha',0.3,'DisplayName','ground truth');
%plot([nan,nan],[nan,nan],'DisplayName','estimated')
%hold on
%plot([nan,nan],[nan,nan],'DisplayName','ground truth')
plot([nan,nan],[nan,nan],'Color',color(1,:),'DisplayName','fitted')
plot([nan,nan],[nan,nan],'Color',color(2,:),'DisplayName','inferred')
legend('show','box','off','location','best')
ylim([5,10])
axis off

%%

set(findall(gcf,'-property','FontSize'),'FontSize',ffont)

%% %


xlim([min(list_alpha),max(list_alpha)])
% legend('show')
set(gca,'xScale','log')
set(gca,'view',[90,90])
set(gca,'xDir','reverse')

%                           set(gca,'yDir','reverse')
hold on
xlim([min(list_D),max(list_D)])
set(findall(gcf,'-property','FontSize'),'FontSize',ffont)
drawnow

end