function [fit_alpha_D_p,list_D,list_alpha,resnorm] = fit_joint_PDF_alpha_D_from_trajectory(Trajectory,dt,lagtime_max,n_peak,n_bin,n_bin_scan)
%%
% Input:
% Trajectory : N x M array of M trajectories of N steps
% lagtime_max: lag_time max for  fitting (tau_max in the article)
% n_peak: number of components for fitting the joint pdf of estimated
% \hat{alpha}, \hat{D}

% Output:


%%
ffont=15; % fontsize of figures
LW=2; % line width of the curves and meshes

color=lines;
color(1:3,:)=[216,27,96;30,136,229;255,193,7]/255; % colors for displaying fitted components
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
parameter_fit.n_bin_scan=n_bin_scan;
parameter_fit.N=N;
parameter_fit.lagtime_max=lagtime_max;
parameter_fit.dt=dt;
%% fitting
[fit_H_D_p,Y_alpha_stat_scan,mat_res,fitted_mat,x_Y,x_alpha,resnorm] = fun_fit_joint_PDF_alpha_Y_est_nocorr_v3(list_alpha,list_D,parameter_fit);

fit_alpha_D_p=fit_H_D_p;
fit_alpha_D_p(1,:)=fit_H_D_p(1,:)*2;

X=(x_Y(1:end-1)+x_Y(2:end))/2;%Y

Y=(x_alpha(1:end-1)+x_alpha(2:end))/2;%alpha

%% plotting

histogram2('XBinEdges',x_Y,'YBinEdges',x_alpha,'BinCounts',mat_res,'FaceColor','flat');
hold on

xlabel('$$\ln(2\hat{D})$$','Interpreter','latex')
ylabel('$$\hat{\alpha}$$','Interpreter','latex')
zlabel('joint pdf')


hold on
n_bin_plot=2*numel(x_Y);
x_Y_plot=linspace((min(list_Y(~isinf(list_Y)))),(max(list_Y)),n_bin_plot);

%             x_alpha_plot=linspace(min(list_alpha),max(list_alpha),n_bin_plot);

x_alpha_plot=linspace(-0.5,1.8,n_bin_plot);

[x_plot1,y_plot1]=meshgrid(x_Y_plot,x_alpha_plot);

xData(:,:,1)=x_plot1;
xData(:,:,2)=y_plot1;


fitted_mat= joint_pdf_alpha_Y_n_point_interp_wrapper_v3(fit_H_D_p(:)',xData,Y_alpha_stat_scan);




fitted_mat(fitted_mat<min(mat_res(mat_res>0)))=nan;
mesh(x_Y_plot,x_alpha_plot,fitted_mat,'facealpha',0,'EdgeColor',[255,29,141]/255,'LineWidth',2)


view(75,45)

set(findall(gcf,'-property','FontSize'),'FontSize',ffont)
drawnow


%%
figure('Color','white')

[Y_alpha_stat_fitted] = covar_Y_alpha_fBm_noise_mean_var_nocorr_v1(2*fit_H_D_p(1,:),fit_H_D_p(2,:),N,lagtime_max,dt);


%%
x_D=logspace(log10(min(list_D(list_D>0))),log10(max(list_D)),n_bin);

% alpha versus D
subplot(3,3,[4,5,7,8]);

semilogy(list_alpha,list_D,'marker','+','Color','black','linestyle','none','DisplayName','Sim')

% [x_D_est,x_alpha_est]=meshgrid(x_D_plot,x_alpha_plot);

n_bin_plot=150;
x_D_plot=logspace(log10(min(list_D(list_D>0))),log10(max(list_D)),n_bin_plot);

%             x_alpha_plot=linspace(min(list_alpha),max(list_alpha),n_bin_plot);
x_alpha_plot=linspace(-0.5,1.8,n_bin_plot);

[x_plot,y_plot]=meshgrid(x_D_plot,x_alpha_plot);

hold on
for peak=1:n_peak
    fit_factor=1;
    [fitted_mat_contour] = joint_pdf_alpha_D_n_point_interp_v2(x_plot,y_plot,[fit_H_D_p(1,peak),fit_H_D_p(2,peak)],fit_H_D_p(3,peak),Y_alpha_stat_scan,1);
    % fitted_mat_contour=joint_pdf_alpha_Y_n_point_interp_wrapper_v2([list_fit_H_D_p(1,peak),list_fit_H_D_p(2,peak),list_fit_H_D_p(3,peak)],xdata_plot,Y_alpha_stat_scan,fit_factor)/fit_factor;
    %  contour(x_alpha,x_Y,fitted_mat',2,'color',color(peak,:),'levels',10^-1,'linewidth',2)%, axis equal
    contour(x_alpha_plot,x_D_plot,fitted_mat_contour'>5*10^-1,2,'color',color(peak,:),'linewidth',2,'DisplayName',['Comp. ',num2str(peak)]);%, axis equal
end

% xlim([min(list_alpha),max(list_alpha)])
ylim([min(list_D),max(list_D)])
xlim([-0.5,1.8])
xlabel('$$\hat{\alpha}$$','Interpreter','latex')
ylabel('$$\hat{D}$$','Interpreter','latex')
legend('show','Location','southeast')
legend('box','off')
% histogram alpha
subplot(3,3,1:2)
h_alpha=histogram(list_alpha,x_alpha,'Normalization','pdf','FaceColor','none');
h_alpha.DisplayName='Sim';
hold on
sum_y=zeros(1,numel(x_alpha_plot));
for ni=1:n_peak
    y_fit_i=fit_H_D_p(3,ni)*normpdf(x_alpha_plot,Y_alpha_stat_fitted.alpha_mean_th(ni),Y_alpha_stat_fitted.alpha_var_th(ni)^0.5);
    sum_y=sum_y+y_fit_i;
    plot(x_alpha_plot,y_fit_i,'LineWidth',LW,'Color',color(ni,:),'DisplayName',['Comp. ',num2str(ni)])
end
plot(x_alpha_plot,sum_y,'LineWidth',LW,'Color','black','linestyle','--','DisplayName','Sum')
xlim([-0.5,1.8])

legend('show','Location','northeast')
legend('box','off')
% histogram D
subplot(3,3,[6,9])


hD=histogram(list_D,x_D,'Normalization','pdf','FaceColor','none');
hD.DisplayName='sim';

hold on
sum_y=zeros(1,numel(x_D_plot));
for ni=1:n_peak
    y_fit_i=2*fit_H_D_p(3,ni)*lognpdf(x_D_plot*2,Y_alpha_stat_fitted.Y_mean_th(ni),Y_alpha_stat_fitted.Y_var_th(ni)^0.5);
    sum_y=sum_y+y_fit_i;
    plot(x_D_plot,y_fit_i,'LineWidth',LW,'Color',color(ni,:),'DisplayName',['Comp. ',num2str(ni)])
end
plot(x_D_plot,sum_y,'LineWidth',LW,'Color','black','linestyle','--','DisplayName','Sum')

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