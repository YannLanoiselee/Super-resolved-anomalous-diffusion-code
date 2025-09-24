function [pdf_joint] = copula_beta_Rice_D_alpha_true_mean_var_ncomp(xb,xg,beta_mean,beta_std,beta_min,beta_max,mu,sigma,rho,p)
% see https://en.wikipedia.org/wiki/Copula_(probability_theory)
% rows are beta and columns are D
%% beta
% Xb=1/(beta_max-beta_min)*betacdf((xb-beta_min)/(beta_max-beta_min),a,b);% CDF beta


delta_beta=(beta_max-beta_min);
pdf_joint=zeros(numel(xb,xg));
for n_comp=1:numel(beta_mean)
mu_beta=(beta_mean(n_comp)-beta_min)/delta_beta;
phi=delta_beta^2*mu_beta*(1-mu_beta)/beta_std(n_comp)^2-1;
a=mu_beta*phi;
b=(1-mu_beta)*phi;

pdf_beta=@(x)exp(-betaln(a,b))*x.^(a-1).*(1-x).^(b-1);
if a>100 && b>100
    mean_X=a/(a+b)*(beta_max-beta_min)+beta_min;
std_X=(beta_max-beta_min)/(a+b)*sqrt(a*b/(a+b+1));
Xb=normcdf(xb,mean_X,std_X);
pdfb=normpdf(xb,mean_X,std_X);
else
xb_rescaled=(xb-beta_min)/(beta_max-beta_min);
% Xb=betacdf(xb_rescaled,a,b);% CDF beta
%pdfb=1/(beta_max-beta_min)*betapdf(xb_rescaled,a,b); % PDF beta
pdfb=1/(beta_max-beta_min)*pdf_beta(xb_rescaled); % PDF beta
Xb=cumtrapz(xb_rescaled,pdfb);% CDF beta
end




%pdfb=betapdf((xb-beta_min)/(beta_max-beta_min),a,b); % PDF beta



% 
% %% lognormal
% Xg = logncdf(xg,mu(n_comp),sigma(n_comp)); % CDF lognormal
% pdfg=lognpdf(xg,mu(n_comp),sigma(n_comp)); % PDF  lognormal Gaussian
% % %% normal
% % Xg = normcdf(xg,mu,sigma); % CDF Gaussian
% % pdfg=normpdf(xg,mu,sigma); % PDF Gaussian
% 
%% RIcian
rad = makedist('Rician','s',mu(n_comp),'sigma',sigma(n_comp));
pdfg = pdf(rad,xg);
Xg = cdf(rad,xg);
% pdfg = pdf('Rician',xg,mu(n_comp),sigma(n_comp));
% Xg = cdf('Rician',xg,mu(n_comp),sigma(n_comp));
%%
[u1,u2] = meshgrid(Xb,Xg); % grid for computing copula

[p1,p2] = meshgrid(pdfb,pdfg);

y = copulapdf('Gaussian',[u1(:),u2(:)],rho(n_comp)); % Gaussian copula

pdf_joint=pdf_joint+p(n_comp)*reshape(y.*p1(:).*p2(:),numel(xb),numel(xg)); % joint PDF 
end
end

