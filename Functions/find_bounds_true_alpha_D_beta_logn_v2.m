function [x_beta,x_logn] = find_bounds_true_alpha_D_beta_logn_v2(beta_mean,beta_std,beta_min,beta_max,mu,sigma,n_bin_beta,n_bin_logn)
% x_beta=linspace(max([min(beta_mean-5*beta_std),beta_min]),min([max(beta_mean+5*beta_std),beta_max]),n_bin_beta);
% 
% x_logn=exp(linspace(min(mu-5*sigma),max(mu+5*sigma),n_bin_logn));

% 
%  min_b=   max([beta_mean-5*beta_std;beta_min*ones(1,numel(beta_mean))],[],1);
% max_b=   min([beta_mean+5*beta_std;beta_max*ones(1,numel(beta_mean))],[],1);
% 
% test=min_b-max_b';

x_beta=linspace(max([min(beta_mean-5*beta_std),beta_min]),min([max(beta_mean+5*beta_std),beta_max]),n_bin_beta);

x_logn=exp(linspace(min(mu-5*sigma),max(mu+5*sigma),n_bin_logn));
end