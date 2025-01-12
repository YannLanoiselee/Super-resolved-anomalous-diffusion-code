function [cov_mat] = cov_mat_fBm_dt(N,H,D,dt)
k=abs((1:N)'-(1:N));
% cov_mat=AutoCovariance(k,H,D);
%     function [ Result ] = AutoCovariance(k,H,D)
%         %Function that gives autocovariance of a fractionnal gaussian
%         %noise
%         % with k l'indice du saut et H l'exposant de Holder
%         Result=D*((abs(k-1)).^(2*H)-2*(abs(k)).^(2*H)+(abs(k+1)).^(2*H));
%     end
if H~=0
cov_mat=D*dt^(2*H)*(repmat((1:N),N,1).^(2*H)+repmat((1:N),N,1)'.^(2*H)-abs(k).^(2*H));
else
cov_mat=D*dt^(2*H)*eye(N,N);
    
end


end