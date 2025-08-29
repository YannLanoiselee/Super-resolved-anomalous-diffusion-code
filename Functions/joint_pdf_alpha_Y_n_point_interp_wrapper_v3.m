function [mat] = joint_pdf_alpha_Y_n_point_interp_wrapper_v3(x,xData,Y_alpha_stat)
%wrapper for fitting joint pdf of alpha and Dbased on the precomputed
%moments of  estimated Y and alpha conditional on the true Y and alpha
list_H_D=reshape([x(1:3:end),x(2:3:end)],numel(x)/3,2);
p=x(3:3:end);
[mat] = joint_pdf_alpha_Y_n_point_interp_v4(xData(:,:,1),xData(:,:,2),list_H_D,p,Y_alpha_stat);
end