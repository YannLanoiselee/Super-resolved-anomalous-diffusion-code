function [A] = R_matrix(N,n)
% Matrix to get TASMD through quadratic form
%% diagonal terms
k=1:N;
mi=zeros(N,1);
if n<=N/2
% mi(k<=n | k>=N-n)=1;
% mi(n<k & k<N-n)=2;
mi(k<=n | k>N-n)=1;
mi(n<k & k<=N-n)=2;

else
% mi(k<=N-n | k>=n)=1;
% mi(N-n<k & k<n)=0;
mi(k<=N-n | k>n)=1;
mi(N-n<k & k<=n)=0;

end
%% out of diagonal terms
idx=abs(k-k')==n;
%%
A=(diag(mi)-idx)/(N-n);
end