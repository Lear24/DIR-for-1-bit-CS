function dP=gradient_vP(beta,ome_i,m)
% dP=gradient_vP(beta,xi_i, m,d)
% The gradient of the original problem in a distributed setting
% 
% beta: (1×d)
% ome_i: (d×d) matrix the transmitting vector, assigned to client i.
% m: number of client



dP=-2*(ome_i -  (beta'*beta/norm(beta)^2) * ome_i )*beta'/norm(beta)^2;
dP=dP/m;  %其实传输的向量本身就以及包含了这个系数


end