function [vdP,dP] = gradient_UP(W,U_i,m,d)
% [vdP,dP] = gradient_UP(W,U_i,m,d)
% Calculate the gradient of penalty function  based on the transfer matrix U 
% 
% U_i: U_i as determined by the transmitted vector U. not change in local
% interation.
% 
% W: local coefficient.
% m: number of client
% d: dimension of coefficient



d1 = d+1;
[~,I2] = util_matrix_I(d);
tr_i = trace(I2*W*I2);
[V,D] = eig(W);
D(D < 0.0001 & D > 0)=0;
D(D > -0.0001 & D < 0)=0;
A = V*D^(1/2);

dP = zeros(d1,d1);


% tr(W,U_i)/(tr(W)-1)
% dP = dP - 2*U_i*A/tr_i/2/m;
% dP = dP + 2*trace(W*U_i)*A/tr_i/tr_i/2/m;

% the gradient according to centering objective function.
dP = dP - 2*U_i*A/tr_i/m;
dP = dP + 2*trace(W*U_i)*A/tr_i/tr_i/m;
dP(end) = 0;


vdP =zeros(d1*d1,1);
for i = 1:d1*d1
    vdP(i) = dP(i);
end


end