function [W,norm_W,norm_W_1]=util_O2W_coefficient(beta)
% [W,u,u1]=util_O2W_coefficient(beta)
% 
% On certain clien, tranfer the coefficent in problem (O) to the
% corresponding coefficent in problem (W).
% 
% beta: coefficient in (O).
% 
% W: the  (d+1 × d+1) variable define in (W).
% norm_W: the norm of beta, 
% norm_W_1: trace(W'W)-1, equivalently the norm of beta when W is rank-1 matrix.



norm_W = norm(beta,2)^2;



W = [beta';1]*[beta,1];

[V,D] = eig(W);
D(D > -0.0001 & D < 0)=0;
D(D < 0.0001 & D > 0)=0;
L=V*(D^0.5);
norm_W_1 = trace(L'*L)-1;
end