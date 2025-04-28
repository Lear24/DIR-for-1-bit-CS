function  P = fun_penalty_UP(W,U_i,m,d)
% P = fun_penalty_UP(W,U_i,m,d)
% Calculate the penalty function values based on the transfer matrix U 
% 
% U_i: U_i as determined by the transmitted vector U. not change in local
% interation.
% 
% W: local coefficient.
% m: number of client
% d: dimension of coefficient


[~,I2] = util_matrix_I(d);
P =1/2 - trace(W*U_i)/trace(I2*W*I2)/2/m -1/2/m;
end