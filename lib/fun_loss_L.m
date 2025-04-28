function L=fun_loss_L(beta,X,Y,n)
% L=fun_loss_L(beta,X,Y,n)
% The loss function in original problem (average by sample size)
% beta: (1×d)the coefficient in vector form
% X: (n×d)
% Y: (n×1)

L = norm(Y-X*beta',2)^2;
L=L/n;

end