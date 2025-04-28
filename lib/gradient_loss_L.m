function dL=gradient_loss_L(beta,X,Y,n)
% dL=gradient_loss_L(beta,X,Y,n)
% The gradient ofloss function in original problem (O) (average by sample size)
% beta: (d*1) the coefficient in vector form
% X: (n×d)
% Y: (n×1)
% 
% dL:(d×1)

dL = -2*X'*(Y -X*beta');
dL=dL/n;

end