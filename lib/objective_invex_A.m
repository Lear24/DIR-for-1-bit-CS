function [J,grad]= objective_invex_A(para_A,U_i,X,Y,m,n,d,lambda)
% [J,grad]= onebitcostfunction(beta,beta_cell,X_i,Y_i,m,d,m_i,lambda)
% The objective funtion and gradient of problem (A)
% In distributed setting, only transmit *U*.
% 
% para: the vector form of A.
% 
% X: sample
% Y: response 
% n: sample size of X - scalar, local client
% 
% W: local coefficient.
% m: number of client
% d: dimension of coefficient
% 
% lambda: given regularization parameter
% 
% J: The loss fuction + penalty value - the objective function value
% grad: the gradient of A in form of vector.



% 将参数转化为需要的形式
[W,~] = util_A2W(para_A,d);


J=fun_loss_H(W,X,Y,n,d)+ lambda*fun_penalty_UP(W,U_i,m,d);

grad = gradient_loss_H(W,X,Y,n,d)+lambda*gradient_UP(W,U_i,m,d);



end