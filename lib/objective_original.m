function [J,grad]= objective_original(beta,beta_cell,X_i,Y_i,m,n,d,m_i,lambda)
% [J,grad]= objective_original(beta,beta_cell,X_i,Y_i,m,n,d,m_i,lambda)
% beta: d*1
% X: n*d
% Y: n*1


beta_cell(m_i)=mat2cell(beta,1,d);

J=fun_loss_L(beta,X_i,Y_i,n)+lambda*fun_penalty_vecP(beta_cell,m,m_i);

    ome=generate_transfer_vector(beta_cell,m,d);
    ome_i = ome - beta'*beta/norm(beta)^2;
grad = gradient_loss_L(beta,X_i,Y_i,n)+ lambda*gradient_vP(beta,ome_i,m);

end