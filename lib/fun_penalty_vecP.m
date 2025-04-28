function P=fun_penalty_vecP(beta_cell,m,m_i)
% P=fun_penalty_vecP(beta_cell,m,m_i)
% The loss function in original problem (average by 2 times of client number )
% 
% beta_cell - gathering of coefficients in all clients (m-cell with d×1)
% m - number of clients
% 
% scope of penalty :[0,1]
% 0 - similarity with cosine value 1
% 1 - similarity with cosine value 0, penalty.

P=0;
beta=cell2mat(beta_cell(m_i));
for i=1:m
    beta_i=cell2mat(beta_cell(i));
    P = P+ 1 - (beta*beta_i'/norm(beta)/norm(beta_i))^2; % scope of penalty :[0,1]
end

P=P/2/m;

end