function  P = fun_penalty_matP(W_cell,m,d,m_i)
% P = fun_penalty_matP(W_cell,m,d,m_i)
% Given the parameters on each client, calculate the penalty function value
% 
% W_cell: (m×1)cell with (d+1 × d+1) matrix
% m: number of client
% d: dimension of coefficient
% m_i: the penalty $P_i = \sum{P_{ij}}$
% 
% scope of penalty :[0,1]
% 0 - similarity with cosine value 1
% 1 - similarity with cosine value 0, penalty.




P=0;
[~,I2] = util_matrix_I(d);
W_mi=cell2mat(W_cell(m_i));
norm_AllW = util_norm_AllW(W_cell,m,d);
norm_mi = norm_AllW(m_i);
    for i=1:m
        W_i=cell2mat(W_cell(i));
        norm_i = norm_AllW(i);
%         P = P+ (sqrt(trace((W_mi)*(W_i)))-1)^2;
%         P = P+trace((W_mi*I2)*(W_i*I2));
        P = P+1-trace(W_mi*(I2*W_i*I2))/norm_mi/norm_i;
    end
    P=P/m/2;
end